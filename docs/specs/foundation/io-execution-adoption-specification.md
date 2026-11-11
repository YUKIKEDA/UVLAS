# 入出力・実行・再開の実装仕様

2026-11-06。版1.0。作業4の採用仕様。[数値基盤](numerical-foundation-adoption-specification.md)とJ1の状態契約を継承する。これは未実装の製品入力形式を説明する設計であり、利用可能なCLI/APIの説明書ではない。

## IO4-01 ケース形式と実行

入口はUTF-8 JSON `case.json`、schema `uvlas.case/1`。コメント・重複キー・未知キー・NaN/Infinity・型の暗黙変換を拒否する。独自parserを作らずnlohmann/json 3.12.0を版固定で使用し、構文解析後に全schema制約を検査する。IDは10進文字列のuint64、配列indexは別の0始まりint64。パスはmanifestからの相対参照を基本とし、元ファイルと解決後の内容hashを保存する。

必須キーは`schema, units, model, flow, analysis, observations, numerics, resources, output`。任意キーは`initial, restart, provenance`。Modelは`geometry, materials, shells, beams, joints, supports, masses, aero_surfaces`を持ち、空集合も明示する。未解決の既定値は計算前の`resolved-case.json`へすべて展開する。材料・減衰・流体密度・支持剛性に都合のよい暗黙既定値を置かない。

| ブロック | 主要な必須項目 |
|---|---|
| units | length, mass, time, angle。内部SIへ変換。温度等を含まない入力へ追加単位を推測しない |
| geometry | source種別、元ファイルhash、チャート原点／基底／長さ、面内多角形、次数付き参照係数、厚さ・層境界、滑らかな共有面と折れの識別 |
| materials | 3D弾性の全必要定数または認定済み等価モデル、密度、材料軸、層順、初期ひずみ。正定値性・適用条件 |
| flow | ρ、既知一様流／時空間入力、採用gustの時刻・位置基準と補間、各揚力面の物理ID・後縁 |
| analysis | static／transient／flutter、対象branch、開始・終了・h・後流政策、探索区間・キュー政策、物理／算法初期条件 |
| observations | 点・断面・層・両側・厚さ位置・基準点／軸、応力測度、比較窓、要求精度と非ゼロ尺度 |
| numerics | patch分割、求積、コア、後流Rとh、各求解上限・尺度、精度政策、固有抽出範囲、誤差配分 |
| resources | memory profile、スレッド数、最大メモリ・時間・反復、checkpoint用一時容量 |
| output | 保存先、観測頻度、幾何・層別出力、checkpoint間隔、再開政策 |

CLI設計は `uvlas validate case.json`、`uvlas prepare case.json`、`uvlas run case.json`、`uvlas resume checkpoint-dir`、`uvlas inspect result-dir`。validateは構文・参照、prepareは参照基底／幾何／rank／容量まで検査する。runは同じprepare成果物をhash照合して再利用できる。入力・準備・通常実行・再開を別の物理モデルへ分岐させない。

終了codeは0=正常完了、2=入力／未対応、3=不適格幾何／拘束、4=数値未収束／精度未解決、5=資源上限、6=保存／互換性、130=利用者中断。終了codeだけでなく`status.json`にFailureと最後の確定状態を記録する。未収束結果に正常完了フラグを付けない。

## IO4-02 Pythonと配布

CPython 3.12、NumPy 2.x、pybind11 3.1.0を最初の検証対象とする。C++23 coreと同じModel／処理を呼び、Pythonに別ソルバーを実装しない。schema辞書による構築、ファイル読込、実行、再開、所有する結果の取得を提供する。NumPy配列の寿命と可変性は数値基盤の規約を守る。

Python入力は呼出し中に不変なコピーを基準とし、計算開始後のPython側書換えがModelへ伝播しない。利用者callbackは数値反復内で直接呼ばず、既知入力を計算前に認定済み補間へ変換する。キャンセル・進捗は確定境界で扱う。単一interpreterの通常GIL版を対象とし、free-threadedや複数interpreterの安全性を未検証で主張しない。

ソース、Windows x64／Linux x86_64 CLI、CPython 3.12 wheelを最初の配布対象として構築する。Linuxの最小glibc・manylinux tagは実装開始直後の配布チケットで実ランタイムに合わせて固定し、WSL上の成功だけでmanylinux互換を名乗らない。MSVC／libstdc++／Pythonランタイムと第三者通知を収集し、開発機のPATHやコンパイラに依存しない新規環境で導入する。初回公開前にCLI／wheel／ソースの同一ケース結果・失敗・再開を照合する。

## IO4-03 Nastran系の限定読込み

MSC Nastran 2022.1 Quick Reference GuideのBulk Dataを基準の意味仕様とし、元deckの生成元・dialectを指定する。small／large field、free field、継続行、D指数と省略E指数を処理し、重複ID・未解決参照・循環座標系を拒否する。単位はdeckから推測せずsidecarで必須指定する。

| カード群 | 受け入れる意味と変換 |
|---|---|
| GRID、CORD2R/C/S | 参照位置・成分座標系を基本Cartesianへ解決。位置とvectorの円筒／球変換を区別 |
| CQUAD4、CTRIA3 | 外部シェル面の位相・向き・材料軸／offsetを取得。PGの要素番号と同一視しない |
| MAT1、PSHELL | 同じ等方材の物理厚さ／密度が定義できる場合の3Dモデルへ変換。別々の膜／曲げ／せん断材、非標準曲げ倍率、MID4等は同じ3Dモデルで再現可能と認定されなければ拒否 |
| MAT8、PCOMP | 層・面内材料軸・厚さを読込み。不足する3D定数をsidecarで必須補完し、任意のE3／ν13等を黙って推定しない |
| CBEAM、PBEAM | 断面・基底・偏心、軸方向分布から採用梁のK6／質量へ変換。せん断中心・重心を分離。拘束反り、端release、未対応なstress recoveryは個別拒否 |
| CONM2 | 質量、指定基準点・基底の偏心と慣性を変換。慣性成分の符号・並びを基準仕様で解決し平行軸移動を一度だけ適用 |
| RBE2 | 指定物理点・部材frameの剛体関係として構成。外部の回転DOF番号をANCF係数へ直接割当てない |
| SPC、SPC1、SPCADD | 外部の並進／回転支持の物理的意味を接続・frameへ変換。曖昧な回転支持は拒否 |

deckを「構造抽出」として読む場合、SOL／METHOD等の解析制御、CAERO／SPLINE等の空力設定は**明示した抽出対象外allowlist**へ記録し、UVLAS側の独立なanalysis／aero設定を必須とする。未知の剛性・質量・拘束カードをallowlistへ一般化して読み飛ばさない。無視した行番号・カード・理由をPreparationReportへ残す。

参照： [MSC公式2022.1 QRG](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)。PSHELLが膜・曲げ・せん断の別特性を持つことは [Autodesk公式PSHELL](https://help.autodesk.com/cloudhelp/2024/ENU/NSTRN-Reference/files/GUID-075E0B52-7716-4BCB-9A3D-CE81C9290123.htm)でも確認できる。dialect間の同名カードの追加欄を同一と仮定しない。

## IO4-04 Gmshと外部面からの参照構築

GmshのMSH **4.1 ASCII**を最初の対応形式とする。Nodes／Elements／Entities／PhysicalNamesを解決し、非連続node tagを保持する。幾何点・一次線・一次三角形・一次四角形を対象とし、partition／ghost／periodic／高次要素／binaryは明示的に未対応とする。physical groupとsidecarで部材・材料・支持・接続を関連付ける。ファイル形式の根拠は [Gmsh公式MSH仕様](https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format)。Gmsh本体やそのAPIを実行時依存にしない。

外部面は**面の幾何と位相**であり、その離散DOFをPGへ直接移すものではない。まず連結した部材を識別し、重複・非多様体・向き・自己交差・閉じていない接続を検査する。三角形は採用PGの面内多角形領域として扱うことができる。矩形へ面積を埋めたり、別の三角形要素を無断採用しない。痩せた切断領域でrank／質量の適格性を満たせない場合は再分割して再認定する。

平面部材は同一平面の参照グラフへ厳密に写す。滑らかな曲面は明示された参照多項式／曲面sidecarを優先する。外部点だけから作る場合は、次数4以下の多項式patchと共有境界上の位置／法線／全厚の一致を拘束した最小二乗近似を構築し、元点・辺・法線の誤差と分割感度を報告する。有限の点列から未観測の曲率を正確と認定しない。滑らかな部材と折れを利用者の部材区分に従い区別する。

境界の定義が不足し適格な共通グラフを作れなければ、参照幾何sidecar不足として拒否する。元モデルから勝手に剛性を変えるsmooth処理をしない。外部モデルのoffset・層順・法線を反転した場合は、物理材料・質量・接続対応も変換する。

変換誤差の受入は、幾何の距離だけで決めない。少なくとも質量・重心・慣性、剛体仕事、独立荷重の変位／反力、対象固有mode／周波数を変換前後で比較し、各観測の全数値予算の中へ幾何・変換誤差を計上する。外部solverの結果がない場合は同じ幾何を解析的に定義した小ケースと独立FEMで検証する。元deckとPGが同じ離散式になるとは主張しない。

## IO4-05 結果

`status.json`、`resolved-case.json`、`provenance.json`、観測`csv`、可視化`vtu/pvd`を基準とする。CSVの列には物理量・単位・基底・参照点・層／面側を結び付け、巨大uint64 IDをfloatへ落とさない。時系列は物理時刻・確定番号・segmentを記録する。

VTUの出力格子は可視化用の標本で、計算のPG自由度や独立な三角形シェルではない。変形前後の位置、法線、物理厚さ、層別の測度指定応力・ひずみを保存し、層界面の両側値を一つに平均しない。断面力は基準点／基底、反力は支持の物理共役量を含む。未認定応力や未対応な梁局所回復をゼロとして出さない。

## IO4-06 Checkpoint

形式はディレクトリ：UTF-8 `manifest.json`と型付き無圧縮binary chunks。schema `uvlas.checkpoint/1`、byte order=little、dtype=f64/i64/u64、shape、byte数、SHA256、モデル・basis・材料・拘束・算法・精度・版のfingerprintを必須とする。各型を手動エンコードし、C++構造体のpaddingやpointerをdumpしない。TwoWordはhi／lo別chunkで正確に保存する。

現在・前状態、BDF2に必要な実ステンシル、現在受理済みD_hv／D_hΓ、全荷重、物理IDとTE所有、全後流位置／μ／履歴、輸送政策、時刻・h・segment、チャート／姿勢原点、既知入力の補間係数、探索枝・mode対応・完了点・未処理queue・不確かさ・資源残を保存する。分解cacheは保存しなくてよいが、復元時の再構築費用を計上する。

保存先と同じfilesystemの新規一時ディレクトリにchunkを書き、flush・close・hash照合後、manifestを最後に完成させて新しい確定番号のディレクトリへrenameする。既存の確定checkpointを上書きしない。利用可能な最新番号は完成manifestの検査で判定する。未完成一時保存から正常復元しない。保存中に親状態を変更しない。

読込時はschema、型・寸法・byte数・hash・モデル／算法／精度、時刻とID、J1の親子・履歴を検査してから状態を構成する。path traversalや重複chunk、過大サイズも拒否する。同一実行版・数値環境では連続／再開bitwise一致を要求し、異OS／異compilerは許容差内の同一物理・数値条件として別に試験する。互換性のない版の移行を無検証に自動実行しない。

## IO4-07 実装で閉じる検査

各対応カードの手計算・変換fixture、GmshとNastranから同じ物理モデルを構成する比較、曲面近似・triangle活動領域・層順／offset、CLI／Python一致、配列寿命、例外、破損／欠落／非互換／中断保存、同一版の再開を必須とする。S10のカード棚卸しだけをimport合格としない。配布物の導入と実務全体の受入は本体実装後のgateであり、今回のmetadata／binding試作の成立範囲へ加えない。
