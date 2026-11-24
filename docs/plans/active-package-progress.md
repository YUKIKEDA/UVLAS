# IM-01・IM-03・IM-06の進捗と残作業

実装中の3パッケージは **20／43段階（約47%）が完了**した。残り23段階のうち20段階が実装・検証、3段階が終了審査。下表は現在の成果と次の完了条件を示す。

| パッケージ | 完了／全体 | 残り | 統合済みの成果 | 次に完了させる成果と判定条件 |
|---|---:|---:|---|---|
| [IM-01 構造](#im-01-構造) | 8／12（67%） | 4 | 準備済みパッチから参照量と両語の射影を受け取り、厚さ段差と非整合面の原拘束・双対を照合できる。 [準備と全厚拘束の検証](../verification/IM-01/pg-prepared-patch-trace.md) | 製品求積則で連続積分への収束を確かめ、同じ原拘束から固定ピボットの疎チャートを作る。実接続で R-02 の容量を測る。 |
| [IM-03 空力](#im-03-空力) | 6／13（46%） | 7 | 局所P1核、大域幾何、後縁所有、共有辺の誘導・全微分を実装した。各離散照会点と束縛辺の中点で非一様な既知流れを評価し、同じ試行の循環・荷重とその接線・転置に使う。 | PGの物質点写像と指定SE(3)運動を接続し、循環時間微分、自由後流、J1の状態確定と再開を検証する。 |
| [IM-06 入出力](#im-06-入出力) | 6／18（33%） | 12 | 型付きschemaとSI変換、マニフェスト直下の全参照の取得、Gmshメッシュの読取り、チェックポイント容器を実装した。容器は既存番号を保護して保存し、検証を通った最新番号のチャンクを復元する。 [チェックポイントの検証](../verification/IM-06/checkpoint-manifest.md) | Gmshのsidecarで部材、材料、支持、接続、参照面への対応を定め、役割ごとの参照を解決する。実行再開にはJ1と探索の受理済み状態が必要。 |

割合は[作業分解](next-package-work-breakdown.md)の完了段階数を全段階数で割って求める。各段階の完了は、採用仕様、後続が受け取る成果、独立した検証根拠で確かめる。段階ごとに規模が異なるため、この数字から残りの工数は見積もれない。

## IM-01 構造

[#17](https://github.com/YUKIKEDA/UVLAS/issues/17) の終了審査まで4段階。準備済みパッチの位置・勾配・測度と、非整合面の全位置拘束を使える。製品求積則の認定と、原拘束を使う疎チャート・R-02 の容量測定が残る。

| 段階 | 状態 | 到達点と残作業 |
|---|---|---|
| 厳密代数 | 完了 | 有理数・疎多項式の準備。 [#91](https://github.com/YUKIKEDA/UVLAS/pull/91) |
| グラフとセル | 完了 | 有理グラフ、half-cell、矩形の活動map。 [#103](https://github.com/YUKIKEDA/UVLAS/pull/103) |
| 全原辞書と full V | 完了 | アンカー補完、必須場と並進の厳密な写像。 [#113](https://github.com/YUKIKEDA/UVLAS/pull/113) |
| 多角形と全積分分割 | 完了 | 被覆を検査した分割と資源計画。 [分割](https://github.com/YUKIKEDA/UVLAS/pull/128)・[資源計画](https://github.com/YUKIKEDA/UVLAS/pull/149) |
| 局所独立化 | 完了 | 局所rank、基底、写像と来歴を所有する。再利用時にグラフ・領域・活動面積・辞書を照合し、B/T/S/W/Zを全活動行で厳密に再検証する。全12系列を独立参照と照合した。 [準備の検証](../verification/IM-01/pg-rank-bounded-preparation/verification.md)・[再検証と資源](../verification/IM-01/pg-rank-bounded-preparation/validation-verification.md) |
| 参照幾何と微分 | 完了 | グラフ・辞書・共通分割を所有する入口から参照Jacobian、可変法線、実法線厚さ、高階・片側微分と物理測度を返す。全活動片の局所正則性を厳密に検査する。 [実装 `09e983f`](https://github.com/YUKIKEDA/UVLAS/commit/09e983ff5d63d9ab95e0b2c75883687b64f219b1)・[検証](../verification/IM-01/pg-reference-geometry/verification.md) |
| 保証付き条件付け | 実装中 | 一定厚と区分一次厚さの平板は、内部格子n=1、2、4で切断領域の幾何内積を厳密に積分する。元の区間に対する合同誤差と必須場を検査する変換の生成入口を用意し、変厚切断n=1で認定した。放物柱面では、130列の基底に固定した1,008点の求積を適用し、幾何内積の区間と合同変換を認定した。認定結果は同じ点と重みを所有する。製品求積則の適格性と連続積分への収束が残る。 [切断領域](../verification/IM-01/pg-flat-gram-cut.md)・[区分一次厚さ](../verification/IM-01/pg-flat-gram-affine-thickness.md)・[変換生成](../verification/IM-01/pg-conditioning-generated-transform.md)・[曲面の固定求積](../verification/IM-01/pg-discrete-gram.md) |
| PreparedPatchと参照作用 | 完了 | 認定済み基底と物理参照を所有する。体積点、層界面の上下、外周側面で、参照位置、正の測度、形状値、Cartesian 勾配を取得し、成分別の誤差枠で Binary64／TwoWord に射影する。固定則の点・重みを保持し、同じ Workspace で体積位置と転置を作用させる。 [物理参照量の射影](../verification/IM-01/pg-prepared-reference-projection.md) |
| 全位置traceと拘束 | 完了 | 非整合面を共通細分し、厚さ段差、指定層と部分面、向きを検査する。原拘束と固定支持の指定値を保持し、厳密な重複行の尺度と双対仕事を照合した。時間依存の指定値と複数接続の閉路処理は固定ピボットの疎チャートで扱う。 [全厚拘束の検証](../verification/IM-01/pg-prepared-patch-trace.md) |
| 固定pivot疎chart | 残 | 元拘束と疎解法から写像と微分・転置作用を実装し、拘束と仕事を照合する |
| 早期R-02 | 残 | 実拘束グラフができた時点で、疎分解のfillと準備・再利用の容量を測る |
| 終了審査 | 残 | V-PGで曲面・切断・非整合を通し、IM-02へ引き渡す |

## IM-03 空力

[#19](https://github.com/YUKIKEDA/UVLAS/issues/19) の終了審査まで7段階。共有辺・四辺形の誘導作用を循環計算へ渡せる。既知流れの実装は独立して進められる。循環とその微分が、全荷重・後流へ進む共通の先行条件。

| 段階 | 状態 | 到達点と残作業 |
|---|---|---|
| 厳密な局所述語 | 完了 | 退化、内外、距離とfar判定。 [#90](https://github.com/YUKIKEDA/UVLAS/pull/90) |
| 所有幾何・ID・query | 完了 | 幾何snapshot、整数接続、厳密な物質位置式。 [#102](https://github.com/YUKIKEDA/UVLAS/pull/102) |
| 包絡と局所packet | 完了 | 位置・幾何包絡と局所適格性。 [#112](https://github.com/YUKIKEDA/UVLAS/pull/112) |
| P1核と全13変数の微分 | 完了 | 3基底の局所誘導核と全Jacobian。 [#129](https://github.com/YUKIKEDA/UVLAS/pull/129) |
| 共有辺・四辺形・全作用 | 完了 | 固定分割に従って全出現に散布し、核半径と強度トレースが一致する逆向きの解析辺を相殺する。値とJVPを独立面積分、VJPを双対内積で照合した。 [共有辺と四辺形の検証](../verification/IM-03/p1-shared-action-v1/verification.md) |
| 大域幾何と後縁所有 | 完了 | 大域接触・不適格の検査、alias運動学、後縁の放出所有。 [#153](https://github.com/YUKIKEDA/UVLAS/pull/153) |
| 既知流れと指定運動 | 実装中 | 一様流と解析突風の値・偏微分と、指定移流を厳密な時刻と片側で評価する。指定運動では Frame の局所軸と PCHIP の導出勾配を使い、世界の原点・姿勢・物質点を二階微分まで評価する。選んだ運動と流れを同じ過渡節点へ結び、周期延長でイベントをまたぐ区間を拒否した。取付モデル、全入力のイベント時刻表、J1の確定・再開が残る。 [既知流れ](../verification/IM-03/case-known-flow-sine-v1/verification.md)・[指定運動](../verification/IM-03/case-prescribed-rigid-motion.md)・[共通時計](../verification/IM-03/case-clocked-prescribed-inputs.md) |
| 循環と全微分 | 実装中 | 接触済みの離散照会点ごとに非一様な既知流れを評価し、実数の非貫通式を解く。同じ試行の原残差を再評価し、位置と既知流れの時刻への接線・転置を独立参照と照合した。PGの物質点写像、幾何と後流輸送も含む時間微分、接触切替え、複素更新が残る。 [非一様流での循環と微分の検証](../verification/IM-03/discrete-known-flow-circulation-v1/verification.md) |
| 全ベクトル荷重 | 実装中 | 同じ試行の循環と束縛辺の中点ごとの非一様な既知流れから、辺の力、節点力、合力、基準点モーメントを評価する。循環を介する寄与と中点からの直接寄与を接線・転置の双方で照合した。循環時間微分は独立入力として三角形の非定常荷重に使う。循環時間微分の時間法とPG構造への荷重の返却が残る。 [非一様流での中点荷重と微分の検証](../verification/IM-03/discrete-known-flow-loads-v1/verification.md) |
| 後流値・移流・試行 | 実装中 | 型付きSIケースの指定移流速度は定数・多項式・PCHIPに対応する。接触済みの一つの離散状態では、束縛渦輪、先頭・旧後流、既知流れを合わせた自由後流の粒子速度と、その接線・転置を独立参照で照合した。登録済みの後縁別名は正準点を使う。親確定状態からの自由・指定輸送、放出IDの発行、全履歴と試行の一括受理・棄却、時間積分を通じた全微分とJ1再開が残る。 [指定移流速度](../verification/IM-03/case-prescribed-wake-velocity-v1/verification.md)・[自由後流の粒子速度](../verification/IM-03/discrete-free-wake-velocity-v1/verification.md) |
| J1状態への接続 | 残 | IM-04側とstep契約を共有し、BDF2と新生台形、輸送別の初期微分、受理済み循環の時間微分、全履歴、一回限りの状態確定を検証する。強連成Newtonの実装はIM-04が担当 |
| 参照と数値予算 | 一部検証済み | [局所核の独立参照](https://github.com/YUKIKEDA/UVLAS/pull/129)を確認済み。全経路を組み立て、V-AERO・V-MULTIの観測量と数値予算を照合する。参照準備は並行可能 |
| 費用と終了審査 | 残 | 全履歴の相互作用数、資源、並列時の結果一致とR-03を確認し、パッケージを審査する |

## IM-06 入出力

[#22](https://github.com/YUKIKEDA/UVLAS/issues/22) の終了審査まで12段階。安全に取得した同じバイト列を、外部入力の形式ごとの読取りに渡す。NastranとGmsh、CLIとPythonは、それぞれの共通基盤が揃い次第並行できる。

| 段階 | 状態 | 到達点と残作業 |
|---|---|---|
| 厳格なscalar変換 | 完了 | 強いIDとIndexへの検査付き変換。 [#89](https://github.com/YUKIKEDA/UVLAS/pull/89) |
| parser資源計画 | 完了 | 固定JSON依存の例外経路を補修し、見積もりと受入条件を実装。 [依存補修](https://github.com/YUKIKEDA/UVLAS/pull/101)・[資源計画](https://github.com/YUKIKEDA/UVLAS/pull/104) |
| 所有文書 | 完了 | 元バイト列、digest、型、十進字句を保持。 [#104](https://github.com/YUKIKEDA/UVLAS/pull/104) |
| 全schema | 完了 | 全caseの型、入れ子、参照と相互条件の検証。 [#130](https://github.com/YUKIKEDA/UVLAS/pull/130) |
| SI・精度変換 | 完了 | 元の厳密値、SI倍率、成分別の精度管理。 [#154](https://github.com/YUKIKEDA/UVLAS/pull/154) |
| 参照snapshot | 実装中 | 宣言収集、パスの正規化、[単一ハンドルの照合](https://github.com/YUKIKEDA/UVLAS/pull/177)は統合済み。Windowsの根相対adapterは [#184](https://github.com/YUKIKEDA/UVLAS/pull/184)、Linuxの根相対adapterは [#195](https://github.com/YUKIKEDA/UVLAS/pull/195) で統合済み。[#202](https://github.com/YUKIKEDA/UVLAS/pull/202) で、マニフェストの元バイト列とdigestを単一のファイルハンドルから取得する部品を統合した。[#219](https://github.com/YUKIKEDA/UVLAS/pull/219) で、Fileは元マニフェストのハンドルと実体ID、Memoryは呼出し側が開いた根を検査済みケースとともに保持するようにした。[#229](https://github.com/YUKIKEDA/UVLAS/pull/229) で、参照ファイルを開く前に、宣言順と位置を保持して参照先を正規化・索引化する部品を統合した。[#236](https://github.com/YUKIKEDA/UVLAS/pull/236) で、指定された参照先1件の実体を内容の読み取り前に照合し、同じハンドルから取得した元バイト列を保持する部品を統合した。[#245](https://github.com/YUKIKEDA/UVLAS/pull/245) で、宣言順に依存をたどり、循環と共有経路の深さ上限を検査する内部処理を統合した。[#252](https://github.com/YUKIKEDA/UVLAS/pull/252) で、先頭の根宣言から取得したバイト列とダイジェストを照合し、依存走査へ渡す内部処理を統合した。マニフェスト直下の全8役割を同じ親の参照表で保持し、取得済みのバイト列からGmshメッシュを読める。宣言順に外部ファイルを取得し、元バイト列、依存辺、宣言位置、依存先を先に並べた完了順を所有する。形式ごとの子依存の抽出と、未認定の役割を含む入力の認定が残る。 [参照グラフの検証](../verification/IM-06/reference-graph.md) |
| 励振と既知入力 | 実装中 | 型付きSIケースから突風なしの一様流を準備し、定数・多項式・PCHIPの値と時間微分を評価する。指定した成分別の誤差予算でBinary64／TwoWordに投影した値を認定できる。参照snapshot、他の励振、J1への接続が残る |
| parametric入力とprepare | 残 | 型付きcaseを共通のRawModelへ渡す。物理準備には対応するIM-01／02／03の成果が必要 |
| Nastran | 実装中 | 取得済みdeckから数値欄と継続を読み、GRIDとCORD2R/C/Sの参照を検査する。CTRIA3とCQUAD4は節点順、PSHELL／PCOMPの参照、材料軸、オフセット、厚さ指定を元の欄とともに保持する。MAT1・PSHELLとMAT8・PCOMPの材料参照、積層の宣言順と継承元を検査し、SPC・SPC1・SPCADDの支持集合とCONM2の重心・慣性欄を保持する。CBEAM の節点・方向・端部指定と、RBE2 の従属節点・自由度・熱膨張欄を保持する。同じデッキの CBEAM から PBEAM の断面・軸方向ステーションと MAT1 をたどれる。CLI は照合済みデッキを走査し、採用カードだけの単一ファイルであることを検査する。単位変換、梁断面の補間と K6／質量への対応、剛体拘束の組立て、sidecar が残る。 [座標・シェル](../verification/IM-06/nastran-shell-cards.md)・[材料・積層](../verification/IM-06/nastran-laminate-cards.md)・[支持集合](../verification/IM-06/nastran-support-cards.md)・[集中質量](../verification/IM-06/nastran-lumped-mass-cards.md)・[はり](../verification/IM-06/nastran-beam-cards.md)・[剛体結合](../verification/IM-06/nastran-rigid-link-cards.md)・[はり断面](../verification/IM-06/nastran-pbeam-cards.md)・[デッキの参照検査](../verification/IM-06/nastran-bulk-leaf-verification.md) |
| Gmsh | 一部実装 | 検査済みのスナップショットからMSH 4.1 ASCIIを読み、疎な節点タグ、所属、物理名、向きと元の座標字句を保持する。CLI は同じバイト列からメッシュ本体に子参照がないことを確かめる。明示した面集合を検査し、共有辺の向きと頂点のつながりを照合して面の連結成分と自由境界を返す。sidecar の型、部材への割当、幾何の自己交差と物理モデルへの対応が残る。 [入力契約](../specs/foundation/gmsh-mesh-input-contract.md)・[面位相](../verification/IM-06/gmsh-surface-topology.md) |
| 共通参照への変換 | 残 | 外部面を共通の曲面・厚さ・接続へ変換し、質量・仕事などの誤差を評価する。対応するIM-01／02の成果が必要 |
| 結果と診断 | 実装中 | 状態と来歴のJSON、確定観測のCSV、適格な可視化標本のVTU・PVDを内部で生成できる。VTU・PVDは層界面の両側と再始動区間を保つ。型検査・SI変換済みケースから、元入力、保持値、誤差半径、現時点で解決した既定値を出力し、その指紋を来歴に結び付けられる。物理回復から可視化標本を作り、物理準備で確定する設定とともに実行結果の保存へつなぐ工程が残る。 [解決済みケース](../verification/IM-06/resolved-case-document.md)・[観測CSV](../verification/IM-06/observation-csv.md)・[VTU・PVD](../verification/IM-06/visualization-series.md) |
| checkpoint容器 | 完了 | 型付きチャンクと指紋からマニフェストを作り、既存番号を保持したまま確定する。中断した保存を除外し、型、長さ、SHA-256、指紋が一致する最新番号を復元する。 [マニフェストと確定ディレクトリの検証](../verification/IM-06/checkpoint-manifest.md) |
| 実状態とresume | 残 | 全BDF2ステンシル、受理済みの微分・荷重・後流履歴、探索状態などを保存する。IM-04／05の実状態を使い、同版・同環境での再開結果のビット一致と非互換の拒否を確認する |
| CLI | 実装中 | `validate` は全スキーマと対応する外部参照を検査し、未対応の役割を拒否する。直下の外部ファイルは SHA-256 を照合し、MSH 4.1 ASCII は節とタグも検査する。`inspect` は保存済み状態と来歴を照合し、来歴に指紋がある場合は解決済みケースも確認する。Python と検査処理を共有する。形式別の子依存の検査、prepare／run／resume、実行終了時の状態記録が残る。 [開発用CLI](../build.md#開発用cliの終了規約) |
| Python | 実装中 | ファイル・辞書の検査と結果確認に CLI と同じ C++ 処理を使い、検査済み入力と参照を所有する。物理準備・実行・再開、受理済み解析結果への接続が残る。 [入口の検証](../verification/IM-06/shared-case-frontends.md) |
| ランタイムと新規環境 | 基盤あり | Windows と Linux で非公開ランタイムを移設し、開発ツールの PATH を外した CLI と新規 Python 環境でケース検査と所有権を確認した。配布用 CLI／wheel の対応環境と、導入後の実行・再開の受入が残る。 [移設後の検証](../verification/IM-06/shared-case-frontends.md) |
| 終了審査 | 残 | V-USE・V-STATEで3入力経路とCLI／Pythonの物理結果、両OSの導入・中断・再開を照合する |

schemaのLinux sanitizerの2試験は [#143](https://github.com/YUKIKEDA/UVLAS/pull/143) の修正を含むmainで再計測し、49.51秒と53.22秒で通過した。Linux sanitizerだけ300秒、Linux Debug・ReleaseとWindows Debug・Releaseは60秒の登録を確認した。 [#142](https://github.com/YUKIKEDA/UVLAS/issues/142) 数値表では同じ段階の子Issueを重ねて加えない。

## 初回実務適用までの位置

IM-00の共通基盤は終了審査を通過した。IM-01の成果をIM-02の材料・梁・接続・構造全体へ渡し、IM-02とIM-03の成果をIM-04の強連成へつなぐ。その後にIM-05の安定性・探索、IM-07の実務規模の数値・性能検証、IM-08の実験・実務受入がある。IM-06は並行しながら、実行と再開の段階で揃った物理状態を接続する。 [全体の依存とゲート](implementation-start-review.md#3-後続パッケージと依存)
