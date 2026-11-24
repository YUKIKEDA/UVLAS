# IM-01・IM-03・IM-06の進捗と残作業

実装中の3パッケージは **18／43段階（約42%）が完了**した。残り25段階のうち22段階が実装・検証、3段階が終了審査。下表は現在の成果と次の完了条件を示す。

| パッケージ | 完了／全体 | 残り | 統合済みの成果 | 次に完了させる成果と判定条件 |
|---|---:|---:|---|---|
| [IM-01 構造](#im-01-構造) | 6／12（50%） | 6 | 参照幾何に加え、辞書と分割から局所rankを準備し、基底と全座標を現在の入力で厳密に再検証できる。全12系列で、元の結果を保持したまま同じ48 GiB枠で照合した。 [局所rank再検証](../verification/IM-01/pg-rank-bounded-preparation/validation-verification.md) | 局所基底と参照幾何から物理内積を組み立て、保証付き条件付けの誤差上限と必須場の保持を検査する。 |
| [IM-03 空力](#im-03-空力) | 6／13（46%） | 7 | 局所P1核、大域幾何、後縁所有、共有辺の誘導・全微分を実装した。定数・多項式・PCHIPの一様流は厳密なSI入力から速度と微分を返し、復元時に入力と係数を照合する。 [一様流の検証](../verification/IM-03/exact-uniform-flow-v1/verification.md) | ケースからの受渡しと必要な精度投影を閉じ、突風と指定運動を同じ時刻・座標規則で照合する。 |
| [IM-06 入出力](#im-06-入出力) | 6／18（33%） | 12 | 型付きschemaとSI変換、マニフェスト直下の全参照の取得、Gmshメッシュの読取り、チェックポイント容器を実装した。容器は既存番号を保護して保存し、検証を通った最新番号のチャンクを復元する。 [チェックポイントの検証](../verification/IM-06/checkpoint-manifest.md) | Gmshのsidecarで部材、材料、支持、接続、参照面への対応を定め、役割ごとの参照を解決する。実行再開にはJ1と探索の受理済み状態が必要。 |

割合は[作業分解](next-package-work-breakdown.md)の完了段階数を全段階数で割って求める。各段階の完了は、採用仕様、後続が受け取る成果、独立した検証根拠で確かめる。段階ごとに規模が異なるため、この数字から残りの工数は見積もれない。

## IM-01 構造

[#17](https://github.com/YUKIKEDA/UVLAS/issues/17) の終了審査まで6段階。参照幾何から位置・微分・測度を受け取り、局所基底とその再構成を現在の入力で確かめられる。次の段階は保証付き条件付けである。

| 段階 | 状態 | 到達点と残作業 |
|---|---|---|
| 厳密代数 | 完了 | 有理数・疎多項式の準備。 [#91](https://github.com/YUKIKEDA/UVLAS/pull/91) |
| グラフとセル | 完了 | 有理グラフ、half-cell、矩形の活動map。 [#103](https://github.com/YUKIKEDA/UVLAS/pull/103) |
| 全原辞書と full V | 完了 | アンカー補完、必須場と並進の厳密な写像。 [#113](https://github.com/YUKIKEDA/UVLAS/pull/113) |
| 多角形と全積分分割 | 完了 | 被覆を検査した分割と資源計画。 [分割](https://github.com/YUKIKEDA/UVLAS/pull/128)・[資源計画](https://github.com/YUKIKEDA/UVLAS/pull/149) |
| 局所独立化 | 完了 | 局所rank、基底、写像と来歴を所有する。再利用時にグラフ・領域・活動面積・辞書を照合し、B/T/S/W/Zを全活動行で厳密に再検証する。全12系列を独立参照と照合した。 [準備の検証](../verification/IM-01/pg-rank-bounded-preparation/verification.md)・[再検証と資源](../verification/IM-01/pg-rank-bounded-preparation/validation-verification.md) |
| 参照幾何と微分 | 完了 | グラフ・辞書・共通分割を所有する入口から参照Jacobian、可変法線、実法線厚さ、高階・片側微分と物理測度を返す。全活動片の局所正則性を厳密に検査する。 [実装 `09e983f`](https://github.com/YUKIKEDA/UVLAS/commit/09e983ff5d63d9ab95e0b2c75883687b64f219b1)・[検証](../verification/IM-01/pg-reference-geometry/verification.md) |
| 保証付き条件付け | 実装中 | 一定高さ・一定厚さ・全被覆の平板は、内部格子n=1、2、4の幾何内積を解析積分し、外向き区間を返す。n=1で合同変換の認定と層・面側の基底評価を照合した。曲面・変厚・切断領域への保証付き積分と変換候補の製品生成が残る。 [平板内積の検証](../verification/IM-01/pg-flat-gram/verification.md) |
| PreparedPatchと参照作用 | 実装中 | 認定済み基底の全列について、位置、Jacobian、測度、値と物理勾配を体積・層界面の片側・側面の点列から一括で取得できる。体積質量の準備に渡す基底値と正の体積Jacobianは、Binary64またはTwoWordの格納誤差を認定して保持できる。次は層密度と求積を結ぶ質量作用、材料・荷重の参照作用と反復用Workspaceへつなぐ。 [点列評価](../verification/IM-01/pg-conditioned-reference-batch.md)・[体積質量に渡す参照値](../verification/IM-01/pg-volume-mass-reference.md) |
| 全位置traceと拘束 | 残 | 辞書・分割・基底を使い、非整合界面の全厚位置、元拘束、双対量を保持する |
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
| 既知流れと指定運動 | 実装中 | [PCHIP係数](https://github.com/YUKIKEDA/UVLAS/pull/179)は統合済み。指定区間の値・時間微分も [#187](https://github.com/YUKIKEDA/UVLAS/pull/187) で統合済み。閉定義域内の節点検索と左右区間候補の取得は [#197](https://github.com/YUKIKEDA/UVLAS/pull/197) で統合済み。[#205](https://github.com/YUKIKEDA/UVLAS/pull/205) で、範囲外のRejectとHold、およびHold端点の内外候補の取得を統合した。[#211](https://github.com/YUKIKEDA/UVLAS/pull/211) で、宣言した周期と節点幅の一致を検査し、負時刻を含む周期時刻の厳密な正規化を統合した。[#218](https://github.com/YUKIKEDA/UVLAS/pull/218) で、要求区間、Hold端点のゼロ傾き、周期端の値と時間微分の一致を厳密に検査する部品を統合した。[#228](https://github.com/YUKIKEDA/UVLAS/pull/228) で、微分を評価する側の選択と、Automaticで必要になる導関数の厳密な一致判定を統合した。[#237](https://github.com/YUKIKEDA/UVLAS/pull/237) で、準備入力を所有し、発行元・世代・実内容を照合する部品を統合した。[#244](https://github.com/YUKIKEDA/UVLAS/pull/244) で、係数準備と境界検査の確保量・処理量を係数の確保前に見積もる部品を統合した。[#251](https://github.com/YUKIKEDA/UVLAS/pull/251) で、確定済みの準備入力を読み取り専用で共有する部品を統合した。[#263](https://github.com/YUKIKEDA/UVLAS/pull/263) で、曲線の一括受入・評価と実内容を照合する復元までを検証した。既知流れ・指定運動への接続、空間・時間微分と再開時の一致の検証が残る |
| 循環と全微分 | 実装中 | 接触済みの離散照会点と明示した法線・面速度から、先頭後流を含む実数の循環方程式を組み立てて解き、非貫通の原式で循環を検査する。固定幾何の循環作用は接線と転置も求解し、投影前の作用で保持した解の残差を検査する。PG物質写像に基づくコロケーション、幾何・時刻を含む全微分、複素循環更新が残る。 [離散循環の検証](../verification/IM-03/discrete-circulation-v1/verification.md)・[接線・転置の検証](../verification/IM-03/discrete-circulation-linearization-v1/verification.md) |
| 全ベクトル荷重 | 残 | 循環・共有核・幾何から全荷重と微分を実装する。構造への双対写像にはIM-01／02の成果が必要 |
| 後流値・移流・試行 | 残 | 同じ循環・荷重状態から自由／指定輸送、全履歴、試行の親不変性と全微分を検証する |
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
| 参照snapshot | 実装中 | 宣言収集、パスの正規化、[単一ハンドルの照合](https://github.com/YUKIKEDA/UVLAS/pull/177)は統合済み。Windowsの根相対adapterは [#184](https://github.com/YUKIKEDA/UVLAS/pull/184)、Linuxの根相対adapterは [#195](https://github.com/YUKIKEDA/UVLAS/pull/195) で統合済み。[#202](https://github.com/YUKIKEDA/UVLAS/pull/202) で、マニフェストの元バイト列とdigestを単一のファイルハンドルから取得する部品を統合した。[#219](https://github.com/YUKIKEDA/UVLAS/pull/219) で、Fileは元マニフェストのハンドルと実体ID、Memoryは呼出し側が開いた根を検査済みケースとともに保持するようにした。[#229](https://github.com/YUKIKEDA/UVLAS/pull/229) で、参照ファイルを開く前に、宣言順と位置を保持して参照先を正規化・索引化する部品を統合した。[#236](https://github.com/YUKIKEDA/UVLAS/pull/236) で、指定された参照先1件の実体を内容の読み取り前に照合し、同じハンドルから取得した元バイト列を保持する部品を統合した。[#245](https://github.com/YUKIKEDA/UVLAS/pull/245) で、宣言順に依存をたどり、循環と共有経路の深さ上限を検査する内部処理を統合した。[#252](https://github.com/YUKIKEDA/UVLAS/pull/252) で、先頭の根宣言から取得したバイト列とダイジェストを照合し、依存走査へ渡す内部処理を統合した。マニフェスト直下の全8役割を同じ親の参照表で保持し、取得済みのバイト列からGmshメッシュを読める。各形式の子依存を抽出し、全依存グラフを閉じる工程が残る。 [直下参照の検証](../verification/IM-06/declared-reference-inventory.md) |
| 励振と既知入力 | 実装中 | 型付きSIケースから突風なしの一様流を準備し、定数・多項式・PCHIPの値と時間微分を評価する。指定した成分別の誤差予算でBinary64／TwoWordに投影した値を認定できる。参照snapshot、他の励振、J1への接続が残る |
| parametric入力とprepare | 残 | 型付きcaseを共通のRawModelへ渡す。物理準備には対応するIM-01／02／03の成果が必要 |
| Nastran | 残 | schema・SI・snapshotを入力に、採用カード、座標系、sidecarを解析する |
| Gmsh | 一部実装 | 検査済みのスナップショットからMSH 4.1 ASCIIを読み、疎な節点タグ、所属、物理名、向きと元の座標字句を保持する。CLIは同じバイト列からメッシュ本体に子参照がないことを確かめる。sidecarの型と物理モデルへの対応が残る。 [入力契約](../specs/foundation/gmsh-mesh-input-contract.md) |
| 共通参照への変換 | 残 | 外部面を共通の曲面・厚さ・接続へ変換し、質量・仕事などの誤差を評価する。対応するIM-01／02の成果が必要 |
| 結果と診断 | 実装中 | 状態と来歴のJSON、確定観測のCSV、適格な可視化標本のVTU・PVDを内部で生成できる。VTU・PVDは層界面の両側と再始動区間を保つ。型検査・SI変換済みケースから、元入力、保持値、誤差半径、現時点で解決した既定値を出力し、その指紋を来歴に結び付けられる。物理回復から可視化標本を作り、物理準備で確定する設定とともに実行結果の保存へつなぐ工程が残る。 [解決済みケース](../verification/IM-06/resolved-case-document.md)・[観測CSV](../verification/IM-06/observation-csv.md)・[VTU・PVD](../verification/IM-06/visualization-series.md) |
| checkpoint容器 | 完了 | 型付きチャンクと指紋からマニフェストを作り、既存番号を保持したまま確定する。中断した保存を除外し、型、長さ、SHA-256、指紋が一致する最新番号を復元する。 [マニフェストと確定ディレクトリの検証](../verification/IM-06/checkpoint-manifest.md) |
| 実状態とresume | 残 | 全BDF2ステンシル、受理済みの微分・荷重・後流履歴、探索状態などを保存する。IM-04／05の実状態を使い、同版・同環境での再開結果のビット一致と非互換の拒否を確認する |
| CLI | 実装中 | `validate` は全スキーマと参照を検査する。外部ファイルはマニフェスト直下の取得と SHA-256 の照合を行い、MSH 4.1 ASCII は節とタグも検査する。未対応の宣言は位置を診断する。形式別の子依存の検査、prepare／run／resume／inspect、実行終了時の状態記録が残る。 [開発用CLI](../build.md#開発用cliの終了規約)・[検証](../verification/IM-06/im06-cli-case-validation.json) |
| Python | 残 | 同じC++処理にPythonの辞書・ファイル・実行・再開と所有権を接続する |
| ランタイムと新規環境 | 基盤あり | CLIとbindingの基盤は [IM-00](im00-package-acceptance.md) にある。実際に配布するCLI／wheelの依存・ABIと、開発環境に依存しない導入・実行を確認する |
| 終了審査 | 残 | V-USE・V-STATEで3入力経路とCLI／Pythonの物理結果、両OSの導入・中断・再開を照合する |

schemaのLinux sanitizerの2試験は [#143](https://github.com/YUKIKEDA/UVLAS/pull/143) の修正を含むmainで再計測し、49.51秒と53.22秒で通過した。Linux sanitizerだけ300秒、Linux Debug・ReleaseとWindows Debug・Releaseは60秒の登録を確認した。 [#142](https://github.com/YUKIKEDA/UVLAS/issues/142) 数値表では同じ段階の子Issueを重ねて加えない。

## 初回実務適用までの位置

IM-00の共通基盤は終了審査を通過した。IM-01の成果をIM-02の材料・梁・接続・構造全体へ渡し、IM-02とIM-03の成果をIM-04の強連成へつなぐ。その後にIM-05の安定性・探索、IM-07の実務規模の数値・性能検証、IM-08の実験・実務受入がある。IM-06は並行しながら、実行と再開の段階で揃った物理状態を接続する。 [全体の依存とゲート](implementation-start-review.md#3-後続パッケージと依存)
