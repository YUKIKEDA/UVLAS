# IM-01の終了審査とIM-02へ渡す状態

審査対象は統合ローカル版 `2f70cff11ca9b7b55a6e1e1919a5cef9e03ecac0`。2026-11-25の参照時点では、originのHEADとmainは `60a356c7031e9ba5895527c4349209f92902a825`。公開PRの [#351](https://github.com/YUKIKEDA/UVLAS/pull/351) と [#361](https://github.com/YUKIKEDA/UVLAS/pull/361) はどちらも未統合で、両PRの公開headは審査対象の祖先にある。局所独立化の子 [#133](https://github.com/YUKIKEDA/UVLAS/issues/133) と親 [#17](https://github.com/YUKIKEDA/UVLAS/issues/17) は開いているため、IM-01の終了判定は保留とする。

## 採用分担表との照合

| IM-01の完了条件 | 実行証拠と判定範囲 |
| --- | --- |
| 各凍結入力の参照幾何・物理測度、格子と厚さ | 凍結R4の28配置と後退テーパR5の4配置で、全体積片、層側と外周の参照位置、法線、厚さと測度を独立式と照合した。観測した格子はR4の1/2/4、R5の1/2。層界面の両側は同じ形状の二層対照を用いた。凍結Caseは34件・36 Shellで、形状と厚さを型付きSI入口へ渡した。[R4幾何](pg-r4-domain-geometry.md)・[R5幾何](pg-r5-domain-geometry.md)・[Case製品参照](vpg-case-product-reference.md) |
| 局所基底、必須場、条件付け、資源 | #133の6条件は局所開始審査で追跡済み。独立Fractionの9系列と6変異、旧領域分割の48 GiB事前拒否、採用方式の全再構成と実体照合を確認した。凍結32配置の有限則は全385,968必須場成分を厳密比較し、合同誤差上界は最大 5.51883335e−38、最大累積課金11,428,331,040 B、最大常駐量6,522,732,544 B。格子系列は両OSで分担し、mixed quarticの同一配置を両OSで照合した。[局所独立化の開始審査](pg-local-rank-start-review.md)・[採用した資源契約](pg-rank-bounded-preparation/validation-verification.md)・[有限則32配置](pg-vpg-finite-conditioning/verification.md) |
| 固定・剛体・時変の全厚拘束、チャート、仕事、反力、実接続とfill | R4二曲面の整合・非整合を実Caseから準備した。非整合は元288行、rank 171、元CSR70,428項、L/U11,106/34,860で両OSが一致した。曲面の指定支持は同じ所有チャートで位置・速度・加速度と片側を照合した。大回転の剛体支持は別の実全厚で元行・世代・反力仕事を照合した。R5の全6面と体積の独立荷重仕事、R5の480点での非アフィン2場、実PGの2／3／4パッチ系列の準備と再利用（4パッチは閉路）を記録した。[R4接続](pg-case-full-chart/verification.md)・[Caseの時変支持](pg-case-moving-chart/verification.md)・[大回転支持](pg-rigid-support-motion.md)・[独立外力](pg-r5-swept-taper-reference.md)・[非アフィン場](pg-case-reference-handoff/physical-reference.md)・[疎接続系列](pg-sparse-chart-series.md) |
| 両OS、来歴、独立な正負対照と引継ぎ | 凍結Case展開のハッシュを両OSで照合した。形状・原行の独立積和、誤向き、反射枠、同種ID重複、世代不一致、片側、容量、取消の対照がある。R5の支持付きCaseでは、入力を破棄した後に同じパッチ・有限重み・元拘束を消費した。R5の非アフィン場は全4配置で両OS各441,702アサーション。[凍結Case来歴](vpg-case-reference.md)・[R5引継ぎ](pg-case-reference-handoff/verification.md)・[R5物理検証](pg-case-reference-handoff/physical-reference.md) |

## 実閉路の資源境界

R-02の4パッチ閉路は2か所に非整合辺を持つ。実PGの辞書、全厚支持と共通区間から元672行・551列組・CSR28,116項を作り、rank 248、L/U4,330/19,173を得た。固定容量は23,503項で成功し、23,502項で拒否された。既存MSVCの3巡、同一因子の前進／転置、一括補正と元有理の実空間監査は別の記録に保持する。計画の累積確保、外部limbの見積り、プロセスの実測ピークを別々に扱う。 [準備と再利用](pg-sparse-chart-series.md)・[容量境界](pg-sparse-chart-series/factor-capacity-review.md)・[補正後の実空間上界](pg-chart-refinement-series/physical-trace-verification.md)

開始審査のR-01はS10の測定・入力・同期と実験不確かさであり、IM-08の判定に残る。R-02の早期審査では、実PGグラフの疎なrank・fillと準備・再利用容量を見る。実務規模のトポロジー、時間積分を含む性能、実験不確かさは、この閉路系列から外挿せず、それぞれの実務受入で判定する。 [R-01とR-02](../../plans/implementation-start-review.md)・[作業分解の早期R-02](../../plans/next-package-work-breakdown.md)・[A5-06](../../requirements/practical-acceptance-plan.md)

### 失敗した所有物と準備済み状態を再利用する

同じ閉路の元拘束を23,502項で準備して失敗させた。返った入力で全28,116係数の列と有理値、全672元行の右辺・来歴、各拘束の物理IDと列配置を照合した。この入力を作り直さず、因子上限だけを23,503項にして再準備した。基準チャートの識別、rank、L/Uの実項数が復元された。

最初に成功したチャートには、計画の確保上限ゼロでのliftと、事前取消での転置を続けて与えた。いずれも確保課金ゼロで拒否された。その後、同じチャートと再準備したチャートで非零の自由座標・力を作用させた。失敗前の全位置551組、自由力303組、元行乗数672組と両者が厳密一致し、非零の仕事および再準備後の全元行残差も照合した。

| Release環境 | 合格アサーション | 不足拒否までの課金 [B] | 再準備 (s) | 全工程の累積課金 [B] | プロセスピーク [B] |
| --- | ---: | ---: | ---: | ---: | ---: |
| MSVC | 8,195 | 68,801,368 | 30.386 | 2,026,577,616 | 937,066,496 |
| GCC 14 | 8,195 | 68,801,368 | 14.664 | 2,044,563,520 | 944,881,664 |

48 GiBの同一計画に計上した試行を各環境で1回ずつ示す。この再試行は容量と所有の照合である。両計時の差をコンパイラー比較や実務規模の予測に使わない。Clang-Tidyでは、凍結ドライバーの既存optional参照に9件残った。同条件の元ドライバーには13件あり、追加の指摘はなかった。

再実行するには、[凍結ドライバー](pg-sparse-chart-series/run-series.cpp.txt)を別ディレクトリの `run-series.cpp.txt` にバイトを保ってコピーする。そのディレクトリで [回復差分](pg-sparse-chart-series/run-recovery.patch) を `git -c core.autocrlf=false apply --no-index --unidiff-zero -p0 <差分の絶対パス>` で適用する。製品ライブラリは審査対象の固定版をReleaseでビルドする。 [元のコンパイル条件](pg-sparse-chart-series.md#再実行)を使い、差分側に試験用mainがあるためCatch2Mainのリンクだけを省く。GCCは `-x c++`、実行時スタック64 MiB、Catch2とstructure/coreを使う。元・差分・合成したソースのハッシュと実測の内訳は[検証台帳](im01-pg-r02-recovery.json)に記録した。

同じ審査対象で、R-02の準備と再利用、容量境界、凍結Case製品参照、R4/R5物理検証の4台帳から実行元と配布物35件を照合した。各台帳が指定する実行時Git blobのSHAをソース記録と、審査対象のGit blobのSHAを配布物記録と比較し、すべて一致した。

## IM-02へ渡す型と残る判定

- `prepare_case_pg_chart` または `prepare_case_pg_motion_chart` の所有結果は `CasePgChart`。`CasePgShellReference` が物理Shellの `PgPreparedPatch`、`PgConstraintPatchLayout`、`PgSparseChart` と `PgOwnedDiscreteVolumePoint` の必須範囲を同じ所有状態から借りる。点、層、片側、$N/G/X_0$、物理測度の厳密結果は `PgConditionedReferenceBatch`。設定した精度に収まる場合の投影結果は `PgProjectedReferenceBatch` で、固定点の順序を識別する。体積則の重みと射影した $N/J$ は `PgPreparedVolumeReference` が所有する。[Case引継ぎ契約](../../specs/structure/pg-case-chart-contract.md)
- チャートの `PgGlobalConstraints` は元CSR、行来歴と右辺、列配置を保持する。`PgSparseChart` は位置、$J$、$J^\mathsf T$、元行乗数を同じ因子から返す。時変の場合は `CaseStructuralSupport` と `CaseStructuralSupportSample` が同じ因子・元行順に時刻と片側を結び、位置・速度・加速度を渡す。[疎チャート](../../specs/structure/pg-sparse-chart-contract.md)・[構造支持の型](../../specs/structure/pg-case-structural-support-contract.md)
- 実Caseで同じ準備済み基底と拘束を消費できる。任意の変位について投影誤差を物理作用全体へ配分する検査と、製品作用の独立求積はIM-02で行う。有限則の合同成功だけでは、連続積分の精度や正の積則の全域到達を評価できない。凍結R5の17×17×3正積則は旧実行でGram準備中に停止し、固定統合版では高速化の先頭16点の端点一致までが記録されている。一般曲面での完走と実用的な資源・積分誤差は未確認。 [曲面積則の実測と範囲](pg-case-reference-handoff/verification.md)
- IM-02のV-PGでは同じ凍結入力と独立荷重から、材料・整合質量・変位・仕事・層別応力・固有対まで観測する。IM-00の [#16](https://github.com/YUKIKEDA/UVLAS/issues/16) は終了済み。IM-01で残る子Issueの終了と審査の承認後、[#18](https://github.com/YUKIKEDA/UVLAS/issues/18) の着手分解を経て本体へ進む。 [V-PGの分担](../../plans/pg-verification-gate-ownership.md)

数値補正を使う場合は、`PgChartCandidate::prepare` に `PgSparseChart&&` で所有を移し、`refine_lift` と `refine_transpose` に呼出し側の明示予算を渡す。実閉路の補正と全厚位置差の証拠は、この入口に原拘束から組み立てたチャートを渡したもの。`CasePgChart::chart()` は定数借用なので、同じCaseから同一のチャートを数値候補へ移す製品入口は未接続。Case所有状態からの厳密作用はそのまま使え、数値補正を後続で必須にする場合は所有の受渡しを実証する。 [数値補正契約](../../specs/structure/pg-chart-original-refinement-contract.md)・[R-02の実空間照合](pg-chart-refinement-series/physical-trace-verification.md)
