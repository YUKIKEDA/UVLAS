# 公開検証資源の調査台帳

調査日：2026-10-29。状態：S10の説明書・小規模データ・モデルを取得して採用審査中。採用済みの検証ケースではない。

詳細は [S10採用審査](../../rd/acceptance/s10/records/s10-resource-audit-2026-10-29.md)。APIで公開ライセンスCC-BY-4.0と取得物のチェックサムを確認した。モデルには梁・集中質量・剛体結合・三角形があり、未対応のまま読み飛ばして利用しない。測定不確かさ、原データ、構造変換等の審査は未完了。

要求ごとの必要証拠と、S10のGVT・静的・解放応答・フラッタ候補への分割は [要求・証拠対応](../../rd/planning/requirements-evidence-map.md) に整理した。比較の目的を分け、独立評価と校正の重複を避ける。

## 方針

インターネットで入手可能な公開資料・データと比較コードを利用する。公開実験データによる妥当性確認、独立コード比較、解析解による実装検証を区別する。一つのベンチマークで製品の全対応範囲を検証済みとはしない。

## 候補

| 資源 | 公開元で確認できた内容 | 想定用途 | 採用前の確認 |
|---|---|---|---|
| Pazy wing / AePW3 Large Deflection Working Group | 大変形する柔軟翼の静的・動的空力弾性ベンチマーク | 変形平衡、変形に伴う振動特性、フラッタの実験比較候補 | 供試体の版、構造詳細、試験条件、不確かさ、対象流れの適用性、利用条件 |
| 10度後退Pazy翼 S10 公開データ | 地上振動、静的、初期条件応答、フラッタ等のデータとNASTRANモデル | 後退翼について静的から動的まで一貫した比較候補 | データ本文、モデルの材料・接続・質量・支持、取得条件、ライセンス。LCOのデータ存在だけでLCOを初回対応範囲に追加しない |
| SHARPy公式例 | Goland翼のフラッタ解析を含む例 | UVLM・連成・安定性の独立コード比較候補 | 採用版、構造モデル差、減衰、後流、空力条件。梁構造の例をANCFシェルの同一定式化として扱わない |
| Project Chrono ANCF shell 資料・実装 | ANCF積層シェル、ANS/EAS等の定式化資料と例 | 要素・材料・内力・質量・接線の比較候補 | 対象版、要素の仮定、積分・安定化、初期曲面・複雑接続への制限、利用条件 |

### アクセス条件の一次確認

NASA AePW3の案内ページには、データ保管場所へのアクセスについて担当者への連絡が必要と記載されている。したがって、その全データを無連絡で取得可能とは扱わない。公開ページの存在と原データの取得可能性を分けて管理する。現時点では連絡を行わず、直接公開されているデータを優先して調査する。

SHARPyの例一覧ではGoland翼のフラッタとT-tailの例を確認した。複数揚力面の比較候補として後者も調査するが、翼構造と空力条件の一致を採用前に確認する。

## 一次資料

- Technion Aeroelasticity Lab, Pazy wing：https://aerospace.technion.ac.il/lab/aeroelasticity/research/pazy-wing-a-benchmark-aeroelastic-very-flexible-wing/
- NASA AePW3 Large Deflection Working Group：https://nescacademy.nasa.gov/workshops/AePW3/public/wg/largedeflection
- Technion, Experimental Data for the 10 Degree Swept Pazy Wing (S10), v1：https://zenodo.org/records/19114860
- SHARPy公式例：https://ic-sharpy.readthedocs.io/en/main/content/examples.html
- Chrono ANCFシェル技術資料：https://www.projectchrono.org/assets/validations/FEA/ancfShellValidation.pdf
- Chrono ANCF 3423 API：https://api.chrono.projectchrono.org/classchrono_1_1fea_1_1_ch_element_shell_a_n_c_f__3423.html

## 次の作業

ケースごとに入力の再構成可能性、実験条件・不確かさ、数値の抽出可能性、使用条件を精査する。モデルとデータの対応が取れた候補について検証仕様を作り、結果を見る前に指標と合否条件を決める。フラッタ等の実験比較候補に加え、積層材、初期曲面、シェル交差接続、弾性・減衰支持の部品・構造検証を別途確保する。

## 有限翼の独立コード候補の原本監査（2026-10-30）

[AVL 3.52公式配布](https://web.mit.edu/drela/Public/web/avl/) のソースを取得し、版・ハッシュと7ファイルを [比較入力マニフェスト](../../rd/aero/finite-wing/records/finite-wing-reference-inputs.json) へ記録した。[有限翼比較仕様](../../rd/aero/finite-wing/records/finite-wing-load-reference-specification.md) では後退テーパ翼と上反角翼の入力を用意した。Web説明書3.40と同梱説明書3.50表記を区別し、3.52の翼面流下脚力OFF・微小コア・軸設定をソースで監査した。実行体・参照数値・残差と帯荷重抽出は未確保。独立コード比較準備までであり、実験妥当性確認や採用済みケースではない。

### AVL実行体・出力取得（2026-10-30）

[実行プロトコルと結果](../../rd/aero/avl-reference/records/avl-reference-execution-protocol.md) に公式Windows3.52実行体のハッシュと2翼・3修正試行を保存。入力読込後のフラグ変化を見落とした初回未達を残し、最終フラグ・軸・流速・密度・高精度FT/FS/FE・反復一致を確認した。FSBは従来精度、FEはdCpと座標等であり全循環／全要素力ではない。定量残差、ソースの再現ビルド、コア監査コピーは未確保。

### AVL原本の倍精度ビルド・計測（2026-10-30）

[原本と計測の照合](../../rd/aero/avl-reference/records/avl-source-measurement-experiment.md) で、専用GNU14.4.0環境によるソースビルド、元の行列・循環・非貫通残差、渦線力と作用点の取得まで完了した。初回の計測不適合を保存し、修正2版40条件＋原本20条件で14判定が成立。原本と公式Windowsの粗格子荷重も一致した。ただしコア0では横力・ヨーモーメントの対称性が破れ、採用保留。固定20帯・全格子系列とUVLAS比較は未実施。

### コア0の原因切分けと限定修正候補（2026-10-31）

[直線近傍核・鏡映検査](../../rd/aero/avl-reference/records/avl-collinear-kernel-experiment.md) で、不要な横力／ヨーを有限辺外の桁落ちに切り分けた。独立高精度参照とFortran回帰を経た修正候補AVL3.52/DP/zero-core-exterior-stable-v1を別版で保持する。原本微小コアの通常出力は不変、修正0コアは今回の20条件で対称性・残差・荷重和を満たした。元の未修正0コアは不適格のまま保存。全格子系列・固定20帯・UVLAS比較は未実施。

### 参照全格子系列・固定帯保存（2026-10-31）

[全100条件の参照検査](../../rd/aero/avl-reference/records/avl-fixed-band-convergence-experiment.md) を完了。修正候補の最大3072未知量まで弦・スパン・対角系列を評価し、帯集計の合力／モーメント保存と残差を確認した。ただし弦系列の総量漸近性、点集中による固定帯分布の収束は未達。製造分布の追加診断も保存し、参照精度受入とUVLAS比較は保留する。

### 一次荷重・帯境界適合の追加検査（2026-10-31）

[一次荷重と保存再構成](../../rd/aero/avl-reference/records/avl-linear-load-wrench-experiment.md) により、有限迎角の影響と観測の影響を分離した。[帯境界適合格子の100条件](../../rd/aero/avl-reference/records/avl-band-aligned-grid-experiment.md) は同じ物理翼・実行体を使い、元点の観測と支持辺再構成が各帯で一致した。それでも総量／分布の漸近性は未達のため参照受入は保留。元格子族と追加候補の未達をすべて保持する。


### 参照の線形求解感度と残る調査（2026-10-31）

[8系の追加診断](../../rd/aero/avl-reference/records/avl-solve-sensitivity-experiment.md) は修正版11条件が成立し、通常出力不変、Decimal80残差、条件推定・随伴・摂動と荷重台帳を確保した。求解による荷重変化は相対最大1.86e−15であり、前回の格子非漸近性は残る。公式LAPACK 3.12.1のDLACN2・DASUM/DCOPYとライセンスを外部参照用に保存し、製品へリンクしていない。次は元係数の独立高精度照合と循環の空間離散化を調べる。参照精度の受入とUVLAS比較は保留。


### 保存係数の独立高精度照合（2026-10-31）

[8系296,832係数の診断](../../rd/aero/avl-reference/records/avl-aic-generation-results.md) では、制御点の対角にも有限辺を含める元の規約で核を再現した。80桁との差は最大1.092e−13 m⁻¹、160桁自己監査・独立線積分と小系荷重も事前条件内。追加AVL実行や物理翼変更はなく、全係数を持つ小系と標本照合の大系を明確に区別する。大系の全体誤差上界、参照収束の受入、UVLAS比較は未確保。次は馬蹄渦と渦輪の基底を揃えた物理空間の循環場比較を進める。


### 比較基底・共通領域の診断（2026-10-31）

[馬蹄渦／渦輪と共通細分](../../rd/aero/avl-reference/records/avl-field-overlay-results.md) の処理検査が成立し、既存40格子を追加AVL計算なしで比較した。全域P0場、後縁トレース、固定320平均を区別し、後縁の弦方向4系列だけが漸近候補となった。元荷重の未達と参照受入保留は継続。次は製造投影による表現誤差と物理循環差の分離、トレースの符号付き積分の診断を進める。


### 表現試験と観測の相殺（2026-10-31）

[製造投影・後縁符号相殺](../../rd/aero/avl-reference/records/avl-projection-trace-results.md) を追加。既存40前線へ2場を投影し、24系列で既知誤差が減少しても方向候補0となった。実データ36比較の正負と領域分解は元の揚力差を再現した。これは受入条件を緩める証拠でも、製造場を実流れの誤差へ流用する根拠でもない。次にτ/Kの保存観測と固定帯モーメントを検証し、参照誤差予算を具体化する。


### 保存観測と追加スパン参照（2026-10-31）

[τ/K観測試験](../../rd/aero/avl-reference/records/avl-conservative-observation-results.md) が成立し、スパン側を [N320の2翼](../../rd/aero/avl-reference/records/avl-span320-results.md) で追加した。別容量8000版は既存4回帰で同一出力となり、新規7680未知量も残差・保存条件内。ただし次数／方向の未達が残り、参照は未受入。さらに細かいケースには [対称静的参照の縮約](../../rd/aero/avl-reference/symmetric-reference-reduction-design.md) を原本監査・全翼照合してから検討する。製品の一般非対称用途やフラッタの参照へ流用しない。

### 対称参照の全翼照合・N640（2026-10-31）

[縮約の実施結果](../../rd/aero/avl-reference/records/avl-symmetry-reduction-results.md) は初回14プロセスの幾何未達を保存し、独立分点とSPACERのπ修正後18プロセスで全翼等価性を確認した。[追加16プロセスと保存観測](../../rd/aero/avl-reference/records/avl-symmetry-fine-observation-results.md) で物理N640までを同じ版から確保した。M24N640の半翼7680未知数は約91〜94秒・ピーク約3.42GB。揚力／帯分布の条件付き推定は予算内だが、モーメントのスパン条件未達により参照は未受入。N1280候補は未知数15360・半翼帯640で両容量を超える。別の資源・求解経路・回帰条件を定義するまで未実施とする。原本、初回失敗、修正版、派生全翼のハッシュは [成果物監査](../../rd/aero/avl-reference/records/avl-symmetry-artifact-audit.json) で照合した。

### N1280候補の容量・資源経路（2026-10-31）

[容量回帰結果](../../rd/aero/avl-reference/records/avl-span1280-capacity-results.md) を追加。既存の同梱逐次LUを維持し、3容量定数だけを変えた別版で8ケース・66出力がバイト一致した。ビルド30.16秒、回帰261.08秒、回帰ピーク3,433,209,856 bytes。これは実未知数7680までの証拠であり、15360未知数のN1280は未実行。割当14,336,000,000 bytesとLinux／cgroup資源監査を保持し、次段階は開始有効空き20GiB・RSS16GiB・各1200秒／計2400秒で別に管理する。参照全体受入とモーメント予算は未完了。

### N1280の完了と参照精度の保留（2026-10-31）

[N1280結果](../../rd/aero/avl-reference/records/avl-span1280-results.md) を追加。最初の入力は局所KSMAX=500で停止したため記録を保持し、別版のKSMAX1000・全8入力LOAD・8回帰66出力一致を検証した。新規8計算は合計1655.53秒・保存159,092,311 bytes、最細各15360未知数は680.87／715.29秒・ピーク13,466,492,928 bytes。14場の保存観測6条件と成果物監査が成立したが、両翼のスパン5観測量は全て次数変化が未達。合算誤差予算・参照採用は保留する。同じ方式のN2560は密配列だけで約53.41GiBとなるため、新たな参照計算経路と判定根拠が必要である。

### スパン次数変動の既存データ診断（2026-10-31）

[診断結果](../../rd/aero/avl-reference/records/avl-span-order-diagnostic-results.md) の8処理条件・製造9件・拒否4件と70桁Decimalを用いた監査が成立した。追加AVL計算0件。全18総量差の帯間相殺比1、最終差の翼端二乗寄与68.56〜79.62%、4つの2項モデルの分布予測残差約72〜264%を保持する。原因を断定せず、尾部・参照採用は未評価。M12・N1280の配分変更4入力を別仕様で作成・監査したが未実行であり、M24の誤差予算を代替しない。

### 固定総数の配分感度と微小応答（2026-10-31）

[4計算結果](../../rd/aero/avl-reference/records/avl-tip-redistribution-results.md) を追加。M12・N1280・各7680未知数で総計412.65秒、ピーク3,433,582,592 bytes、原本ZIPと派生場合計109,336,508 bytes。7処理条件・6場対上限・独立監査が成立した。tip-poorの応答は翼端側に集中したが、tip-richの応答は微小で方向も対称ではなかった。補正偶力の成分診断も保存した。残差から微小応答の前進誤差を保証せず、次は半翼の循環・荷重求解感度の設計を具体化する。M24予算・参照採用・UVLAS比較は未完了。

### 対称求解感度の実施と次の2因子対照（2026-10-31）

[求解感度結果](../../rd/aero/avl-reference/records/avl-symmetric-solve-sensitivity-results.md) を追加。小系2件で元行列・通常／転置解・Decimal80残差を検証し、主対象6件を計637.84秒・ピーク4,377,976,832 bytesで実行した。全8件の元状態不変、原本回帰とZIP内66ファイル、56列の独立荷重重み照合、循環差の独立Decimal評価を監査した。M12配分応答に対する感度包絡は最大8.99e−7倍で、今回の方法変更は応答を説明しなかった。前進誤差全体・係数生成誤差・M24尾部は未保証。6144／9216未知数の独立配分4入力を準備したが未実行。
