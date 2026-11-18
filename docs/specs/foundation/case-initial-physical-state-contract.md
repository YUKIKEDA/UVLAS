# 初期摂動の完全な物理場・状態と点値の整合確認

Issue #110、PR #111のFD-08に対する追加規約である。初期変位・速度を完全な物理場または物理状態として指定し、点値はその確認に使う。点から未指定の自由成分を補間・最小ノルム化して選ぶ機能は設けない。製品の目的・要素・精度・J1初期残差を変更しない。

本書は[フィールド辞書](case-schema-fields.md#fd-08-初期状態と再開要求)の`Perturbation physical`行と、その`displacements`/`velocities`に使う`PointValue`行に優先する。旧[設計例](case-schema-examples.json)の`fragments.Perturbation.physical`も、新しい完全指定例で置換して読む。旧ファイルと証拠はそのまま保持する。`Initial cold`、`Initial prepared`、`Perturbation none`、再開の契約は維持する。

本書と[新しい設計例](case-initial-physical-state-examples.json)は入力・写像・検証の契約である。typed validator、SI変換、PGの準備、連成初期化が実装・認定済みであることを表さない。

## 初期状態を決める入力

以下のkeyは全て必須であり、追加keyを拒否する。`Number`、`Rational`、ID、単位と元tokenの保持は[全schema契約](case-schema-contract.md#cs-02-字句型と閉じたobject)を継承する。

| 型 | key:型 | 条件 |
|---|---|---|
| Perturbation physical | kind:Tag{physical}, shells:List<ShellPhysicalState>, beams:List<BeamState>, rigids:List<RigidState>, displacements:List<PointCheck>, velocities:List<PointCheck> | 3種の状態集合で対象modelを完全被覆する。点確認の2配列は空可 |
| ShellPhysicalState | shell:Ref<Shell>, displacement:FullShellField, velocity:FullShellField | 両方とも同じ平衡からの増分。前者L、後者L/T。零場も明示する |
| FullShellField global_polynomial | kind:Tag{global_polynomial}, r:VectorPolynomial, p:VectorPolynomial, c:VectorPolynomial | 無次元chart座標(s,t)の多項式。3場を省略しない |
| FullShellField half_cells | kind:Tag{half_cells}, cells:List<CellField> | 正面積の活動half-cellを各1回だけ被覆する |
| CellField | cell:Vec<Index,2>, r:VectorPolynomial, p:VectorPolynomial, c:VectorPolynomial | cell=(c_s,c_t)の局所(u,v)係数。3場を省略しない |
| VectorPolynomial | terms:List<VectorTerm> | 空termsはexact零。元の項順を保持し、同じpowersの重複を拒否する |
| VectorTerm | powers:Vec<Count,2>, value:Vec<ExactQuantity,3> | valueは世界Cartesianのx/y/z順。係数ごとに使用側のLまたはL/Tを持つ |
| ExactQuantity | Number または Rational | JSONのnumber kindとobject kindで区別する。Rationalは既存のnumerator/denominatorだけを持つ。自由なkindや単位keyを足さない |
| PointCheck | point:InitialPoint, value:Vec<Number,3>, tolerance:Vec<Number,3> | value/toleranceはdisplacementsでL、velocitiesでL/T。各tolerance≥0、同じ配列内の重複pointを拒否する |
| InitialPoint shell | target:Target{kind=shell}, coordinates:Vec<Number,3>, side_s:Tag{interior,lower,upper}, side_t:Tag{interior,lower,upper} | coordinates=(s,t,ζ)、全て無次元。s/tは[0,1]、ζは[−1,1]。両片側tag必須 |
| InitialPoint beam | target:Target{kind=beam}, coordinates:Vec<Number,3> | 先頭は梁の参照物質座標s[L]、後2成分はexact零。区間内への所属はprepareで検査 |
| InitialPoint rigid | target:Target{kind=rigid}, coordinates:Vec<Number,3> | 3成分とも剛体の物体座標[L]。対象はMass.kind=rigidのみ |

`shells`は全Shell、`beams`は全Beam、`rigids`は`Mass.kind=rigid`の全Massを各1回含む。対象種別違い・欠落・重複を拒否する。外部sourceの解決後に対象集合が定まる場合は、prepareでその実集合も照合する。空配列は対象modelのその集合が空の場合だけ許す。粒子質量に独立状態を追加せず、hostの同じ物理点から導く。

`BeamState`と`RigidState`は旧FD-08の完全な型を使い、**摂動後の絶対状態**を与える。梁ではmodelの2節点と同順で両端の位置・姿勢・並進速度・角速度を全て指定する。剛体でも同じ全量を指定する。取り付けで従属する剛体の値も省略せず、hostとの整合を検査する。入力の完全性は独立DOFの追加を意味しない。

`InitialPoint`は`target.kind`で許可keyを決める。beam/rigidに`side`や`side_s`/`side_t`を許さない。shellの片側は[PGグラフ規約](../structure/pg-graph-cell-contract.md#所有と片側)と同じであり、各軸の厳密cell内部ではinterior、内部knotではlower/upperを指定する。0ではupper、1ではlowerだけが内向きである。両tagの型はschemaで、選択cellと向き・物理領域への所属はprepareで検査する。ζは厚さ座標であり、両chart軸の片側tagとは独立である。位置・速度は層界面でも共通の場なので、層stressのlower/upperを指定しない。この新しい点型の適用はPointCheckだけであり、旧MaterialPointの他consumerを本規約で置換しない。

## シェルの全厚さ位置場を指定する

幾何、活動領域、厚さ、参照尺度、アンカー、基底の版は、対象Shellの同じ準備参照から得る。利用者に別のg、厚さ、法線を指定させない。 $$ g=1+f_X^2+f_Y^2>0 $$ 、 $$ \widehat z=H(s,t)\zeta/(2L_{\mathrm{ref}}) $$ とする。

$$
\delta x(s,t,\zeta)=\delta R(s,t)
 +\widehat z\frac{\delta P(s,t)}{g(s,t)^{3/2}}
 +\frac{\widehat z^2}{2}\delta C(s,t),
\qquad
\delta v(s,t,\zeta)=\delta R_v(s,t)
 +\widehat z\frac{\delta P_v(s,t)}{g(s,t)^{3/2}}
 +\frac{\widehat z^2}{2}\delta C_v(s,t).
$$

各R/P/Cの係数は、変位場でL、速度場でL/Tである。Pは厚さ勾配の分子であり、directorそのものではない。物理初期場は $$ x_0=x_{\mathrm{eq}}+\delta x $$ 、 $$ v_0=v_{\mathrm{eq}}+\delta v $$ 。速度場に時間微分を掛けて変位場から生成したり、独立な速度入力を変位入力で置換したりしない。

`global_polynomial`は(s,t)、`half_cells`は(u,v)の昇冪係数を表す。局所座標は $$ s=(2c_s+1+u)/(4n) $$ 、 $$ t=(2c_t+1+v)/(4n) $$ 。各cellの分母gも同じgraphをこの座標へ引き戻す。総次数上限はRが13、Pが21、Cが2である。これは受入可能な形の必要条件であり、有限空間への所属を保証しない。例えばCのs²は次数上限内でも、既存の区分双一次hat空間へ所属しないため拒否する。

half-cell集合は、背景cellと物理領域の共通部分が正面積となる全cellである。辺だけの接触は活動cellに数えず、除外の根拠を保持する。利用者指定のactive boolを証明へ使わない。cellの欠落・余分なcell・重複・範囲外を検査し、片側係数を平均しない。同じ場に必要な連続性は、次節の全係数の所属検査でも確認する。入力の欠落を零cellで埋めない。

全項数、scalar数、指数、cell件数、typed payloadとscratchを既存`SchemaLimits`およびprepareの資源計画へ合算する。整数積と添字はchecked計算とする。元NumberとRationalをBinary64へ丸めてから有理係数へ復元せず、SIの正確な長さ・時間換算と元値を保持する。

## 有限場から準備済み係数へ写す

全世界成分について、同じ全活動cellのR/P/C係数行を共通順へ展開する。入力場の列をb、既存exact独立基底Bの係数行列をA_Bとする。原辞書の重複列を、未指定自由成分の選択へ使わない。

$$
A_B\beta_\alpha=b_\alpha,\qquad \alpha\in\{x,y,z\}.
$$

既存[PG係数契約](../structure/pg-coefficient-contract.md#4-rankと必須場を係数上で保存する)の独立化順とexact再構成を使う。全cellの係数恒等式が成立し、Bが独立であればβは一意である。このβを採用済みの条件付け変換へ移し、同じ参照の準備済みn_s×3係数を作る。変位増分と平衡係数の加算、および速度増分の加算でもN4の保持精度・誤差包絡を維持する。

点collocation、点最小二乗、最小ノルム、暗黙の近似射影によって非所属場を受理しない。所属しないことが確定した入力は`InvalidInput`。必要なrankを確定できなければ`RankUnresolved`、変換・保持精度を認定できなければ`PrecisionUnresolved`。確保・作業上限へ達した場合は、既存の`ResourceExceeded`を保つ。未条件付けraw辞書の完成だけを、この変換の完成と数えない。

## 梁と剛体の姿勢・速度を保持する

位置と並進速度は世界座標、回転行列Rは物体座標から世界座標への写像である。angular_velocityは物体座標の物理角速度を入力し、rootのangle/timeからrad/sへ変換する。回転行列へangleの倍率を掛けない。

$$
\dot R=R[\omega_{\mathrm{body}}]_\times,\qquad
\omega_{\mathrm{body}}=J_r(\theta)\dot\theta,\qquad
\dot x_c=v_p+R(\omega_{\mathrm{body}}\times c).
$$

梁の速度を局所回転チャートへ渡すときは、[梁採用仕様](../structure/beam-adopted-specification.md)の同じアンカーと右Jacobianを使う。非零姿勢で角速度配列をそのままチャート速度へコピーしない。SO(3)、相対Log分岐、接続、剛体取付の整合を検査し、回転行列を正規直交化して入力を修復しない。チャート不適格は既存の`ChartSingular`等の原因を保持する。

## 点値は最終保持状態の確認に使う

各PointCheckは、同じ物質点・同じ時刻での「摂動後−平衡」の世界Cartesian成分を指定する。displacementsとvelocitiesを別に扱う。点から係数、梁の姿勢、剛体の角速度を選ばない。

入力した完全場と点指定のexactな整合、および最終保持q/vからの評価の資格を分ける。後者の評価差をd、SI変換・状態変換・場評価・差引き・比較に伴う保守的な誤差上界をUとすると、成分ごとに次を満たす必要がある。

$$
|d_\alpha|+U_\alpha\le\mathrm{tolerance}_\alpha,\qquad U_\alpha\ge0.
$$

宣言した完全場自体が点指定と矛盾すると確定すれば`InvalidInput`。変換・保持・評価の誤差が明示された予算を満たすと認定できなければ`PrecisionUnresolved`。零toleranceはexactな証明がある場合だけ成功でき、一般の有限精度経路で成功を約束しない。大きな絶対位置の差引きでも低位語を保持し、入力場を直接評価した値だけで最終状態の検査を済ませない。

| 原因 | 既存FailureCode |
|---|---|
| key/shape/ID/完全被覆/単位/点予算の不正、確定した非所属場、非SO(3)の入力、位置・速度・取付・点値の矛盾、片側指定の不適合 | InvalidInput |
| Model・参照・平衡・チャートの同一性不一致 | StaleGeneration |
| 有効な姿勢からのチャートが特異、相対Logの分岐不適格 | ChartSingular |
| 必要な活動rankを確定できない | RankUnresolved |
| SI/係数/保持/評価の精度、明示PointCheck予算を認定できない | PrecisionUnresolved |
| 事前容量・作業・確保の不足、取消 | ResourceExceeded、Cancelledをそれぞれ保持 |

既存の幾何準備など内側の処理が返す具体的な既知Failureは、その原因と詳細を保つ。例えば参照幾何そのものの不適格を、一律に新しい入力分類へ上書きしない。

## 初期候補を構成し検査する順序

1. 新型の全key・kind・shape・Number/Rational・次数・件数を検査し、型付きIDを解決する。物理状態の全対象被覆とPointCheckの参照先も検査する。
2. 同じModel・平衡枝・支持・流れ・荷重・後流政策の参照を固定する。外部source、活動cell、SI変換を確定し、元入力の出所を保持する。
3. 各シェルの完全場をexactに所属確認して一意な係数へ写し、梁・剛体の完全な絶対状態と組み合わせる。採用精度で保持したq/vと、同じ参照・チャートの独立座標を構成する。
4. 変更後の全状態で、元の位置拘束・速度拘束、全位置接続、剛体取付、姿勢と幾何の正則性を検査する。既存J1の許容尺度を適用し、不整合なq/vを無断変更しない。
5. 全PointCheckについて完全入力場の整合と、最終保持状態の評価差・誤差包絡を検査する。
6. 連成では、平衡から引き継ぐΓ、実後流、政策と物質ID、および変更後のq/vから作る付着幾何を同じ初期状態として検査する。Γを据え置いた状態が非貫通を満たさなければ拒否し、本規約に新たなΓの補完法を足さない。cold/prepared_historyで明示したΓ・全後流を扱う既存経路は維持する。
7. [J1採用仕様](../coupled/coupled-adopted-specification.md#3-始動イベント刻みと末端)と[最終微分規約](../coupled/j1-initialization-output-derivative-addendum.md)に従い、同じ輸送政策でa₀/Γdot₀/λ₀と算法始動値を構成する。非零摂動後は既存の適合するIVPとして初期化し、離散固定点を理由にa₀/Γdot₀を零へ固定しない。
8. 全検査を満たした候補だけを既存`commit_initial`へ渡す。失敗時はModel・元入力・平衡・既存確定状態を変更しない。候補の完全状態、参照・基底・チャート・平衡のID、精度と誤差、初期化の出所を所有して記録する。

この順序は、完全な初期状態を選ぶ規則である。点指定を満たす別の平衡を探したり、支持・接続条件を緩和したりする探索ではない。

## 設計例と製品で確かめること

[設計例JSON](case-initial-physical-state-examples.json)の`cases`は、入力・exact数式・姿勢と速度の独立照合用である。`context`と`expected`は設計資料のmetadataであり、製品のcaseに加えるkeyではない。`negative_cases`の`patches`はそのbase_caseに対する設計検査用JSON Pointerの置換/追加/削除であり、実行時の入力patch機能ではない。

| 例・境界 | 確認内容 |
|---|---|
| flat_clamped_global | g=1、Rz=s²/100、速度Rz=s²/50、P/Cは明示零。s=0の全厚位置と速度を保ち、(1,1/2,0)で1/100と1/50。各点予算は1e−12 |
| flat_full_thickness | R=0でもP/Cを別々に非零指定。ζ=±1で異なる変位・速度となり、中央面だけで完全状態を決められないことを確認 |
| flat_half_cells | 同じglobal場をn=1の全4cellへexactに引き戻し、global表現との全係数一致、欠落・重複の拒否を確認 |
| beam_rigid_nonzero | 梁根元I、先端R_x(π/2)、先端body角速度(0,0,1)、同じ先端frameの剛体を明示。物体点(1,0,0)の世界速度は(0,0,1) |
| 旧点指定だけ・全対象欠落 | 全状態を選ぶ情報の不足として拒否。未指定成分を零または最小ノルムで埋めない |
| 非所属C=s²・不連続cell・支持違反 | exact空間と元拘束への不適合を拒否し、入力値を近似・修復しない |
| 点値矛盾・負予算・片側tag欠落/外向き指定・非SO(3)・取付不整合 | 検査対象と失敗段階を特定し、旧状態が不変であることも試験する |
| rank/精度/資源・零予算 | 未解決と不適合を区別し、零点予算で正の評価誤差を成功扱いしない |

点値だけの旧型を拒否する負例には、レビューの反例が対応する。 $$ w_1=\delta s^2 $$ と $$ w_2=\delta s^2+\varepsilon s^2(t-1/2) $$ は同じ中央面点値と根元の値・傾きを持つが異なる完全場である。本契約では利用者が完全場を明示するため、この自由成分を実装者が選ぶ余地を残さない。

`root_cases.coupled_transient_nonzero`は旧root設計例をコピーし、初期摂動を新型へ置換したschema設計例である。構造入力場のexactな照合成功を、非貫通・実後流を含む連成初期化の成功へ読み替えない。製品の開始判定には、この型・写像・全検査と独立参照の審査が必要である。
