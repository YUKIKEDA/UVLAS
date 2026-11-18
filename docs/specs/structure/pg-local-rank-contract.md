# PG局所独立化：exact rank・全再構成・包含

2026-11-18。設計Issue #133、IM01-5。[係数契約](pg-coefficient-contract.md)4〜7節を、実装開始に必要な型・有限手順・資源へ具体化する追加規約。SH-D03/05/07、N4-02/08を適用する。元仕様・辞書・分割・凍結証拠は変更しない。数学契約と有限比較手順を固定する。5節の累積所有方式は7節の実入力下界により製品用には不採用とする。製品実装開始ゲートは閉じたままであり、資源方式はIM01-7へ渡す。製品rank本体の実装済み・実務認定を意味しない。

## 1. 入力、範囲、結合

入力はfactory生成済みの`const PgCoefficientDictionary&`と`const PgIntegrationPartition&`、全項必須のlimits、Workspace、ResourcePlanである。報告用activity span、callerのbool列、hashだけから資格を作らない。両ownerの実graph内容（x₀/w/a/f/n/convention/全cell）を比較し、相違はStaleGeneration、移動元はInvalidInput。比較にはcanonical有理数と実termを用いる。

活動cellは分割factoryの全cell記録から取得する。正面積だけを採り、零面積接触を含めない。完全性・cell順・重複/範囲・areaとactiveの整合を検査する。どんな小さい正面積も閾値で消さない。分割の厚さ・層・補助cutはrankの係数を変えない。rankのsource witnessはgraph、元Dの順序付き有理頂点、全cellの正面積根拠、係数recipe/独立化版を所有する。物理Q/o/Lref/Hや支持の資格へ広げない。

返すのは無条件付けの局所係数空間である。物理N/G、Gram、区間条件付け、求積、M/K、物理支持の核、PreparedPatch、全モデルの疎rankは後続。matrixは有界な1patchだけに作り、全体10万自由度へ密RREFを広げない。

## 2. 係数行、候補、保持する関係

$$ F $$ はfull Vの実列数、 $$ F_0 $$ は原列数、 $$ K_a $$ は活動cell数。n=1/2/4、原列上限122/246/530、full V上限228/451/933を維持する。係数行keyは(field,cell_s,cell_t,a,b)、fieldはR/P/C順、残り昇順。活動cell内の原列の非零係数の和集合だけを行にする。R/P/Cの総次数上限13/21/2から

$$
m\le(105+253+6)K_a=364K_a.
$$

全分子は同じgraphのgを使う。正面積開集合で零の多項式は恒等零なので、Dで切った点値や別分母の近似行列を使わない。厚さの三次数も独立であり、Rだけのrankを全位置空間のrankと報告しない。

候補の原辞書座標を $$ L=[G,I_F,J] $$ とする。Gの33列は既存required_mapの順（X/Y/f、j-majorの全域s曲率24場、法線3、並進X/Y/Z）。選択は33列→full V原列順。後ろのJは選択へ加えず、包含だけを検証する。3並進はスカラー係数が重複してもCartesian対象αを別タグで保持する。独立な必須場はBの先頭を占め、従属した必須場も座標を残す。

Jは原辞書各列から、そのanchor R/P追加列と既存Cのanchor_c_mapを引いた座標である。R/P両方が非零なら両追加列を引く。Cは既存Cへの構成mapを引き、新しいC列を作らない。元sourceの原列ordinalとfieldを使い、独立列ordinalで代用しない。Jは $$ F\times F_0 $$、候補全数は $$ N=33+F+F_0 $$。raw mapの重複項は所定の元順で足し、bit計画前に相殺を仮定しない。

$$
C=A L,\qquad B=C[:,b],\qquad
A=B S,\quad A G=B W,\quad A J=B Z,\quad B=A T.
$$

Aは活動係数行×full V列。Tは選ばれた候補の元辞書座標、S/W/ZはB座標。全関係を有理数として零残差で確認する。選択した列を後で並べ替えない。C専用列とR/P列の分離を検査し、rankは $$ r=r_{RP}+r_C $$。必須場の独立数、選ばれた候補ordinal、全従属列の座標を所有する。Bareissのpivot行履歴は比較手順の証跡であり、将来の別資源方式へ要求する公開結果ではない。

### 全背景anchorと活動商空間を混同しない

全背景の $$ T_aL_a=I $$ と $$ \ker(T_a|_V)=(I-L_aT_a)S_0 $$ は元辞書の恒等式である。全t half-cellのs=0係数traceで $$ T_a A_{\rm full}J=0 $$ を検査し、Jが元列−anchor持上げであることを由来・C mapと照合する。これにより元の核のDへの制限が $$ A J=B Z $$ で包含される。

D外だけに支えを持つ列は活動rankで消えても、全背景のanchor traceでは非零となり得る。したがって、全背景T_aが活動商空間へ降りることや、Zが物理支持の完全な核であることは主張しない。Dがs=0へ届かない例を負対照にする。

## 3. 整数化と決定的なfraction-free消去

入力有理数は既約p/q、q>0。各活動係数行に加え、全背景anchor trace行を別集合として作る。trace行keyは(field,cell_t,b)、上限 $$ \tau\le(14+22+3)2n=78n $$。u=-1への代入では各原termに $$ (-1)^a $$ を掛けて足す。trace行はrank判定へ混ぜない。

各行で入力項の分母積を $$ d_i $$、各L列のraw map分母積を $$ e_j $$ とし、分母1は掛けない。零列はe=1。分母積をLCMと称さず、全項を数えて作る。

$$
M_{ij}=d_i\, C_{ij}\, e_j
      =\sum_h(d_i A_{ih})(e_j L_{hj})\in\mathbb Z.
$$

各積/除算の整数性を確認する。行と列の非零尺度だけを掛けるので、rankと候補選択順は変わらない。候補scaleを戻さずにS/W/Zを返さない。

次の一つの手順を固定する。外側の候補順と内側の行順は2節の順である。浮動QR、確率的modular rank、閾値pivotを入れない。

1. 候補jの活動行整数をvへcloneする。既存pivotの順で $$ t_i=v[p_i] $$ を保持する。
2. pivot列の変換後整数ベクトルE_iと非零pivot $$ u_i=E_i[p_i] $$ を使い、未pivot行kを
$$
v_k\leftarrow (u_i v_k-t_i E_{ik})/u_{i-1},\qquad u_{-1}=1
$$
へ更新する。分子2積と差を作ってからexactに除す。商の分母が1でなければRankUnresolved。過去pivot行は零のまま。
3. 残りが全零なら従属。非零なら最初の非零行を次pivotとし、候補自体をBへ採る。Eには変換後vを、上三角Uには各t_iと新pivotを保存する。pivot行は発見順で、数値row番号へ再sortしない。
4. GとI_Fの後、Jの候補を同じ係数でreduceする。新pivotが必要なら包含証明失敗であり、Bを拡張せずRankUnresolved。
5. 各従属候補は、その時点のrank qと共通分母 $$ \delta_q=u_{q-1} $$（q=0は1）で上三角系を解く。整数分子nを後退順に
$$
n_i=(\delta_q t_i-\sum_{k>i}U_{ik}n_k)/U_{ii}
$$
とする。各除算の整数性を検査する。選択候補は単位座標。最終 $$ \delta=u_{r-1} $$（r=0は1）へ $$ n\leftarrow n\delta/\delta_q $$ をexactに換算し、後から増えた座標を零で埋める。
6. 全候補の整数座標N_cについて
$$
M_D[:,b]N_c=\delta M_D,\quad
B_{ik}=M_{i,b_k}/(d_i e_{b_k}),\quad
Q_{kj}=e_{b_k}(N_c)_{kj}/(e_j\delta)
$$
を活動行だけで検査し、QからW/S/Zを取り出す。M_Dは先頭m活動行、残りのM_aは全背景trace行であり、M_aの全候補をBで再構成する条件は課さない。B/Tとrequired/kernel mapの橋渡しも検査する。全背景traceのJ列はexact零でなければならない。

この消去はSylvesterの小行列式恒等式による整数保存消去を使う。[Bareiss (1968), Mathematics of Computation 22, 565–578](https://www.ams.org/mcom/1968-22-103/S0025-5718-1968-0226829-0/S0025-5718-1968-0226829-0.pdf)。ここでの候補優先順、skip-column、trace検査、mapの尺度復元は本契約の規約であり、論文が製品のこのI/Fを認定したとはしない。

## 4. 具体I/Fと所有

予定名は実装責任を示す。公開APIは英語の利用者向け説明を付ける。

| 型/処理 | 必須の内容 |
|---|---|
| PgLocalRank::plan | const辞書/分割、PgRankLimits→所有する寸法/bit/opcode/visit/payload/external予測。確保せず、Workspace/Planを変更しない |
| PgLocalRank::prepare | 同じ入力、Limits、Workspace&、ResourcePlan&、任意のcaller所有PgRankDiagnostic*→core::Result<PgLocalRank>、Operation::PrepareModel、noexcept |
| PgRankLimits | ExactLimits、最大活動/trace行、候補/独立列、行列slot/map entry、numeric opcode/topology visit、stop_token。全項明示。0は無制限でない |
| PgRankOrigin | RequiredTargetまたはOriginalColumnの閉variant。Cartesian対象とscalar indexを混ぜない |
| RowKey / BasisColumn / OriginalColumn / RequiredTarget / AnchorKernelColumn | 別型。行は(field,cell,powers)、列はそれぞれの空間を表す |
| basis coefficients | (BasisColumn,RowKey,ExactRational)の非零項、列/行順、field/cellごとのoffset/count。全項と由来を所有 |
| T / S / W / Z | それぞれ専用の型付きrow/columnを持つ非零triplet、shape、offset/count。Sはr×F、Wはr×33、Zはr×F₀、TはF×r |
| source witness | canonical graph x₀/w/a/f/n/convention、元D頂点、全cell area/active根拠、source provenance、G/J、recipe/結果版。Rational/term/native配列を所有 |
| report | m/F/r/r_RP/r_C/必須独立数、整数性/全再構成/traceの検査件数、値bitと作業bit、全call予測/実課金、任意phase timing |
| validate_for | const辞書/分割と明示Limits/Workspace/Planでsource実内容・全再構成を再検査。hashの比較だけで成功させない |

PgLocalRankはfactory専用constructor、default/copy不可、move-only。getterはlvalue限定readonly view、移動元はInvalidInput、borrowはowner移動/代入/破棄まで。入力破棄後もB/maps/診断は利用できる。reportのactivity spanやrank整数だけからqualified ownerを構築できない。

cache identityはgraph/D/全cell/係数recipeと選択規約版を含む。入力を変えた結果を同じ寸法だけで再利用しない。validate_forは所有witnessとの実値比較と現在入力への再構成を行うので、明示予算を必要とする。物理参照や反復点評価ごとにrankを再計算しない。

Workspace leaseはmetadata読取り/確保callback前に取得し、同Workspace再入を拒否する。別Plan/Workspaceの並行callは同じcanonical結果になる。失敗後も親/旧結果は不変、成功済み課金は返金しない。

InvalidInputは未知tag/shape/移動元、StaleGenerationは参照内容不一致、ResourceExceededは全quota/bit/byte/checked overflow/bad_alloc、Cancelledは取消、UnsupportedはRequireProvenBound。必要なexact整数性・包含・再構成を証明できなければRankUnresolved。既知の資源不足をRankUnresolvedへ隠さない。診断はstageとtyped row/column/不変条件を所有し、内側Failure詳細を上書きしない。

## 5. 有限比較手順のbitと資源（製品方式は不採用）

### 5.1 計画の入力値

M₀=m+τ、R=min(m,F)、N=33+F+F₀。Tρは活動行とtrace行へ読む元term出現数（高々全termの2倍）、KはLのraw map entry数。 $$ K\le35F+c^2,\ c=(n+1)^2 $$。形だけの全零rowは所有しない。planは有限row universeを列順のcanonical項から読み、mと分母bitを確保前に数える。各coefficient lookupは既存record/termへの二分探索であり、全row×F×比較上限をvisitへ含める。

bは借用する全値の最大bit。行iの分母積bit上限D_i、列jのmap分母積bit上限E_jは、1でない各分母のbit長の和（空積は1）。a_iは行iでの最大分子bit、h_iは一つの原列の同trace rowへ入る項数（活動行は1、traceは高々22）。 $$ \ell(x)=\lceil\log_2\max(1,x)\rceil $$ とし、

$$
W_A=\max_i(D_i+a_i+\ell(h_i)+1),\quad
W=\max(b,W_A,\max_j(W_A+E_j+b+\ell(K_j)+1)).
$$

K_jは当該mapのraw entry数。上界は相殺前を数える。原列/必須場/anchorの各列scaleを省かない。行列式のLeibniz展開から、階数R+1までのminorを

$$
H=(R+1)(W+\ell(R+1))+2
$$

で囲える。消去E/U、pivot、Cramer分子はH bit以内。後退代入の和、尺度復元、公開mapのcross-multiplyを含め、

$$
V=\max(2H+\ell(R+1)+3,\ 2H+2E_{\max}+1,\
 H+E_{\max}+b+\ell(K)+2,\ 2W+2),\quad B=2V+2.
$$

を全leaf結果/作業の数学上界とする。ExactLimits.maximum_value_bitsはV、work capは既存exact契約のword丸め/guard込みBと事前照合する。observed最大値でこの上界を置き換えない。

### 5.2 opcode上界

数はexact scalarの生成/clone/整数部取得/四則/比較1callを1とする。符号/bitの読取りとnative loopはvisitへ分離する。全積/和はchecked uint64で、Index/size_tへ変換前も検査する。

| 固定工程 | exact call上界 |
|---|---:|
| graph/D/活動witnessとmap clone、符号、±1 | I+2K+2 |
| 行/列分母積、raw整数、scaled map | M₀+N+7Tρ+6K |
| 全候補Mの構成 | 2M₀K |
| 候補clone・Bareiss・pivot rhs | Nm+4NmR+NR |
| 整数後退代入 | 2NR²+3NR |
| 最終共通分母への換算 | 2NR+R |
| 公開map/Bの尺度復元 | 3NR+2mR+K |
| 全整数再構成・公開値の橋渡し | (2R+2)mN+7NR+6mR+2RK+RN |

Iは所有witnessへcloneする入力scalar数。表の和Oを先にquotaへ照合し、各call前にも消費を検査する。どのcell/列/必須場も上限に合わせて省略しない。過去pivotや依存候補の再走査を追加する場合は表も更新する。

### 5.3 payloadと同時生存

native配列は元入力witness、行key M₀、候補origin N、source provenance F、map K、M[M₀,N]、E[m,R]、U[R(R+1)/2]、整数座標[R,N]、行分母[M₀]、列scale[N]、scaled map[K]、residual[m]、rhs[R]、出力B項[mR]、S/W/Z項[RN]、T項[K]、offset/countを各shapeから1回確保する。可変列は所有array+offset/countへ平坦化し、未計画のvector増設をしない。実際の出力は非零数を数えてcapacityを減らしてよいが、forecastから外す理由にしない。

各exact callは高々1個のRationalStorageを新規所有する。累積owned payload上界はnative全arrayのsizeof×capacityにO×sizeof(RationalStorage)を足す。sizeofはprivate純粋layout accessorの実値とし、推測の64Bを製品条件にしない。親の既存課金と残量を確保前に照合する。

同時exact値slotは上記scalar arrayの総数にwitness Iを足し、outer bodyの一時値を最大7個足す。7は公開mapの橋渡しでのp/q2個、左右2段積4個、差1個である。他工程はraw整数5、Bareiss4、後退代入のsumを含め6以下。新しい同時一時値を増やしたら計画も更新する。

外部limbは[既存exact資源式](pg-exact-algebra-resource-contract.md)の固定Boost/et_off/fresh非alias経路を用いる。guard込みlimb数をL、limb byteをhとし、同契約の関数をC_B/R_B/G_B/D_B/K_Bと呼ぶ。入力の実capacity合計H_inと同時slot数S_liveから、外部byte上界を $$ h[H_{in}+2S_{live}C_B+6C_B+R_B+\max(G_B,D_B,K_B)] $$ とする。2S_liveは各値の分子/分母、6C_Bはscalar leafの外側整数、R_Bはresize旧領域、最後の項はnested poolである。累積owned、同時limb forecast、実課金、プロセスRSSを別に返す。EstimateAdmissionであり、RequireProvenBoundは未証明Unsupported。

## 6. 限定設計検証と将来の製品検証

限定probeはFractionの通常有理Gaussを独立oracleとし、整数化/消去/尺度復元は別実装で照合する。手書きliteralも置く。製品から期待rankを生成しない。旧R&Dの浮動QRによる列選択を本優先順のoracleにしない。

- 全n=1/2/4、全活動cell、R/P/C、非単位分母・異なる行/列scale、必須33→原列順、重複/非零極小/不活動列を検査する。
- A=BS、AG=BW、AJ=BZ、B=ATと全公開mapの尺度を許容差0で照合する。旧元rank値やpassedを製品合格へ転記しない。
- s=0へ届かないD、零面積接触、正面積の微小片、同寸法別graph/D、全cellの欠落/重複/追加を含める。
- 必須場順入替え、列scale復元省略、C field省略、小係数切捨て、不活動行混入、anchorを活動商traceとみなす誤りを限定変異として検出する。
- 製品本体ではmoved owner、旧結果不変、全quotaの十分/equal/一つ不足、overflow、bad_alloc、取消、再入、別Workspace並列・timingのcanonical byte一致も必要。

独立レビューは式/実API、全call表、bitと同時slot、独立oracle、正常/負対照を照合する。限定probe成功だけで製品rank本体や実務性能を完了にしない。今回の限定oracleは手書き多項式行列9例で、全PG辞書や製品APIを実行しない。数学手順・尺度復元・anchor反例・6種の誤りの感度だけを確認する。上記の製品系列は後続DoDであり未実行。数学契約が閉じても、7節の資源不採用により製品rank実装の開始条件は満たさない。

## 7. 二つの資源障害と開始ゲート

現分割のn4 U字・2層・部分境界・refinementの凍結計測は、owned forecast 205,783,717,415 B、実累積33,873,559 B、外部forecast 19,978,141,584 B。現在prepareはowned forecastをResourcePlan.remainingへ確保前に照合する。48GiBをそのcapに使えば、その時点で拒否する。これを実プロセスが205.8GBを確保した結果と呼ばない。

本契約は高capで既に得られるqualified partitionからrankを設計・検証できるが、それを48GiB実務準備の成立へ読み替えない。分割admissionの再設計はIM01-6/[#136](https://github.com/YUKIKEDA/UVLAS/issues/136)の独立依存であり、旧V規約や現在の数値を本Issueで変更しない。実務/R-02 gateは本rankの契約・実装と、その資源設計/実装の両方を必要とする。累積call予算と同時プロセス48GiBは同じ指標ではなく、高capの使用を黙示の実務許可にしない。

### rank自身の累積所有方式を不採用とする根拠

既存n4四次曲面の固定記録 `docs/verification/IM-01/im01-3/n4-stage-observation.json` は全term数1,017,024、F=933、F₀=530である。gは各cellで総次数6以下なので高々28×64=1,792項。各活動係数行には原列ごとに高々1項しかないため、全背景活動時には

$$
m\ge\left\lceil\frac{1,017,024-1,792}{933}\right\rceil=1,089,\quad N=1,496,\quad R=933.
$$

5節のBareiss項だけで $$ 4NmR\ge6,079,965,408 $$ callとなる。RationalStorageが64 Bの固定layout条件では累積payload下界は389,117,786,112 B（約362.4 GiB）であり、48 GiBの51,539,607,552 Bを超える。他の工程を数える前から拒否する。これは実確保量でも実RSSでもなく、このadmission式の証明済み下界である。64 Bを仕様の普遍定数にはしない。数学上有限な方式であることと、製品予算内で開始可能であることは別である。

したがって、opcodeごとに可変長exact値を新規所有し累積課金する比較方式を製品用には採らない。型・所有・失敗の4節は意味上の引継ぎであり、資源実装まで確定した公開APIではない。再利用整数slot/scratchと、決定的modular discovery＋exact証明の比較・採否はIM01-7/[#140](https://github.com/YUKIKEDA/UVLAS/issues/140)に属する。未検証候補はここでは採用しない。候補選択順、全再構成、必須場とanchorの包含は変更しない。

分割側の現行48 GiB拒否・確保0 Bの実製品probeは[PR #138](https://github.com/YUKIKEDA/UVLAS/pull/138)の固定head `f278ee1423dd09e23f8c9f6a4d601482396ed8b5`、`docs/verification/IM-01/im01-partition-incidence-design.json`を参照する。本変更で重複実行しない。24 GiBは残予算が小さいため同じ事前比較で拒否し、別実行は不要である。

#133の実装開始までを含むDoDは未完である。本PRは数学契約と資源方式の不採用判断だけをレビューに出し、#133・親#17・R-02を完了扱いにしない。分割側IM01-6とrank側IM01-7の双方を閉じてから、製品実装の開始審査を行う。
