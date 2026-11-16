# PG参照係数の処理契約

2026-11-16。Issue #17の最初の実装閉包に対する追加規約。独立レビューを完了し、実装開始には7節の条件を適用する。IM-01全体の完了を承認する文書ではない。[主シェル採用仕様](main-shell-adopted-specification.md)のSH-D01〜03/05〜07を内部処理へ具体化する。採用した空間、物理範囲、数値予算は変更しない。

上位には[製品要求](../../requirements/product-requirements.md)、[受入精度](../../requirements/acceptance-criteria.md)、[構造統合](structural-adopted-specification.md)、主シェル採用仕様を適用する。式の詳細は[参照係数の定義](../../../rd/structure/main-shell/paired-geometry/paired-reference-coefficient-definition.md)と[局所／大域契約](../../../rd/structure/main-shell/paired-geometry/paired-reference-local-global-contract.md)を継承する。型・失敗・資源は[N4](../foundation/numerical-foundation-adoption-specification.md)、所有と実装済み型は[状態契約](../foundation/state-contract.md)に従う。上位と矛盾する場合は実装を止め、内部都合で上位の意味を変更しない。

## 1. 最初の閉包は未条件付けの係数辞書までとする

対象は有理数・区分多項式の演算、直接持上げの両族、基本r/d/c、参照補完、アンカー補完、活動域上の局所rankと必須場の包含である。結果は未条件付けの係数辞書を所有する。これを、物理評価の準備が全て済んだPreparedPatchや適格なModelとして返さない。

幾何内積による保証付き条件付け、法線管の全域正則性、実際のN/G・体積／面求積、パッチ間全トレース、疎chart・J/Jᵀは、同じ[IM-01](../../plans/implementation-start-review.md)の後続閉包に残る。全3D材料・物理M/K・慣性・非線形部材接続はIM-02、時間・強連成はIM-04が担当する。全必須能力が揃う前の部分公開は行わない。

最初の閉包は現在配置、荷重、材質、支持条件、時間履歴を入力に取らない。物理固定端を選んで基底を小さくする処理も含めない。これらに依存しない参照空間を後段へ渡す。

## 2. exact入力と配列の意味を固定する

### 2.1 有理数と多項式

係数の意味は既約な $$ p/q $$ 、 $$ p\in\mathbb Z $$ 、 $$ q\in\mathbb N_{>0} $$ とする。符号は分子だけに置き、零は0/1へ正規化する。既存のBoost Multiprecisionで表現し、新しい依存を追加しない。型の内部ABIや整数limb配置は本契約では固定しない。

内部入口は整数分子・分母を受ける。十進文字列を受ける補助入口は、符号・10進整数部・小数部・10の整数指数を有限の有理数へ厳密に変換する。非有限、分母零、指数や必要桁数の上限超過を拒否する。Binary64からの明示入口はその語の二進値を厳密に復号し、元の十進値と同一だったと推定しない。JSONや外部格子の一般入力変換はIM-06の責任である。

二変数多項式は非負整数指数対(a,b)と有理係数の有限和とする。同一指数の入力はexact加算し、零項を除く。標準順はa昇順、同じaではb昇順。零多項式は空の項列、零でない多項式の次数は実際の非零項から得る。指数や項数の上限を固定配列長から推測しない。

加減算・積・スカラー倍・代入・微分・有理点評価はexactに行う。積の指数加算、出力項数、分子分母の増大を確保前に検査する。微分の軸はs/tまたは局所座標を型・引数で区別し、未知の軸を既定値へ丸めない。積分は変数と下端を明示し、任意定数を黙って選ばない。

### 2.2 参照とセル

グラフは無次元座標 $$ (s,t)\in[0,1]^2 $$ 上の

$$
\widehat r_0=(x_0+w s+a t,t,f(s,t)),\qquad w>0
$$

とする。fは総次数4以下の有理多項式、x₀/w/aも有理数として所有する。物理Q/o/L_refと厚さは後段の参照準備が保持するが、辞書の同一性へ入力全体の内容識別を渡せるようにする。Qを有理数だから正確な直交行列だとみなす処理は行わない。

内部格子nは採用済み1/2/4。基本B-splineの区切りはj/n、持上げと係数評価の共通細分は各方向2nセルである。セル番号 $$ c_s,c_t\in\{0,\ldots,2n-1\} $$ の局所座標を

$$
u=4ns-(2c_s+1),\quad v=4nt-(2c_t+1),\quad (u,v)\in[-1,1]^2
$$

とする。従って $$ \partial_s=4n\partial_u $$ 、 $$ \partial_t=4n\partial_v $$ 。セルの線形順は $$ c_s(2n)+c_t $$ 。内部境界の評価には側タグを要求し、片側の係数微分を平均しない。物理アンカーtraceはc_s=0の片側多項式をu=−1へ制限して得る。そのt方向の区分係数を全sセルへ拡張する。別セルの局所u=−1はそのセルの左辺であり、アンカーと混同しない。全域の同じ多項式を別セル座標からアンカーへ代入する場合に限りu=−(2c_s+1)を使える。

活動域は主仕様の矩形または正面積多角形を維持する。辞書は背景全セルを持ち、rank判定へ使うセル集合は活動域と正面積で交わるものに限る。境界に触るだけのセルを活動セルへ数えず、小さい正面積を閾値で零へ落とさない。

### 2.3 活動mapは領域から専用factoryが生成する

`QualifiedActivityMap`は専用の`prepare_activity_map`だけが構築する所有型とする。constructorを非公開にし、呼出側がbool配列、cell番号一覧、自己申告のhash／成功フラグから資格済みmapを組み立てる経路を設けない。factoryの入力は元の領域定義、n、参照内容identity、資源上限である。元領域・全cellの交差判定・その根拠・使用した判定器の版を所有し、判定器のscratchへ参照しない。

矩形全体の場合は、この係数閉包自身のfactoryが全(2n)²セルを標準順に生成する。各交差は背景cellそのものであり、面積はexactに1/(4n²)。利用者から活動cellを選ばせない。

多角形の場合、同じIM-01の参照領域処理が供給するexact交差判定器を必須依存とする。判定器は元領域Dと各背景cell Bの実際の共通部分を求め、符号付き境界の内外規則と向きを保って、重複を除いた正面積領域のexact面積と零面積判定を返す。複数の連結成分が生じても全てを数え、頂点／辺接触だけを正面積へしない。独自の幾何方式を本書で追加する指定ではなく、SH-D07の共通部分・分割の処理責任である。

factoryは判定器に全背景cellを一度ずつ渡す。返値は各cellの番号、交差根拠と非負のexact面積を持ち、面積が正であることからfactoryが活動状態を導く。callerの活動フラグを信頼しない。全cell数／shape、重複、欠落、範囲外index、標準順への対応、元領域／参照／n／セル定義の内容hashを検査する。活動cellの追加と漏れをともに拒否する。全交差がDを重複面積なしで被覆することもexact面積の和と元領域面積で検査する。面積和だけを交差の局所的正しさの証明には使わず、各返値が実際のD∩Bである責任は交差判定器が持つ。

上位SHは多角形頂点を有理数に限定していない。有理数／有限十進の頂点はその値を使い、Binary64の頂点は明示変換でその二進語のexact有理値を使う。これは表現された領域をそのまま判定するためであり、頂点の丸め直し・snap・近接点統合を行う許可ではない。外部の実曲面や測定位置をこの有限表現へ変換する誤差は、幾何生成／IM-06側の既存責任である。不確かさ区間だけがあり境界の一致／内外を確定できない場合、中央値やtolで活動を決めない。

元領域の不正、shape／cell番号の不正はInvalidInput、領域の自己交差・重複等の幾何不適格はGeometryInvalid、参照内容不一致はStaleGenerationとする。判定に必要な幾何が足りない、または正面積／零面積の証拠が資源内で確定できない場合はRankUnresolvedを返し、補助診断に活動域根拠未認定と問題cellを記す。上限到達が判明した場合はResourceExceededを保つ。全矩形のmapや成功rankを代用しない。

多角形を含むrank処理の開始条件には、この交差判定器の契約と独立試験の完了を含める。未実装の依存を利用者boolmapで埋めない。矩形の係数／factoryを先に閉じても、多角形対応やIM-01全体の完了とはしない。ここでの活動判定は積分用の求積精度や物理Mの資格を与えるものではない。

## 3. 同じ位置空間を有限係数から構成する

### 3.1 共通分母と基本場

物理平面座標への微分を

$$
f_X=f_{,s}/w,\qquad f_Y=f_{,t}-a f_{,s}/w,\qquad g=1+f_X^2+f_Y^2
$$

とする。gは正であり、辞書のスカラー場を有理多項式の三組 $$ (R,P,C) $$ として保持する。

$$
\Phi(s,t,\widehat z)=R+\widehat z\,P/g^{3/2}+\widehat z^2C/2.
$$

Pは厚さ勾配の分子であり、Pそのものをdirectorとして評価しない。基本dの $$ B/\sqrt g $$ は $$ P=gB $$ により同じ表現へ入る。gの平方根や正規化法線を有理係数へ近似してrankを判定しない。

基本rは各方向の開いた二次B-splineで、端knotsは0/1を各3回、内部j/nを各1回とし、そのtensor積をi_s(n+2)+i_t順に並べる。参照補完R=fを末尾へ加え、原r群の数は $$ n_r=(n+2)^2+1 $$ 。X/Yは基本rからexactに再現する。

基本dの二次B-splineは内部j/nを各2回とし、tensor積Bを $$ P=gB $$ としてi_s(2n+1)+i_t順に並べる。その後に $$ P=(-g f_X,-g f_Y,g) $$ の3列をこの順で加える。原d群の数は $$ n_d=(2n+1)^2+3 $$ 。この3列は単位法線のCartesian成分を再現する。

基本cはj/n上の連続区分一次hat関数のtensor積で、順序はi_s(n+1)+i_t、数は $$ (n+1)^2 $$ 。cは独立な二次厚さ場である。材料層の積分分割や子patchの厚さ分割とは別の自由度を意味する。

### 3.2 s/t両族を直接積分する

s方向の原点セルoと次数kについて、曲率はその区間だけで $$ L_k(4ns-(2o+1)) $$ 、それ以外で零とする。o=0,…,2n−1、k=0,…,7である。Legendre多項式は $$ L_0=1,L_1=u $$ と通常の三項漸化式でexactに構成する。節点の片側は前節の側タグで決める。

$$
I_m(s)=\int_0^s \sigma^m k(\sigma)\,d\sigma,\quad b=sI_0-I_1,
$$
$$
\widetilde v_s=-fI_0+\int_0^s f(\sigma,t)k(\sigma)\,d\sigma,\quad
\widetilde v_t=-s\int_0^s f_{,t}(\sigma,t)k(\sigma)\,d\sigma+\int_0^s \sigma f_{,t}(\sigma,t)k(\sigma)\,d\sigma,\quad
\widetilde v_z=b.
$$

積分は先行セルの完了積分を定数として次セルへ渡す。セルごとに積分定数を零へ戻さない。t族はfとs/t、セル軸を入れ替えて同じ式を評価し、面内成分の順を元へ戻す。その後、両族とも

$$
v_X=\widetilde v_s/w,\qquad v_Y=\widetilde v_t-a\widetilde v_s/w,\qquad v_Z=\widetilde v_z
$$

を適用する。曲面全体の係数積分であり、柱面の曲率値を点別に置換する方式へ変えない。

各v成分の厚さ勾配は同じvから導く。 $$ v_{,X}=v_{,s}/w $$ 、 $$ v_{,Y}=v_{,t}-a v_{,s}/w $$ 、

$$
t_X=-f_Xv_{X,X}-f_Yv_{Y,X}+v_{Z,X},\qquad
t_Y=-f_Xv_{X,Y}-f_Yv_{Y,Y}+v_{Z,Y}
$$

と置き、その分子を

$$
p_X=-(1+f_Y^2)t_X+f_Xf_Yt_Y,\quad
p_Y=f_Xf_Yt_X-(1+f_X^2)t_Y,\quad
p_Z=-f_Xt_X-f_Yt_Y
$$

とする。 $$ d=p/g^{3/2} $$ は[係数定義](../../../rd/structure/main-shell/paired-geometry/paired-reference-coefficient-definition.md)の接平面／法線方程式を満たす。各j=X/Y/Zの $$ (R,P)=(v_j,p_j) $$ をスカラー列として採り、全3Cartesian変位方向へ展開する。材料応答を使って列を選ばない。

対の列順はj、軸（s=0,t=1）、o、kの辞書式順、番号は $$ n_r+n_d+j(32n)+\mathrm{axis}(16n)+8o+k $$ 。対の数は96n。原辞書はr→d→対→cの順とし、各列の由来タグを所有する。独立化後の添字と原辞書添字を混用しない。

### 3.3 アンカー補完とfull V

アンカーはs=0、固定多項式 $$ \ell(s)=1-s $$ とする。原辞書Sの各列のR/P/C分子をこの面へ制限し、非零のR分子、次に非零のP分子の順で $$ (\ell R|_{s=0},0,0) $$ 、 $$ (0,\ell P|_{s=0},0) $$ を追加する。各群内は原列順とし、exact重複もまず由来を残す。Pの拡張は元のg(s,t)の分母を使い、アンカーでのみg(0,t)と一致する。

cのアンカー像の同じ拡張は既存のbilinear c空間に含まれるので、新しいc列は加えない。この包含もexactに検査する。結果Vについて

$$
T_aL_a=I\text{ on }T_aS,\qquad
V=S+L_a(T_aS),\qquad (I-L_aT_a)S=\ker(T_a|_V)
$$

を係数で確認する。支持条件による核Hは後段で取り出す。原R/Pのクランプ済み列だけを製品のfull Vとして返さない。

## 4. rankと必須場を係数上で保存する

係数行のキーは(field,cell_s,cell_t,a,b)、fieldはR/P/Cの順とする。行列は必要な非零項だけから構成できる。異なる分母を丸めた点値の行列をexact係数行列と呼ばない。各厚さ次数が独立し、各活動セルに正面積の開集合があるため、そこで零となる多項式分子の係数は零である。この有限多項式の性質を用いて活動セル行のrankを判定する。

局所で $$ A=BS $$ を満たすexactな独立列Bと原列の座標Sを所有する。rankはBの列数であり、再構成誤差は有理数として零。全体10万自由度へ密RREFや密零空間を広げない。

必須の滑らかな場は、X/Y/f、s族の全域曲率 $$ L_k(2s-1) $$ （k=0,…,7）の各3成分、法線3成分の計30列である。順序はX/Y/f、成分jごとにk、法線X/Y/Z。全域曲率を各区間へexact再展開した係数mapを作り、原辞書からの包含を検査する。定数並進も別の必須場として検査する。

独立選択は必須場をこの順でexactに追加し、その後に残る原列を原順で追加する。従属と判定した列の座標を必ず残す。小さい浮動小数点方向を落とす理由にしない。Bの先頭は必須場の独立集合とし、後段の条件付けでもこの部分空間を識別できる情報を渡す。r/cの由来と座標変換を取り違えてrだけのrankを全位置空間のrankと報告しない。

## 5. 所有する返値と失敗保証

内部処理名は役割を示すものでABIの固定ではない。

| 処理 | 入力 → 所有する返値 |
|---|---|
| construct_coefficients | 検査済みグラフ有理係数、n、内容識別、上限 → 全セルのR/P/C分子、g、由来、原列順、アンカー補完 |
| differentiate_coefficients | const辞書、明示軸・階数 → 同じ参照identityとexact微分係数。元辞書を変更しない |
| restrict_anchor | const辞書 → 原traceの分子、元列対応、L_aと包含証拠 |
| prepare_activity_map | 元領域、n、参照内容identity、上限 → 全cell交差証拠と活動集合を持つQualifiedActivityMap。多角形はexact判定器を必須依存とする |
| qualify_coefficient_space | const辞書、QualifiedActivityMap、上限 → 活動係数rank、B/S、必須場map、全再構成のexact零証拠 |

指数は非負だが、配列寸法・総列数は既存のchecked64bit寸法で扱う。配列化する場合は最終添字が最速とし、scalar係数とCartesian成分αを別の型で識別する。Cartesian展開後は上位仕様通り3i+αである。

返値は自分の項列・有理数・由来・証拠を所有し、Workspaceの一時メモリへ参照しない。const借用は同期呼出しだけで、rvalueからの借用を拒否する。同Workspace再入を拒否し、別workerのWorkspaceを共有しない。内容identityにグラフ係数、w/a/x₀、n、セル定義、アンカー、活動域、基底版を含め、同じ寸法を同じ参照とみなさない。

入口は既存core::ResultとOperation::PrepareModelを使い、例外を値へ変換する。形状・未知のタグ・分母零・次数不適合はInvalidInput、非有限入力はNonfinite、参照世代不一致はStaleGeneration、上限超過はResourceExceeded。exact計算が資源内で完了できず独立性を証明できない場合に擬似逆を返さない。分かっているModel/Patch/世代・資源単位をFailureへ補い、内側の既知詳細は上書きしない。失敗時は入力・既存辞書・親確定状態を変更せず、途中の辞書を成功型として公開しない。

## 6. 精度と資源の契約

この閉包の有理四則と係数恒等式に数値許容差を導入しない。零比較はexactである。整数overflowをラップさせず、丸めた小係数を零にしない。B64の点値やEigenのrankを正確rankの代用にしない。

[主仕様5節](main-shell-adopted-specification.md#5-求積と数値精度)の合同変換誤差≤1e−20と係数組合せ区間半径≤1e−25は後段の保証付き準備の条件であり、本exact辞書の完成だけで達成したと記録しない。後段は同じ辞書から保証付きで組合せ、保持表現と明示B64投影を別に認定する。物理Nの無次元性、q[m]、G[1/m]、実厚さ・L_refの換算も後段へ渡す。exact辞書の有理数をTwoWordの低位語へ置き換える仕様ではない。

利用者は最大項数、最大次数、各有理数／作業領域の最大整数bit数、最大局所列数、byte予算を明示する。これらは物理精度を下げる設定ではなく、必要量が超えたら失敗するための上限である。値をclipしない。積と原始関数で増える次数・bit数の保守上限を計算し、未知の実割当量を実測値として報告しない。

所有配列のpayloadは既存ResourcePlanへ正直に課金し、Boost内部の可変整数・分数・一時領域をその課金だけでhard capしたと主張しない。外部／可変長確保と最大bit数・総項数の帳簿を別に持ち、bad_allocをFailureへ変換する。全体のResourceProfile上限に対する資格は後続R-02の実graph・fill・RSS測定まで保留する。

準備工程は係数構成、アンカー補完、exact独立化、包含検証を分け、実際の項数・次数・最大bit数・payload/peak・時間を記録する。反復点評価のために毎回exact辞書を再生成しない。工程計測無効時の数値byte一致を守る。

## 7. 独立試験と開始審査

実装と別の式／別sourceで期待値を作り、原本R&Dは書き換えない。最初の製品試験は次の全系列を含む。

| 系列 | 固定する照合 |
|---|---|
| 有理／多項式 | 符号・零・約分・重複項・消去・非単位分母、積／微分／原始関数の独立手計算、指数／bit／byte上限・分母零 |
| 格子 | n=1/2/4の全セル、端点と内部片側、局所→s/tの4n因子、先行セル積分を残す対照、非単位w=4/5・a=3/20 |
| 持上げ | 総次数4までの全15単項式、s/t両族・全3成分・全8曲率。膜2条件と横せん断／法線条件・混合二階微分をexactに検査し、残る膜成分の非零例を保持 |
| 基本／補完 | B-spline連続性、R/P分母の意味、X₀/法線/並進包含、C独立場、full Vとアンカー核、支持を変えても元辞書不変 |
| rank | 重複／完全不活動列、極小だが非零の有理係数、活動mapの全cell完全性・追加／欠落／重複／範囲外・同寸法別領域hash・正面積／辺接触／極小切断片の独立判定、全原列と必須場の再構成exact零 |
| 境界 | owner破棄後の返値、別参照／同寸法世代拒否、資源不足と既存結果不変、未知精度・非有限入口・並列Workspace |

参照結果は次の読取り済み原本を使う。SHA-256はファイル全byteの値であり、製品試験の成功を表す値ではない。

| id／手掛かり | SHA-256 | 根拠の範囲 |
|---|---|---|
| paired-lift-exact-result、[record](../../../rd/structure/main-shell/paired-geometry/records/paired-lift-exact-result.json) | `92fb8778e93e52fca922f139300c8dbe23957a7b39f95619a7cfef01c1a18608` | 各128場×2チャートの恒等式・zero director負対照。全製品認定ではない |
| mass-paired-reference-result、[record](../../../rd/structure/main-shell/paired-geometry/records/mass-paired-reference-result.json) | `51ddf1467a311cd42d8960de64c2ca5b2dd9e79ee1f7c5237f246f80c60d64bf` | 参照再現・微分・法線管・全traceの代表照合。今回のexact閉包を越える検査も含む |
| general-paired-polynomials、[source](../../../rd/structure/main-shell/paired-geometry/experiments/general_paired_polynomials.py) | `dd82587536402eae97cdf69b1dc8c3520322a2e612bbea1e210e5d6d9673d03c` | lifting/director_numerator/cell_raw/cell_dictionary/coarse_source_mapを読んだ式・順序の根拠 |
| paired-exact-polynomials、[source](../../../rd/structure/main-shell/paired-geometry/experiments/paired_exact_polynomials.py) | `f8321bb2188970885c5daddb7f119a80343ceb5cdc2926a0900a619be12a46e8` | 有理演算、B-spline、セル座標の読取り根拠。二次面限定の旧cell_dictionaryを一般面の主経路へ戻さない |

別担当が原R&Dの恒等式と独立literalを確認し、t族削除・director符号・4n係数・積分定数・アンカー位置・小係数切捨て・必須場mapの欠落を隔離変異で検出する。過去recordのpassed値を製品結果へコピーしない。新fixtureは参照版と原本hashを保持し、再生成は別directoryへ出す。

本契約のレビュー、式／順序／有限係数表現の独立確認、専用活動map factoryと多角形exact判定器の責任・開始条件の照合、資源の事前検査とFailure対応、全系列の独立期待値と負対照が揃ったとき、この閉包の実装を開始できる。C++型名・補助関数配置は実装者が選ぶ。方程式・未知空間・意味のある値を再選択する穴が見つかったときは、実装で埋めず設計欠陥として返す。

閉包終了時に主張できるのは、この係数処理が本契約へ適合したことまでである。区間条件付け／全物理trace／疎chartの実装、V-PG、R-02の実結合準備・fill、物理材料や実務受入は、[実装開始審査](../../plans/implementation-start-review.md)と[受入計画](../../requirements/practical-acceptance-plan.md)の対応gateで引き続き確認する。
