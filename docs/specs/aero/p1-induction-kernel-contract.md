# P1誘導核の演算・全13変数・資源契約

Issue #119、既存IM03-0の「P1核・全13変数」を実装へ渡す追加契約。UI-D01〜06、N4-02/03/06/08/09を維持する。[誘導採用仕様](uvlm-induction-adopted-specification.md)2〜5節の核を変更しない。[局所packet](local-geometry-packet-contract.md)の後に適用し、共有・四辺形・全作用は[既存作業分解](../../plans/next-package-work-breakdown.md)の次段で扱う。

## 1. 全3基底・全13列を所有して返す

内部入口は evaluate_p1_basis(const LocalGeometryPacket&, P1Accuracy, P1Limits, P1Workspace&, ResourcePlan&, stop_token) → AeroResult<P1Basis>。型名は論理境界であり公開API互換性を固定しない。packetの移動元、context不整合、未対応精度を拒否する。snapshot/query ownerを保持し、元語・物質位置式・核版・半径・向き・LocalOnly資格を結果に結び付ける。成功・失敗とも入力と旧結果は不変である。

basis[k,xyz]は3×3、jacobian[k,xyz,j]は3×3×13、行優先。kは三角形頂点0/1/2のP1基底、jはpのxyz、x₀のxyz、x₁のxyz、x₂のxyz、εの順である。世界基底でBは[1/m]、Jは[1/m²]。μ[m²/s]は入力せず全3列を返す。ε=0でも半径列を省かず、非負半径の右微分を返す。

P1AccuracyはBinary64/TwoWordと有限非負の正規化全幅要求tau_basis/tau_jacobianを必須とし既定値を置かない。packetと同じ最大辺長ellを用い、全9/117成分で

$$
\ell\,\mathrm{width}(B_{ka})\le\tau_{\rm basis},\qquad
\ell^2\,\mathrm{width}(J_{kaj})\le\tau_{\rm jacobian}
$$

を認定する。零要求は真の零幅だけ成功。各成分は尺度付き中心・絶対半径を持ち、TwoWord中心は両語の厳密和である。最終丸めを半径へ含め、正の誤差をsubnormalへの丸めで零にしない。指数を分離して幅を比較し、保持不能・幅超過はPrecisionUnresolved。Binary64へ自動投影しない。

P1Basisは固定126成分、元descriptor、枝、離散求積差、work/資源報告を所有する。RoundoffEnclosureとQuadratureDifferenceを別型にし、後者を証明済みerror_upperへ変換しない。getterはlvalue限定でworkspaceの参照を返さない。

## 2. 元の低語と物理微分を保持する

元の厳密dyadic式からp−x₀、x₁−x₀、x₂−x₀、εを作る。巨大原点＋小さいlowをhiだけの座標差へ置換しない。packetの包絡中心だけで核を評価せず、元位置式を再消費する。導出pと第三物質重みを先に丸めない。

元exactなell²から $$ L=2^{\lfloor\log_2\ell\rfloor} $$ を整数bit長と比較で選び、差/Lとε/Lを評価する。したがって $$ 1\le\ell/L<2 $$ 。Lと原点はこの呼出しの数値条件付けの定数とし、選択算法を微分しない。jetのseedはp、各頂点、半径それぞれの対応列に1/Lを置く。x₀列も局所変数x₀へseedを置き、原点を動く未知量として引き直さない。値と全jetを最後に1/L倍する。この構成は元の13物理変数の微分であり、最大辺同率点やfar閾値の選択微分を含まない。

前進jetは値1＋微分13。各演算の値と導関数を同じ外向き区間演算で囲む。分岐境界の両側を独立面積分で検査し、切替差が事前精度を外れれば成功にしない。

incident点の散布は元物質重みを使い、同じ頂点の全出現で

$$
D_i u=u_{,x_i}+\beta_i u_{,p}+u_{,\epsilon}D_i\epsilon
$$

とする。最後の項は固定半径を明示したconsumerだけ零。局所13列を散布結果で上書きせず、全列とdescriptorを後続へ渡す。

## 3. 近接の安定式と理想核の外側極限

k番基底は $$ \mu_i=\delta_{ik} $$ 。採用仕様の境界閉路0→2→1→0と、面内項の+n向き辺0→1→2→0を区別する。e=x₁−x₀、f=x₂−x₀、G=(e·e)(f·f)−(e·f)²として

$$
g=e\frac{(f\cdot f)(\mu_1-\mu_0)-(e\cdot f)(\mu_2-\mu_0)}G
+f\frac{(e\cdot e)(\mu_2-\mu_0)-(e\cdot f)(\mu_1-\mu_0)}G.
$$

真のG>0でも区間が零を含めばPrecisionUnresolved。細長比閾値でGeometryInvalidへ分類し直さない。

辺a→bのLₑ,t,s,z=s−Lₑ,k²,rₐ,rᵦは採用仕様3節と同じ。同符号の場合は

$$
I_0=\frac{L_e(s+z)}{r_a r_b(s r_b+z r_a)}\quad(sz>0),
$$

他の場合だけ $$ I_0=(s/r_a-z/r_b)/k^2 $$ とする。符号は元dyadic内積から厳密に選ぶ。理想核の辺延長上でk²=0でも前者は有限。直線寄与が零でも幾何微分を消さない。正の重み $$ \mu_{\rm eff}=(\mu_a r_b+\mu_b r_a)/(r_a+r_b) $$ を使う。

面内項は採用仕様4節の $$ {\cal I} $$ 。各辺の対数は $$ \log1p(2L_e/g_e) $$ とし、

$$
g_e=|s|+|L_e-s|-L_e+
 \frac{k^2}{r_a+|s|}+\frac{k^2}{r_b+|L_e-s|}
 =r_a+r_b-L_e.
$$

値は左の安定形、微分は右の滑らかな式 $$ dg_e=dr_a+dr_b-dL_e $$ で作る。s=0やLₑ−s=0でabsの非微分可能性を混入させない。gₑの正性を区間で確認できなければ精度未解決。

h²=d²+ε²>0では $$ \Omega=2\operatorname{atan2}(hJ,\Delta) $$ 、法線係数 $$ (d/h)\Omega $$ を使う。元式でε=d=0でも、有限閉三角形の外側なら適用可能である。この枝は

$$
{\cal I}_{n}=d\,\frac{2J}{\Delta},\qquad
d{\cal I}_{n}=\frac{2J}{\Delta}\,dd \quad(d=\epsilon=0)
$$

という極限を使い、零の値と非零のd方向微分を区別する。εの右微分は零。有限閉三角形そのものの理想自己・非incident接触はpacketで拒否済みであり、この外側極限と混同しない。Δの区間が零を含めば精度未解決。

ε=0,d≠0はh=|d|の滑らかな片側枝で半径列は零。ε>0の面/辺/頂点自己はh>0で全13列を評価する。自己直線辺が零でも、面内項や閉路全体を零にしない。

## 4. 有界な区間演算でsqrt・log1p・atan2を囲む

公開精度はBinary64/TwoWordのみ。内部は固定256有効bitの外向きdyadic区間を使う。端点は符号付き整数仮数とcheckedな2冪指数。加算は小さい項の捨てbitを外向きへ加え、指数差を巨大左shiftにしない。積・除算・丸めの残りは整数商余りで検出する。値とjetの積/商則を全成分で外向きに計算する。

元式と厳密符号は既存checked32768bit整数、区間演算は同じBoost依存の固定checked1024bit整数を使う。新依存・動的精度・丸めmode変更は不要。内部256bitで要求幅を満たせなければ失敗する。指数overflowもPrecisionUnresolved。

sqrtは正の端点を[1,4)へ2冪正規化し、整数比較による257回のbit試行で256fraction-bitの隣接dyadicに囲む。等号だけ零幅。jetは $$ d\sqrt x=dx/(2\sqrt x) $$ 。下端零で導関数を囲めない場合は、3節の元式による極限枝以外で成功にしない。

log1p(x)は1+x>0を確認し、y=1+xを $$ 2^q r,\ 1\le r<2 $$ へ縮約する。小さいxを1への通常加算で失わず区間でyを作る。t=(r−1)/(r+1)として

$$
\log r=2\sum_{j=0}^{127}\frac{t^{2j+1}}{2j+1}+R,\qquad
|R|\le\frac{2|t|^{257}}{257(1-t^2)},\quad |t|\le1/3.
$$

ln2もt=1/3の同式で生成し、q倍の誤差を含める。2冪境界を跨ぐ区間は端点ごとに単調関数として囲む。jetは $$ d\log1p(x)=dx/(1+x) $$ で、縮約qは微分しない。

atan2(y,x)はy>0の連続枝である。符号・大小からr=min(|x|,y)/max(|x|,y)∈[0,1]へ縮約する。 $$ z=r/(1+\sqrt{1+r^2})\le1/2 $$ として

$$
\arctan r=2\sum_{j=0}^{127}(-1)^j\frac{z^{2j+1}}{2j+1}+R,\qquad
|R|\le\frac{2|z|^{257}}{257}.
$$

πはMachin式 $$ 16\arctan(1/5)-4\arctan(1/239) $$ を各128項と交代級数剰余で囲む。ここだけはhalf-angle縮約をせず、各引数そのもののatan級数を使う。π/2やπとの差へ誤差を伝える。象限境界を跨ぐ区間は連続なy>0で包絡する。jetは $$ d\operatorname{atan2}(y,x)=(x\,dy-y\,dx)/(x^2+y^2) $$ 。std関数へ未規定ulp保証を仮定しない。

全項を外向き計算し、解析剰余と演算丸めを別に加える。仮数積512bit、除算用拡大最大769bit、sqrt候補比較最大1024bit内へ正規化する。最終変換はhi、残差からlo、残余の上向き半径の順。checked失敗から有限値を捏造しない。

## 5. 遠方の離散求積差を丸め上界と混同しない

farは元exact述語の $$ d_T^2>16\ell^2 $$ だけで、等号は近接。採用仕様2節と同じ面積分を使い、

$$
\beta=((1-u)(1-v),u,(1-u)v),\quad
x=\sum_i\beta_i x_i,\quad dS=J(1-u)\,du\,dv
$$

でGauss–Legendre12×12と8×8、計208点を評価する。順はu外側/v内側/基底/xyz。各点でnormal/Jを含む13jetを作り、12次を返す。求積座標と重みに依存しないnormal、J、半径と定数のjetは同じfacetと照会の中で一度だけ評価し、両次数で共用する。未解決時に黙って解析枝へ切り替えない。

節点・重みは256bit dyadic区間の固定表。実装時に有理Sturm列で各Legendre根区間の一意性・全根被覆・順序を独立確認し、 $$ w_i=2/((1-z_i^2)P'_n(z_i)^2) $$ を囲む。[0,1]へ移した重みの和1と0〜2n−1次モーメント1/(k+1)の包含も検査する。binary64表を高精度へ拡張しただけの表を使わない。

各基底速度ベクトル、各Jacobian列を別々に比較する。 $$ V_k=\ell B_k,\ H_{kj}=\ell^2J_{kj} $$ とし、12次中心c、8次との差のノルム上界に両次数の丸め半径を加えたδで

$$
\delta\le\max(10^{-10},10^{-8}\|c\|_2)
$$

を全てに要求する。右辺のノルムは下端を使う。零近傍は絶対基準。これは採用protocolの全基底・全列の相対/絶対基準をSI尺度によらず適用する具体化であり、試作の全列一括normによるabs+rel判定より保守側にする。1節の要求幅認定は別に行う。

丸め包絡は固定Gauss12作用素に対する包含である。12対8の差は連続面積分との厳密上界ではなく未解決診断。この結果だけを連続作用の厳密error_upperへ使うconsumerは拒否し、既存IM03-0の空力数値予算gateで独立収束を扱う。

## 6. 全callのpayload・work・失敗保証

容量は式の静的命令列から計算する。近接の全枝を持つP_near、遠方の共通jetを持つP_far,setup、farの1点を持つP_sample、定数/元位置/認定のP_setupを作る。共通の辺長・端点距離・normal/Jは全基底で共有し、各命令は既存結果のindexだけを読む。1命令1Jet13 slotを割り当て、使用済みslotの削減で容量を楽観化しない。各命令列の節点数をNとすると

$$
S=\max(N(P_{\rm near}),N(P_{\rm far,setup})+N(P_{\rm sample})+18)+N(P_{\rm setup})
$$

がjetの同時slot上限である。18はGauss8/12各9累積jet。sampleは208点に再利用し、全点のjetを保持しない。Jet13は14区間、区間は2端点、端点は仮数/指数である。固定scalar scratchは各primitiveの命令列から同じ方法で求める最大slot R、元exact scratchも同時生存する元式命令から求めるEとする。全byteは

$$
S\,\mathrm{sizeof}(\mathrm{Jet13})+
R\,\mathrm{sizeof}(\mathrm{Interval})+
E\,\mathrm{sizeof}(\mathrm{ExactInteger})+
\mathrm{sizeof}(\mathrm{workspace\ header})
$$

でcheckedに集計する。命令数・slot・演算costは製品の同じ式定義をcounting型で展開した静的計画が持ち、上限らしい定数を手で置かない。節点の結果を一度だけ評価し、元入力/選択枝/式版と計画identityを照合する。512slot等の仮置き上限を開始条件の根拠にしない。

workの単位を区間基本演算1回、sqrtのbit候補1回、区間演算の外側にある元exact整数演算/比較1回とする。区間基本演算内部の固定1024bit操作は二重課金しない。jetを丸ごと1stepとは数えない。基本opcodeのcostは以下である。

| jet opcode | 区間演算stepの構成 |
|---|---|
| add/sub/neg | 値1＋微分13＝14 |
| mul | 値1＋微分13×(積2＋和1)＝40 |
| div | q=x/yを再利用し、値の商1＋微分13×(qとの積1＋差1＋商1)＝40 |
| sqrt | 値のsetup/257候補/丸めの展開cost＋2sqrtの積1＋微分の商13 |
| log1p | 値の縮約/128項/剰余の展開cost＋1+xの和1＋微分の商13 |
| atan2 | 値の縮約/128項/剰余の展開cost＋分母の積2/和1＋微分13×(積2/差1/商1) |

値側のcostも同じ基本opcodeで展開し、logのpower更新/項の商/累積和は1項3step、atanも符号反転を含む実命令を加える。sqrtは257候補、Machinの2つのatanは各128項、剰余式と端点変換も命令列へ含める。基本costの和をW(P)として、全call admissionは

$$
W(P_{\rm setup})+
\max\{W(P_{\rm near}),W(P_{\rm far,setup})+208W(P_{\rm sample})+W(P_{\rm accumulation/check})\}
$$

である。nearの各分岐は選ぶ一方のcostの最大を取り、未評価の枝を実課金しない。固定208点の各基底/xyzの加算もaccumulationに展開する。計画のsum/product/byte/indexを全てcheckedに計算する。

結果は126成分と診断/owner headerをPlannedBuffer<P1Storage>(1)で所有し、packet ownerは共有保持する。一回の累積owned payloadはworkspace準備分（必要時）と結果の実sizeof。準備済みworkspaceの評価中にscratchを追加確保しない。結果予約とleaseを構築前に確立し、失敗した自身の予約だけを戻す。max_stepsは全callのchecked累積上限で、不足は既存FailureCode::ResourceExceeded。各求積点・級数項・bit候補前に取消と残量を検査する。admissionは時間/RSSの証明ではない。
固定1024bit backendはallocatorを持たない。既存32768bit演算の外部scratchはpacketと同じforecastを別計上し、同時積を重ねない。owner/control block、allocator overhead、RSSの証明は持たないため、EstimateAdmissionとRequireProvenBoundを区別し、後者に十分な証明がなければUnsupported。

検査順はlease/移動元/context→要求値/表現→checkedサイズ/admission→予約→数値評価→全幅/far照合→所有返却。各反復前の取消も維持する。原因code/detailsを保持し、精度・資源・取消・bad_allocをGeometryInvalidへ変えない。例外をAeroFailureへ翻訳し再入を作用前に拒否。独立workspaceは並行利用でき、同じ入力/版/環境の値・半径・枝・stepはスレッド数によらずbyte一致とする。

## 7. 独立参照・負対照・実装開始判定

参照は境界式を共有せず、採用仕様2節の面積分を別の多倍長経路で積分する。入力・分割・次数・依存版・閾値をfixtureへ凍結し、細分化/精度増加差も基準内に入れる。微分は固定実参照配置のbarycentric積分へcomplex-stepを13列別々に適用し、normal/Jの変化も含める。abs/conjugateや実数専用枝をcomplex-step内部へ持ち込まない。

全3基底の全成分・13列で相対1e−8または絶対1e−10、代数1e−10、幾何共変1e−8を[採用protocol](../../../rd/aero/kernel/records/uvlm-induction-adoption-protocol.md)から維持する。近接の丸め包絡は連続面積分の高精度参照で包含を検査する。farの丸め包絡は同じ固定Gauss12を独立高精度で評価した値と照合し、連続面積分/complex-stepとの差は離散誤差の検査として分ける。

| 系列 | 入力と観測 |
|---|---|
| 全基底/近接 | 凍結protocolの平面/skew/非平面quad各子。ε=0/.02/.08、面外距離1e−2/1e−3/1e−5、正則化距離0。全3基底・全13列 |
| 理想外側極限 | T=((0,0,0),(2,0,0),(0,2,0)),p=(3,3,0),ε=0。法線微分とε右微分0。辺延長上p=(3,0,0)も検査 |
| 自己/incidence | LG-F01〜14、全支持、両quad子、全6置換。ε=0閉面上と非incident接触の拒否、ε>0面/辺/頂点の全作用 |
| 尺度/低語 | 巨大原点＋小さいlow、2^±900、通常doubleでε²がunderflowする例。保持不能時は失敗し零切捨て成功にしない |
| 分岐/保存 | far等号と隣接値、移動/回転/尺度1e−3/1e3、一定μの渦輪復元、向き/係数反転、物質点の全出現鎖則 |
| 演算包絡 | sqrt非平方/等号、logの微小x/2冪境界、atanのΔ符号/零、s=0/L。全幅と半径の境界、subnormal未満の真の誤差 |
| 所有/資源 | quota等号/1不足、取消、移動元、世代変更、再入、独立並列、結果/Boost確保失敗。入力と旧結果をbyte照合 |

隔離control成功後、低語除去、面内項欠落、自己全零、半径列欠落、h=0のd微分零、abs微分混入、far差無視、半径と全幅の混同、枝選択微分混入を変異し、compile/link成功後に検出する。旧R&Dの成功を製品結果へコピーしない。

本契約の独立レビュー後にP1核実装を開始できる。製品適合には全系列・変異・D/R・Linux sanitizer・対象静的検査・work/payload実照合が必要。設計成立と数値認定は別である。共有16変数/実複素μ/JVP/VJP、大域幾何/TE、循環、全荷重、後流、R-03実規模費用は既存IM03-0の後続gateに残す。
