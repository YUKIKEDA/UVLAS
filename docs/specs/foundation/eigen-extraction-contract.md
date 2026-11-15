# 固有抽出は元問題・実返却値・要求範囲で認定する

2026-11-16、Issue #64。N4-07の追加契約。Eigen 5.0.1による小規模完全法、対称一般化問題、Spectra 1.2.0による一般実／複素作用Arnoldiを対象とする。方式と品質を変更せず、実装者が選び直す余地のあった入出力と認定を定義する。製品コードと数値検証は、この設計PRでは追加しない。

本書は実装とレビューの正本であり、公開ABIを固定しない。型名は役割を示す。N4-02/03/06/08/09、[精度保持](precision-contract.md)、[外部メモリ](backend-memory-contract.md)、[参照fixture](reference-fixture-contract.md)を継承する。

## EIG-01：三つの原問題と入力を固定する

| 処理 | 原問題 | 入力と所有する出力 |
|---|---|---|
| `full_spectrum` | 一般の実／複素正方行列の全固有値 | 不変の原A、成分layout・尺度・単位、Context、精度・予算・容量 → 全n固有値、複素Schur packet、認定と資源診断 |
| `symmetric_generalized_spectrum` | 実対称Kと実対称正定値Mに対する全一般化固有値 | 別々の不変K/MとContext、共通layout・尺度、精度・予算・容量 → 全n実固有値、M正規化vector群、元残差・直交性・M資格 |
| `arnoldi_spectrum` | 同一の一般実／複素作用素Aの部分固有値 | 元作用・随伴、世代stamp、次数、固定初期vector、抽出要求、精度・予算・容量 → 所有Ritz値・vectorまたは部分空間、元残差、要求と実際の抽出数・作用数 |

一般問題の原式は次のとおり。左vectorは共役転置に対して定義する。実行列でも固有値・vectorの結果は複素数を許す。

$$
A\phi=\lambda\phi,\qquad A^*\psi=\overline\lambda\psi,
\qquad K\phi=\lambda M\phi.
$$

行列は列優先でshapeとleading dimensionを所有する既存のDense型を使う。入力は同期呼出し中の読取り専用借用とし、再利用するprepared packetは元係数とmetadataのsnapshotを所有する。K/Mの次数・layout・モデル・試行・時刻・尺度政策は一致させる。数値行列世代は個別に保持し、KとMの世代番号を同一にするための付け替えはしない。

既存のlinear::Contextは親・最後の確定状態IDを持たない。既知の場合は呼出側が別のStateDiagnosticContextとしてoptionalのAcceptedStateIdを渡し、処理・結果・失敗が値として保持する。未知はnulloptであり、ID0を代用しない。この識別は報告用で、数値基盤がAcceptedState実体の親子関係を検証した証拠ではない。

全ての数値配列・固有値・誤差上界は有限であることを検査する。NaN、Infinity、shape違いは入口で拒否する。特異A、負または零のK固有値は合法であり、特異な線形求解と混同しない。0×0の完全法は全0個の適格な空packetを返す。次数1/2の完全法も合法。Arnoldiは採用backendの次数条件を満たす必要があり、小次数で完全法へ黙って切り替えない。

## EIG-02：元問題の単位と固定尺度を保存する

一般問題は、各成分の正で有限な物理尺度を並べた対角行列Dと、正で有限な固有値尺度 $$ s_\lambda $$ を入力する。尺度・単位はModelと要求が先に決め、得られたvectorや残差に合わせて変更しない。全成分の次元は $$ [A_{ij}]=[s_\lambda] [\phi_i]/[\phi_j] $$ と整合させる。無次元作用、無次元固有値と座標は次で定義する。

$$
\mathcal A=D^{-1}AD/s_\lambda,\quad
\mu=\lambda/s_\lambda,\quad z=D^{-1}\phi,\quad y=D\psi.
$$

左作用は $$ \mathcal A^*=D A^*D^{-1}/s_\lambda $$ である。無次元行列を丸めて作った場合も、それを元Aと呼ばない。backendの候補生成に使った変換・balance・丸めと、原問題に対する認定を区別する。

対称一般化問題では、Dに加えて、 $$ \phi^T M\phi $$ と同じ単位を持つ正で有限な尺度 $$ m_0 $$ を入力する。構造問題では $$ m_0 $$ はJ·s²、 $$ s_\lambda $$ はs⁻²となる。単なるkgの尺度と取り違えない。

$$
\mathcal K=\frac{D^TKD}{m_0s_\lambda},\qquad
\mathcal M=\frac{D^TMD}{m_0},\qquad
\mathcal K z=\mu\mathcal Mz.
$$

基準結果は実際に返す物理座標の $$ \phi,\psi,\lambda $$ を正本とする。認定ではそれらから無次元量を評価する。backend座標zだけを認定し、D倍した返却値の丸めを無視しない。尺度変換の積・商が精度政策内で包めない場合はPrecisionUnresolvedで拒否する。

## EIG-03：原残差・正規化・後退誤差を別々に判定する

無次元の原残差を次で定義する。原係数、保持両語、実際の固有値を使う。縮約Hessenbergやshift後だけの残差を代用しない。

$$
\rho_R=\frac{D^{-1}(A\phi-\lambda\phi)}{s_\lambda},\quad
\rho_L=\frac{D(A^*\psi-\overline\lambda\psi)}{s_\lambda},\quad
\rho_G=\frac{D^T(K\phi-\lambda M\phi)}{m_0s_\lambda}.
$$

残差誤差を含む成分上界を $$ u_i\geq |\rho_i| $$ とする。原残差の指標は $$ \eta_R=\max_i u_i/\|z\|_\infty $$、 $$ \eta_L=\max_i u_i/\|y\|_\infty $$、 $$ \eta_G=\max_i u_i/\|z\|_\infty $$ とし、分母は正の下界を使う。非零vectorとそのnorm下界を証明できない候補は拒否する。これは成分尺度と固有値尺度に対する絶対残差であり、行列normで割った相対残差とは区別して記録する。

成分別後退誤差は、`NotRequested`または`Required{budget, observation, scope}`の明示要求とする。N4-07の固有残差と別の追加条件であり、線形求解の必須条件をそのまま固有問題へ転用しない。要求された場合の一般右問題では $$ d_i=(|\mathcal A|\,|z|)_i+|\mu|\,|z_i| $$、左問題ではAを共役転置、zをyに替える。一般化問題では $$ d_i=(|\mathcal K|\,|z|)_i+|\mu|(|\mathcal M|\,|z|)_i $$ とする。成分ごとの保守下界に対する $$ \beta=\max_i u_i/d_i $$ を認定する。分母も分子も厳密零ならその比は零、分母零で分子が非零なら未達、正の下界を作れなければ未認定とする。上界を分母として楽観的に割らない。追加要求が未達・未認定ならPrecisionUnresolvedで拒否し、NotRequestedへ変更しない。

明示行列では元係数から絶対積和の下界を作る。opaqueな作用でも全行列を組み立てず、三角不等式 $$ |(\mathcal A z)_i|\leq (|\mathcal A|\,|z|)_i $$ により、作用値の絶対値下界を使える。作用誤差を引いた非負下界とλ項の下界を合成するため、指標はより保守的になり得る。lower boundが零でも分子の厳密零を証明できれば比は零である。認定不能を隠すための行列norm推定や分母の上方補正は行わない。

Schurと部分空間の追加後退誤差を要求された場合は、残差行列 $$ E=\mathcal A Q-QT $$ の各成分上界を、 $$ (|\mathcal A|\,|Q|+|Q|\,|T|)_{ij} $$ の下界で割った最大値とする。opaque作用では前項と同様に $$ |\mathcal A Q| $$ の下界を使える。分母零の規則は同じである。pairの式をSchur列へ流用しない。

| 結果の種類 | 必須の数値条件 |
|---|---|
| 孤立pair | 右原残差・L2正規化、左を要求した場合はその全条件と一対一対応 |
| 一般化全basis | 全列の原残差・M正規化、全basisのM直交性・独立性 |
| 完全Schur／要求部分空間 | 原行列残差、Qの直交性・独立性、Tのshape／上三角性、要求個数 |

適用しないpair正規化などはNotRequestedであり、達成値0で代用しない。成分別後退誤差がNotRequestedならその分母を評価する必要はなく、追加指標だけの未表現・零下界を必須のnorm残差の失敗へ持ち込まない。原残差そのものの非有限・包絡不成立は引き続き拒否する。

入口のAccuracyは原残差、vector正規化、部分空間の残差・直交性をそれぞれ別の非負有限な予算とObservationId/ScopeIdで要求し、成分別後退誤差の追加要求も保存する。原残差の上限はN4-07の1e−10以下で、観測の5%割当からさらに厳しい要求があればそちらを使う。IM00-11の標準小系では原残差・正規化・直交性を各1e−12以下とし、成分別後退誤差はNotRequestedとする。解析的な単純根と凍結PGの固有値照合も1e−12級を維持する。一般のdefective問題に、残差だけから普遍的な固有値前進誤差1e−12を保証しない。

この区別の反例はA=diag(5,1)、λ=5、 $$ v=(1,\epsilon)^T $$ である。原残差のnormは $$ 4|\epsilon| $$ へ減るが、第2成分の後退誤差は非零εに対し2/3のままである。小さいtailを切り捨てて追加条件を通さず、要求があれば実候補を検査して拒否する。追加条件の正対照には厳密な対角固有vectorを使う。共通認定helperには同じε=2⁻⁸⁰の候補を与え、予算1e−12でβ=NotRequestedならnorm残差・正規化が合格、β=Requiredなら拒否することと、tailとgenerationの不変を確認する。これは採用前草稿への反例に基づく訂正で、既存のN4-07の基準を緩和する変更ではない。

積はEFT、加算は補償評価と外向き誤差包絡を使う。複素数は実虚両成分を評価し、複素絶対値とnormも下界・上界を区別する。underflow・非有限・表現不能な包絡は成功にしない。作用素が報告した成分誤差上界と、残差集約の誤差を両方加算する。予算0は厳密零を要求し、丸め幅や隠れた下限を足さない。

## EIG-04：vectorの正規化・位相は物理modeのIDではない

一般右／左vectorはそれぞれ $$ \|z\|_2=1 $$、 $$ \|y\|_2=1 $$ を要求予算内で満たす。位相は無次元成分の最大絶対値を持つ成分を非負実数へ回す。同値なら最小成分indexを選ぶ。実vectorではその成分を正にする。位相変更・正規化・物理座標への変換後の実返却値を再認定する。数学的に不要なbit一致を両OSへ要求しない。

固有値の表示順は実部、虚部の昇順とし、packet内のindexとしてだけ使う。同値な値の内部順やvectorの位相を物理mode IDにしない。実Aの非実根は共役対を保持し、複素Aへ共役対を強制しない。微小実部・虚部、負値・零を閾値でclampせず、平方根・周波数だけへ変換して符号を失わない。

孤立pairを要求した場合は、右と左の原残差をそれぞれ認定する。左右は固有値の一対一対応で結び、近接した値を重複消費しない。対応の政策は要求時に非負有限の絶対幅aと相対幅bを与え、無次元値間の距離が $$ a+b\max(|\mu_i|,|\mu_j|) $$ 以下のときだけ対応候補とする。予算と照合対象をpacketへ保存する。対応が一意でなければclusterとして要求し直すか未認定とする。群全体の照合では内部順を同定せず、重複度を保った完全matchingを要求する。この幅はassociationの政策であり、数学的分離や固有値前進誤差の上界ではない。試験の1e−12を全製品入力へ暗黙に適用しない。

$$
\kappa_{\rm pair}=\frac{\|y\|_2\|z\|_2}{|y^*z|}.
$$

この近似vectorによる感度は診断推定である。厳密な固有値誤差上界や臨界速度の不確かさへ自動変換しない。分母を正と認定できない場合はConditionUnresolvedを保存し、Infinityを有限の認定へ混ぜない。未要求・未評価も別の状態にする。

## EIG-05：完全Schur packetとclusterの資格を分ける

一般完全法は全n個の固有値と、無次元問題の複素Schur分解Q/Tを必須で所有する。実入力でも複素形式へ統一できる。上三角Tの対角から得た値と、物理単位の全固有値との対応を所有する。物理λへの変換丸めを含め、Tの対角を実返却λ/sλとして評価した場合のSchur残差と、要求された追加後退誤差も同じ予算で検査する。元のTの合格だけで、変換した固有値を認定しない。

$$
\mathcal A Q=QT,\qquad Q^*Q=I.
$$

原Aから作る残差の最大絶対行和を $$ \|Q\|_\infty $$ の正の下界で割った上界と、 $$ \|Q^*Q-I\|_\infty $$ の上界を要求予算と比較する。全n列、全成分、上三角性を検査する。列数が正ならGram誤差上界が1未満であることも必須とし、緩い要求予算で従属列を通さない。空packetの残差・Gram誤差は0。全数の存在とこの分解の適格性をFullSpectrumQualifiedと呼ぶ。有限精度の代数分解への資格であり、元の非正規問題の固有値前進誤差や物理的安定性の証明ではない。

clusterの抽出は、所有する完全packetと、そのpacketの固有値indexの重複しない明示集合を入力する別処理とする。indexは0以上n未満、空集合は空の部分空間である。基準のpacket identityと数値世代を照合する。異なるpacketの同じindexを流用しない。

標準的なunitary Schur交換で指定集合を先頭へ移し、k列の $$ Q_C $$ とk×kの $$ T_C $$ を所有する。原Aに対する $$ \mathcal A Q_C-Q_CT_C $$、 $$ Q_C^*Q_C-I_k $$、指定個数kと列の独立性を別に認定する。kが正ならGram誤差上界は要求予算以下かつ1未満とする。全Qの成功で一部clusterの成功を代替せず、単に元Qの該当列を抜かない。交換の不安定・次数不足・誤差予算未達はPrecisionUnresolvedとする。

左部分空間も要求する場合は、 $$ \mathcal A^* $$ に対して $$ T_C^* $$ と同じ根・重複度に対応する不変部分空間を別に生成・認定し、選択集合との対応を所有する。右Qの列や共役をそのまま左へ流用しない。対応を認定できなければPrecisionUnresolved。近い固有値の群分けは選択政策であり、数学的な分離の証明ではない。分離・交差・物理mode対応が必要なIM-05では、追加の分離診断・部分空間角・観測を要求する。部分空間の資格だけから分離を断定しない。

n=0またはk=0ではbackendとcallbackを呼ばず、空の残差・Gram誤差は0とする。非零vectorや正のnorm下界を空集合へ要求しない。

Jordan blockにも独立なn本の固有vectorを強要しない。個別vectorの全basisが要求された場合は、無次元列を並べたZの左逆候補Lに対し $$ \|I-LZ\|_\infty<1 $$ を丸め込みで認定して独立性を確認する。不足を全固有値取得と混同しない。基底内部の回転・位相が違っても、同じ部分空間を別modeと判断しない。

採用EigenのComplexSchur/RealSchur公開入口と、`unsupported/Eigen/src/MatrixFunctions/MatrixFunction.h`のSchur交換実装は読取りで確認した。公開reorder APIとして保証されているとは扱わない。実装Issueでは固定版に対する監査adapterか、同じ標準Givens交換の狭い実装を選び、元問題で再認定する。新しい抽出方式や依存をここで採用しない。

## EIG-06：一般化問題は原Mを認定し、元K/Mへ戻して検査する

K/Mとも全成分を厳密な数値等号で対称と確認する。下三角だけを読み、上三角の不一致を隠さない。近い値の平均化、ridge、擬似逆、負固有値方向の削除は行わない。

原Mの正定値資格には[疎基準契約](sparse-solve-contract.md)のRumpに基づく十分条件を使う。密入力を同じ実内容の検査用canonical CSCへ変換する費用も資源へ計上する。資格は元M snapshotと結び、丸めたscaled MのLLT成功で代替しない。正定値資格が作れなければRankUnresolvedで拒否し、それを原Mが負である証明とは呼ばない。負曲率・零modeの証人がある場合は既存の診断に保存する。

Eigenの対称一般化処理を用いて候補を作る。原K/Mを全成分検査した後のscaled候補は、片側の合同変換を計算して他方へ同値copyできる。元行列は保持し、原問題の認定を行う。原Mの資格と別に、backendへ渡す丸めたMのCholesky成功・有限係数・正対角も検査する。固定Eigenの一般化入口がその資格を全て検査すると仮定しない。丸めによるbackend前提の崩れはPrecisionUnresolvedとし、原Mの負曲率の証明にしない。Kへ正定値性を要求せず、M⁻¹Kを無断で自己共役化しない。全n値は実数のまま符号付きで保持し、返却vectorは次を要求予算内で満たす。

$$
\Phi^TM\Phi/m_0=I,\qquad
\mathcal Z^T\mathcal M\mathcal Z=I,\quad \mathcal Z=D^{-1}\Phi.
$$

M正規化とEuclid正規化を取り違えない。位相はEIG-04の無次元成分で決める。原残差は列ごとに、直交誤差は $$ \|\Phi^TM\Phi/m_0-I\|_\infty $$ の保守上界で認定する。全basisの独立性にはGram誤差の上界が1未満であることも必須とし、予算が1以上でもこの条件を外さない。重根の内部vector順を固定せず、M内積の部分空間として比較する。数値的な零と、厳密零・剛体modeの証明は別の診断にする。

## EIG-07：作用と随伴は同じ凍結問題を評価する

同期作用の入力は、n成分の高位・低位const spanと表現tagを持つ、呼出し中だけの読取り専用viewとする。所有PrecisionBufferを借用してviewを作れるが、backendの一平面から各作用ごとに所有bufferを構築しない。低位はBinary64では空、TwoWordでは高位と同じ長さとし、viewを呼出し後へ保存しない。出力は事前確保したn成分の非重複領域と、各成分の絶対誤差上界とする。callbackは例外なしのResultを返す。対応精度を入口で宣言し、hi/loを両方消費しない作用へTwoWordを渡さない。複素作用は複素線形とし、複素入力に対して実部だけを評価しない。出力と誤差の全成分の記入を検査し、前のscratchの零を未記入の成功へ使わない。

作用stampは元問題のContext、数値世代、layout、尺度、精度政策、変換／shift政策、backendを含む。原作用と随伴は同じstampを共有する。各callbackの前後と最終返却前にstampを照合し、途中変更はStaleGenerationとする。全期間で係数が不変という呼出側の契約も必要で、stampだけで外部メモリの改変を発見できるとは主張しない。明示行列のprepared再利用は世代と実内容の両方を照合する。

随伴は同じ内積と尺度の共役転置である。実装試験では複素内積の随伴一致も検査する。数回の内積試験だけで任意callbackの線形性・随伴性を証明したことにはしない。元P/P*の構成はJ1とIM-05の責任であり、数値基盤へ別の物理微分を実装しない。

backend callbackのadapterは最初の作用失敗を保存し、例外で越境させない。backendを中断できず以降に安全なdummy値を与える場合も、戻った候補を必ず破棄し、最初のFailureを返す。callbackの失敗をEigenCoverageUnresolvedへ一律に置き換えない。

## EIG-08：Arnoldi要求を黙って変えない

要求は元作用に対するLargestMagnitude、または指定shiftの近傍抽出を区別する。初期設定はN4-07のnev=24/ncv=64。要求と実効値を両方保存し、 $$ 1\leq nev\leq n-2 $$、 $$ nev+2\leq ncv\leq n $$ を確保前に検査する。小次元へclipしない。PG24では検証入力として明示的にnev=6/ncv=20を与える。

初期vectorは有限・非零の実／複素配列を呼出側が明示し、normと実内容hashを保存する。単にseed番号だけを保存してbackend版ごとの乱数列へ依存しない。同じ要求内では初期vector・元作用・逆作用を変更しない。前処理を使えるのはshifted線形求解の内部であり、通常のArnoldiを左前処理作用へ替えない。Krylov対象は元作用かEIG-09で明示したspectral変換だけとする。必要な追加試行・追加seedは全体予算へ数える。

backend内の許容差は候補抽出の設定であり、EIG-03の資格ではない。収束件数、nev/ncv、backend反復数、backend報告作用数、adapterが数えた実作用数を別に保存する。実Aの共役対境界では、要求nevを最低数とし、対を切らず最大nev+1個を返せる容量を先に確保する。複素Aへこの余分な対規則を適用しない。収束数が要求に足りない、対応が不明、最終原残差未達なら正常な部分抽出を返さない。

採用Spectraが必ずnev+1個を返すとは仮定しない。実Aで片側だけ得た境界pairは、元pairから共役候補を構成して別に元作用・要求随伴で認定し課金するか、適法な新しい抽出試行で得る。nev=n−2の境界で、backend上限を超えるnev+1を要求しない。対応するpairを枠内で認定できなければ成功packetを返さない。

各実返却pair／部分空間を、元作用と要求された随伴作用で再評価する。実作用に対する複素vectorの認定は実部・虚部の両作用を行い、その2回を実作用数へ課金する。backendの抽出済みvector配列は返却時にcopyして所有し、backendのscratchを借用しない。

Arnoldiのcluster要求は、固定版backendの実働Krylov基底Vと縮約Hをread-only監査adapterで取り出して処理する。実働次数mを使い、Hの複素Schur分解で要求した全k根を先頭へ交換してVで持ち上げ、候補Qを作る。未使用の確保幅ncvを処理へ含めない。従属し得るRitz固有vector列だけから独立なk列を復元できるとは仮定しない。列数を落とすQRや、任意の別vectorによる穴埋めも行わない。

QのGram上界が1未満であることを認定し、元作用から $$ H_{orig}=Q^*\mathcal A Q $$ を作る。その小行列のSchur変換後のQ/TをEIG-05と同じ部分空間残差と、要求された追加後退誤差で認定する。全k値がEIG-04のassociation政策で選択Ritz値と一対一対応することも確認する。shiftで候補を作った場合も元作用でこの処理を行う。大規模の全Aや完全Schurは作らない。候補部分空間の生成・数値資格が未達ならPrecisionUnresolved。入力寸法、通常未収束・件数不足、資源不足、callback失敗はEIG-12の分類を維持する。小問題の完全法cluster処理へ切り替える場合は呼出側が明示する。

## EIG-09：shift-invertは逆作用だけを認定して終わらない

一般作用のshift-invertは無次元座標で $$ \widehat\sigma=\sigma/s_\lambda $$、 $$ B_\sigma=(\mathcal A-\widehat\sigma I)^{-1} $$ とする。逆問題の値をθとすると $$ \lambda=\sigma+s_\lambda/\theta $$ である。一般化問題へこの式を無断転用しない。shiftと逆作用のsnapshot・分解世代・精度を要求へ含め、θの零・逆数や変換の非有限を拒否する。

逆作用callbackは同じshifted原系の線形求解を使い、線形認定と誤差上界を返す。近似逆作用を厳密な固定行列と偽らず、その作用誤差を予算へ加える。逆作用が要求精度を保証できない、または元作用がない場合は未認定である。元λへ戻した実返却値について $$ A\phi-\lambda\phi $$ を再評価する。θ側の残差が小さいだけでは受理しない。

左右を要求する場合は随伴のshiftを $$ \overline\sigma $$ とし、同じ原問題の左残差を検査する。shift近傍で収束しても、元作用の最大絶対値を全域で抽出したとは呼ばない。変換種別と抽出領域を結果に保存する。

## EIG-10：数値表現・結果の資格・物理被覆を型で分ける

候補生成backendの初期能力はBinary64である。精度政策は入力係数・作用入力／出力・保持値・返却値ごとに宣言する。未対応TwoWord入力や作用は計算前にPrecisionUnsupportedで拒否し、hiだけを渡さない。TwoWordを必要とする消費経路を含めて完成したとは主張しない。将来その経路を足す際も本書の原問題・実返却値の資格を守る。

保持値を返す要求と明示Binary64投影は別variantとし、結果generationを別にする。投影を要求した場合は固有値、左右vector、Q/Tの実投影値について、要求された残差・正規化・直交性・必要consumerを全件再評価する。保持値の合格を流用せず、未達はPrecisionUnresolved。現在の保持値が既にBinary64でも要求と世代の照合を省略しない。

成功packetには、FullSpectrumQualified、RequestedSubsetQualified、RequestedSubspaceQualifiedを別の型で記録する。未評価・未達の指標を有限のゼロ値で埋めない。個数、範囲、尺度、精度、元／返却generation、使用backend、実際の全予算・上界を所有する。必要なconsumerの評価器は実返却packetを受け、既存ErrorCertificateの観測・scope・予算bit列・世代・表現と一致した全件のAchievedだけを受理する。

物理的被覆は別の状態PhysicalCoverageUnresolvedを必ず持つ。このIM-00の入口からPhysicalCoverageQualifiedを構築できない。完全スペクトルでも固定された代数作用だけの範囲であり、物理モデルの空間・時間・後流解像度の資格ではない。物理被覆が必須の上位要求に対してはEigenCoverageUnresolvedを返し、安定、最小臨界値、無限大のフラッタ速度を生成しない。

## EIG-11：資源を事前に数え、全認定作用も課金する

Limitsは最大次数、返却数、ncv、Arnoldi backendへ渡す最大反復数と最大抽出試行数、元作用・随伴作用・逆作用の各最大呼出し数、owned payloadと外部backendのadmissionを必須とする。無制限の既定値はない。全寸法積・byte数をchecked 64bitで計算し、外部indexへの変換も検査する。

本契約で任意のbackend反復上限を公開するのはArnoldiである。実際の再始動回数を公開getterから得られない場合はNotExposedとし、任意の独立した再始動hard capを強制できると見せない。backendの報告counterの定義を保持し、computeへ渡した上限と実測のcallback数を同じ値と扱わない。完全法は固定版Eigenの停止規約30nと最大次数を記録する。自己共役一般化は内部固定上限、一般ComplexSchurは本実装の固定設定として30nを使い、設定値をchecked indexで検査する。公開されない実反復数を0や実測値として報告せず、backend未収束はinfoを検査して拒否する。これは壁時計時間のhard capではない。

| 経路 | 少なくとも見積もる同時領域 |
|---|---|
| 一般完全法 | 原A snapshot、backendの行列／Schur作業、複素Q/T、全固有値、左右／clusterを要求した作業と所有結果、認定scratch |
| 対称一般化 | 原K/M、M資格用CSCとCholesky、scaled候補のK/Mと分解、全vector、全結果と認定scratch |
| Arnoldi | n×ncvの実／複素basis、HessenbergとQR／再始動作業、返却最大nev+1、左右／cluster作業、shift分解と作用scratch、認定とcopyの同時領域 |

ResourcePlanの自前配列payloadと外部backend予測は別に判定する。見積りしかない場合にRequireProvenBoundを通さず、RSSや24/48GiBのhard capを認定したと書かない。準備済み分解でも今回の要求で再判定する。結果・証明・認定scratchは数値反復前に確保し、製品adapterの作用loopで新規確保しない。外部backendの反復中確保がある場合は隠さず別に計数・計測する。

backendの探索だけでなく、最終原残差、左右対応、consumer、shifted線形求解、再試行の呼出しも各予算へ数える。各callbackを呼ぶ前に枠を確認し、枠を超える一回を実行しない。枠切れはResourceExceededと既知のResourceUsage、通常の枠内未収束はEigenCoverageUnresolvedとする。backend反復数とadapter作用数を同じ単位へ混ぜない。

## EIG-12：失敗・所有・世代の保証を共通化する

新しい処理識別`ExtractEigen`をOperationへ足し、SolveLinearやSearchFlutterへ偽装しない。FailureCodeは既存の種類を使い、呼出しContextと最後の確定状態を保持する。

抽出入口の失敗型は`EigenFailure{core::Failure cause, AttemptDiagnostics diagnostics}`を所有する`std::expected`とする。内側callbackは既存`core::Result`を維持する。AttemptDiagnosticsは固定サイズのoptional値で、要求数・既知の収束数・backend反復数・抽出試行数・各種実作用数を保持し、未公開値を0にしない。確保失敗時にもこの失敗型を構築できるよう、文字列や動的配列を持たせない。共通Failureを勝手に拡張せず、causeの処理・code・既存Detailsを保存する。数値計算を始める前の失敗では未実行の項目を区別する。

| 条件 | 失敗 |
|---|---|
| shape・負予算・無効な抽出数・零初期vector・非対称K/M | InvalidInput |
| 非有限の入力／作用出力 | Nonfinite |
| 原Mの正定値未認定・必要vector basisの独立性未認定 | RankUnresolved |
| backend分解中の数値的破綻 | FactorizationFailed |
| backendの通常未収束・必要件数不足 | EigenCoverageUnresolved |
| 未対応表現・作用能力 | PrecisionUnsupported |
| 元残差・誤差包絡・投影・正規化・交換が未認定 | PrecisionUnresolved |
| 事前予算不足・作用数上限・確保失敗 | ResourceExceeded |
| 元問題・prepared・callbackの世代不一致 | StaleGeneration |
| 明示キャンセル | Cancelled |

作用callback自身のFailureは処理・診断を保持して返し、extract入口のContextを補う。cause.detailsのparent/last_acceptedは保持し、欠落している場合だけEIG-01の明示診断入力から補完する。渡されていない親IDを推測しない。外部例外は入口で値へ変換し、bad_allocはResourceExceeded、その他のbackend例外はFactorizationFailed、consumer例外はPrecisionUnresolvedとする。未収束時の内部候補を適格packetにせず、収束数・予算の既知診断を失敗の報告へ残す。

失敗・キャンセルは、入力・preparedの確定snapshot・親・既存結果を変更しない。可変Workspaceは呼出側所有で同時使用を拒否する。成功結果は移動・再利用・破棄されたpreparedやcallbackに依存せず有効。再利用キーは全stampと元係数実内容を照合する。並列Workspaceの返却が互いの結果を上書きしないことも試験する。

## EIG-13：独立参照と変異を要求へ対応させる

以下は実装後に行う試験の固定入力と判定であり、この文書だけで実行済みとは扱わない。小系はD=I、 $$ s_\lambda=m_0=1 $$ を基準とし、別に非一様な2の整数乗のDと単位変換を加えて同じ原問題の残差を照合する。一般化の2×2入力にはD=diag(3,0.2)、 $$ s_\lambda=7,m_0=3 $$ の非2冪ケースも加え、実返却値への変換丸めを検査する。原残差・正規化・直交性は1e−12以下。成分別後退誤差の追加要求はEIG-03の独立対照で別に検査する。固有値の集合照合は参照を重複消費しない完全一対一matchingで、各辺に $$ |\lambda-\lambda_{ref}|\leq10^{-12}\max(1,|\lambda_{ref}|) $$ を要求する。解析的fixtureの微小値の符号・零数は別の条件とする。PGの微小な数値値から厳密な符号や零数を証明しない。

| 独立入力／参照 | 確認する資格・負対照／変異 |
|---|---|
| 実6次：[[0,−2],[2,0]] ⊕ [[3,1],[0,−1]] ⊕ [0] ⊕ [−2⁻⁶⁰] | 根は±2i,3,−1,0,−2⁻⁶⁰。共役対・全6個、微小負の符号と相対1e−12、零1個。絶対許容差だけで零へclampする変異を検出 |
| 複素6次：diag(8,3,2,1,1/2,1/4)+UVᵀ、U=i e0、V=(1+i)e0+e1 | literalの先頭行は(7+i,i,0,…)、根は7+i,3,2,1,1/2,1/4。元literal作用・複素線形性・随伴一致を独立照合。V*への誤置換なら9+iとなる |
| 実Arnoldi：±5iの2×2 blockと対角0.2,0.3,0.4,0.5,0.6,0.7 | nev2/ncv8で共役の支配2根。一般完全法とも照合。作用を密行列で代用しない |
| M=[[2,2],[2,5]], K=[[-2,−2],[−2,10]] | Mの首座小行列式2,6、根−1,4、方向e0と(−1,1)。Euclid直交化・KのSPD強制・M⁻¹Kの対称化を検出 |
| M=diag(1,4,16), K=diag(−δ,0,16δ), δ=2⁻⁶⁰ | 根−δ,0,δ。負／零／正の数と微小非零の相対1e−12を独立検査 |
| A=[[1,1,1],[0,3,0],[0,0,1+δ]], δ=0または2⁻⁴⁰ | 非連続cluster{1,1+δ}の空間はspan(e0,e2)、projector diag(1,0,1)、制限行列[[1,1],[0,1+δ]]。Q列の直接一致を要求せず、projector・元部分空間残差・直交性を1e−12で照合。δ0で完全固有basisや有限pair感度を偽認定しない |
| 凍結PG／全自由後流P24と元SciPy固有値 | 全24値を重複度込みで照合、明示nev6/ncv20の部分抽出は最大絶対値の3共役対と一対一照合。完全／部分とも上記1e−12基準と元P残差を満たす。FD由来の固定Pの代数検証であり物理微分精度1e−12とは呼ばない |
| 非対称K/M、特異／不定M、非有限、shape／世代違い | 全成分検査・原M資格・個別K/M世代照合を確認。下三角だけの採用、世代の付け替えを検出 |
| 不変な元作用、途中失敗／例外／stamp変更／低位2⁻⁸⁰ | callback失敗・未対応精度・投影未達を区別。先行の必要consumer合格後に別consumer未達を置く。旧結果・親・入力をsnapshotで比較 |
| 作用上限ちょうど／一回不足、メモリ0、要求24/64をn24へ適用 | 認定用作用も課金し、追加呼出し前に拒否。clip、callback失敗後の正常結果、backend残差だけの合格、予算回数の初期化を検出 |
| n96の対角作用diag(96,…,1)、nev24/ncv64、明示非零初期vector | 既決の初期幅を小次数へ縮めず通し、上位96,…,73を一対一照合。合成容量の試験であり実物理mode被覆の資格ではない |

精密な元残差の参照はBoostの精密値またはcpp_rationalとし、製品残差関数を期待値にも使わない。PGの出典id・元hashとfixtureの来歴を維持し、Eigenを正解値の生成器へ替えない。shift試験は解析的な対角／非正規小系で元λの復元と元残差を確認し、逆問題のθを原固有値として返す変異を検出する。

## EIG-14：共通契約の後、独立した三実装を並行する

| 担当単位 | 先行条件・編集範囲 | 独立した成果物 |
|---|---|---|
| 共通型・認定・資源 | この設計の採用後。eig共通header、Operation、checked容量、認定helper・独立試験。CMake/追跡表は統合担当が編集 | 尺度・stamp・accuracy・所有packet・Failure境界。三backendが同じ資格を使う |
| 一般完全法／cluster | 共通型が固定・必要試験通過後。専用source/header/test | 全Schur、左右pair、明示cluster。小系とPG24の完全参照 |
| 対称一般化 | 共通型と既存SPD資格。専用source/header/test | 原M資格、負値・零を含む全実固有値、元K/M残差とM直交性 |
| action Arnoldi | 共通型と既存作用・線形精度契約。専用source/header/test | 一般実／複素作用、左右、shift adapter、全作用の課金と部分資格 |

三backendは互いの製品実装を待つ必要はない。完全法との照合が必要なArnoldi統合試験だけは完全法の固定結果を待つ。解析的な独立参照・負対照は先に用意する。共通headerやCMakeの同時変更を避け、別branchとbuild先を使う。

この追加契約で既存の数値結果を変更しない。再検証対象は実装される固有抽出とその依存範囲であり、既存の線形・Python・物理kernelを文書変更だけで再実行しない。IM00-11には本実装後の両OS製品検証と残る精密基準群が必要。IM-05には実際の同一平衡P/P*、枝・共通物理観測・部分空間角、ncv／seed／shift拡大、全状態解像度、左右感度と臨界区間の誤差台帳を渡す。IM-07/08の実務性能・実験受入を先取りしない。
