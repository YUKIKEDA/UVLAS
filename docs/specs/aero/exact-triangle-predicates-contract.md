# 有限三角形の厳密述語の実装契約

2026-11-17。Issue #86、IM03-1。[局所幾何契約](local-geometry-contract.md)の有限三角形述語を、準備経路の内部部品として具体化する。採用された自己・接触の物理的意味や近接政策は変更しない。この部品は向き付き三角形と点の数学的関係を返し、面・辺の物理ID、共有辺の打消し、物質位置式、大域幾何の適否を認定しない。

## 入力と所有する結果

内部入口は `src/aero/exact_triangle_predicates.hpp` の `classify_triangle` とする。入力は SI m の $$ (x_0,x_1,x_2,p) $$ を各 xyz の順に並べた12成分。Binary64とTwoWordを別のvariantとし、TwoWordでは同じ長さのhigh/low両面を必須にする。各成分の真値は有限語の厳密和で、正規化やhighだけへの丸めを行わない。符号付き零は同じ数学的零である。

入力語は入口で固定長配列へ複写し、その後の算術・返却値は外部配列を参照しない。借用入力は複写完了まで生存し、不変でなければならない。結果は独立のResourcePlanへ確保したpayloadを所有し、入力、workspace、planの生存期間から独立する。結果への借用getterはlvalueだけに提供し、移動元は空を返す。Failureの既知contextは値で保持する。

この入口は12個の実座標を対象とする。共通単位 $$ q=2^{-2148} $$ を用いても、物質位置式 $$ p=\sum\beta_i x_i $$ の準備が完成したことにはならない。後続A1では、その厳密な位置式から同じ整数述語coreへ接続する経路を開始条件とし、二語に収まらない位置を先にBinary64/TwoWordへ丸めて本入口へ渡してはならない。

非退化結果は平面の上下、直交射影のFace/Edge/Vertex/Outside、有限面への厳密な距離二乗、全最小候補のmask、最長辺二乗とfar関係を所有する。退化は独立した数学的結果型であり、物理的に利用可能という成功認定ではない。導出長さ・法線の近似値はこの部品から返さない。

## 同じ整数単位で行う式

全入力語を $$ x=qX $$ の整数へ復号してから加算する。以下のベクトルは整数単位とし、 $$ e=X_1-X_0,\ f=X_2-X_0,\ r=P-X_0,\ c=e\times f,\ J^2=c\cdot c,\ h=c\cdot r $$ と置く。 $$ J^2=0 $$ だけを退化とする。平面側は $$ h $$ の厳密な符号である。

 $$ E=e\cdot e,\ F=f\cdot f,\ G=e\cdot f,\ u=e\cdot r,\ v=f\cdot r $$ とし、射影の重心分子は $$ N_1=Fu-Gv,\ N_2=Ev-Gu,\ N_0=J^2-N_1-N_2 $$ 。分母は正の $$ J^2 $$ で、各分子の負・零・正からOutside、開いたFace、各Edge、各Vertexを分ける。平面外でもこの射影分類を行う。

距離候補は次の閉集合から作る。各候補は $$ q^2N/D $$ 、 $$ N\ge0,\ D>0 $$ を表す。

| 候補 | 整数分子・分母 | 有効条件 |
|---|---|---|
| 頂点 $$ X_i $$ | $$ N=(P-X_i)\cdot(P-X_i),\ D=1 $$ | 常に |
| 辺 $$ a\to b $$ 、 $$ t=b-a,\ w=P-a,\ L=t\cdot t,\ s=t\cdot w $$ | $$ s\le0 $$ はa、 $$ s\ge L $$ はbの頂点距離。内点は $$ N=(w\cdot w)L-s^2,\ D=L $$ | 常に。非退化三角形の各辺は非零 |
| 面への直交射影 | $$ N=h^2,\ D=J^2 $$ | $$ N_0,N_1,N_2\ge0 $$ |

候補間の比較は交差積の符号による。返却最小値は全有効候補以下で、少なくとも一候補と等しい。距離零は $$ N=0,D=1 $$ へ正規化する。その他の分数へGCD簡約を要求しない。返却単位指数は常に−4296である。

maskのbit0〜2は頂点0〜2、bit3〜5は閉じた辺01/12/20、bit6は閉じた面射影を表す。同じ最近点を与える全候補を含むため、最小次元のsupportの一意ラベルとは異なる。例えば頂点0上の点では、頂点0・辺01・辺20・面が同時に最小となりmaskは105（1+8+32+64）である。後続の物理ID選択は別途、採用済みの最低次元と安定IDの政策に従う。

 $$ \ell^2=\max(E,F,(X_2-X_1)\cdot(X_2-X_1)) $$ とする。farは $$ N>16\ell^2D $$ の場合だけで、等号はBoundary、厳密な未満はNearと返す。Boundaryは近側に属する。ε、絶対長さの接触閾値、相対許容差は導入しない。

## 固定整数の容量根拠

Boost 1.92の `cpp_int_backend<32768,32768,signed_magnitude,checked,void>` と `et_off` を使用する。各語の最大指数を含む復号後、TwoWord和の絶対値は $$ 2^{3173} $$ 未満、差成分は保守的に $$ 2^b,\ b=3176 $$ 未満とする。次数d、符号を無視した単項式数Pの式は $$ d b+\lceil\log_2P\rceil $$ bitを保守容量とする。相殺による小ささへ依存しない。

| 式 | d | P | magnitude bit上限 |
|---|---:|---:|---:|
| dot $$ E,F,G,u,v,w\cdot w $$ | 2 | 3 | 6354 |
| 辺の端点判定 $$ s-L $$ | 2 | 6 | 6355 |
| crossの各成分 | 2 | 2 | 6353 |
| $$ J^2 $$ | 4 | 12 | 12708 |
| $$ h $$ | 3 | 6 | 9531 |
| $$ N_1,N_2 $$ | 4 | 18 | 12709 |
| $$ N_0 $$ | 4 | 48 | 12710 |
| 頂点N / D | 2 / 0 | 3 / 1 | 6354 / 1 |
| 辺N / D | 4 / 2 | 18 / 3 | 12709 / 6354 |
| 面N / D | 6 / 4 | 36 / 12 | 19062 / 12708 |
| 頂点−頂点比較 | 2 | 6 | 6355 |
| 頂点−辺比較 | 4 | 27 | 12709 |
| 頂点−面比較 | 6 | 72 | 19063 |
| 辺−辺比較 | 6 | 108 | 19063 |
| 辺−面比較 | 8 | 324 | 25417 |
| 面−面比較 | 10 | 864 | 31770 |
| far：頂点 $$ N-16\ell^2D $$ | 2 | 51 | 6358 |
| far：辺 | 4 | 162 | 12712 |
| far：面 | 6 | 612 | 19066 |

面−面比較は同時に一つの三角形で不要でも保守表へ含める。係数16はfarのPへ含めた。最大31770bitは32768bit未満である。この根拠は、入力が同じqへ復号され、式と項数が上表どおりである範囲に限る。新しい位置式や分母付き座標へ無条件に流用しない。

## 資源・中断・失敗

workspaceは実際に使用する固定整数配列・候補配列をResourcePlanへ計上し、結果payloadは別のResourcePlanへ計上する。workspaceのatomic leaseを入力読取りより先に取得し、同じworkspaceの再入・同時使用をInvalidInputで拒否する。独立workspace間の並行実行は妨げない。例外と早期失敗でもleaseを解放する。

固定整数であっても乗算の外部scratchが確保され得る。Boostのfixed Karatsuba setupは、最大operand limb数sに対し5s、結果側容量が不足する分岐では $$ sz+5s $$ limbを確保し、 $$ sz\le2L,s\le L,L=\lceil32768/w\rceil $$ から単一bufferのforecastを $$ 7L\,sizeof(limb)=28672 $$ byteとする。逐次評価・et_offの当部品に限定した事前見積りで、allocator overhead、全process/RSS、全外部確保の証明済み上界ではない。

専用factoryだけがこの見積りを生成し、proven boundはAbsentを返す。EstimateAdmissionでは28672byte以上を要求し、それ未満はResourceExceeded。RequireProvenBoundは指定値によらずUnsupportedで、算術・結果確保の前に拒否する。ResourcePlanのpayloadを外部確保実測へ読み替えない。準備時のCpp確保回数とprocess peakは別に記録し、物理反復kernelの無確保を主張しない。

製品と監査対象は `Boost::headers` を使う。必要なcpp_int/cpp_rationalはheader-onlyであり、広い `Boost::multiprecision` targetが持つ未監査のcompiled依存を監査対象へ追加しない。初回configureでBoost::containerの混入を監査gateが拒否したため、必要なheader依存だけへ絞った。compiled依存を監査済み扱いする例外は追加していない。

max_stepsは算術開始前に課金する名前付き準備段階の上限である。12成分復号、退化検査、重心・平面分類、3頂点候補、3辺候補、面候補準備、各有効候補の比較、far検査を各1段階と数える。退化なら13、非退化で面射影が無効なら28、有効なら29。0は無制限ではない。各段階で中断を検査しCancelledを返す。

検査順はlease→表現/shape→全語の非有限検査とsnapshot→中断→外部見積りadmission→結果payload確保→課金付き算術。shapeはInvalidInput、NaN/∞はNonfinite、中断はCancelled、payload/段階不足とbad_allocはResourceExceeded、予期しない算術例外はPrecisionUnresolved。operationはEvaluateAeroとし、入力に与えた既知FailureDetailsを上書きしない。失敗した呼出しで新たに課金したpayloadは既存ResourcePlanの累積計上規約に従う。

外部scratch根拠の固定source SHA256は、Boost `cpp_int/multiply.hpp` が `25413c62b71c24e9a4a42ea210ee206701c2a97b976960a2f9a39299905ab64c`、`cpp_int_config.hpp` が `ca7cfa2c2fcc6f6c5e4c6daaec365d60b49dbd3f85e0a1c5cce0208a031204bf`。前者のfixed setupを対象とし、variable backendのstack300limb分岐を固定backendへ転用しない。

## 独立な判定と開始条件

試験は入力語をcpp_rationalへ直接取り込み、別の最近点構成を使う。辺では有理な線分パラメータをclampし、面では射影点を原座標へ戻して向き付き辺で内外を判定する。各最近点から元の座標差を作り距離を計算し、製品の分子・分母の式を期待値生成へ転用しない。最小値が全候補以下かつ一候補と等しいこと、全最小mask、平面符号、射影feature、最長辺とfarを許容差0で照合する。

局所幾何契約のLG-F01/04/05/06/07/08/09/10/11/12の厳密部分、6頂点順序、符号付き零、全Binary64指数、最小subnormal、最大有限語、highだけでは潰れるlow、尺度・回転を対象にする。3-4-5三角形の内側射影を持つ高さ20は $$ \ell^2=25,d^2=400=16\ell^2 $$ 、直前/直後の語をNear/Farと固定する。T22の(3,3,0)は有限面距離二乗8で、無限平面距離0を使う変異を拒否する。

独立workspace、旧結果の所有、移動、同じworkspace拒否、全既知診断、寸法・非有限・正確な容量境界を検証する。隔離変異はlow破棄、辺clamp/面適格性/交差積/far等号等を壊し、controlとcompile/linkの成功後のassert失敗だけを検出に数える。固定source・依存版・実行条件と実際の未検証OSは検証記録へ分離する。

割当失敗は専用試験processの標準new/delete置換で注入する。workspace準備後の結果payload確保とBoost乗算scratch確保を個別に失敗させ、ResourceExceeded、既知診断、Planの予約解放/成功payload累積、入力・旧結果、lease解放と再実行を確認する。生成し得ない内部算術bugを起こすための製品hookは追加しない。その例外変換は既存core guardの基準試験に対応づける。

この部品の開始条件は上記の式、入力、返値、bit容量、資源区分、独立判定が閉じていることである。後続A1、有限渦核、共有辺、後流、荷重、複素Jacobian、全域幾何・TE、IM-03全体の開始・完了条件を代替しない。
