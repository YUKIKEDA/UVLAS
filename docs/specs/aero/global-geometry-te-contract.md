# 大域幾何・後縁所有の開始契約

2026-11-18。Issue #134、IM03-5、親IM03-0。[統合仕様](uvlm-integrated-specification.md)2節、[後流採用仕様](uvlm-wake-adopted-specification.md)5節の具体化である。接触を許す集合、物理alias、固定quad分割、後縁の政策は変更しない。既存の局所packet/P1結果は引き続きLocalOnlyである。

## 1. 処理と資格

内部処理を次の三つに分ける。型名は公開APIの予約ではない。

| 処理 | 入力・出力 |
|---|---|
| BuildGlobalTopology | OwnedGeometry、全PhysicalFaceIdのBound/Wake役割、edge binding、明示TE → 所有GlobalTopology |
| CheckGlobalGeometry | 同topologyと同一参照・接続の現在OwnedGeometry → 所有GlobalContactState |
| CheckAliasKinematics | 同geometry、共通global DOF配置、全点の速度・全疎Jacobian・明示時間項 → 所有AliasKinematicsState |

結果は非公開constructorで作り、元geometry/topologyのimmutable所有者と入力原本を保持する。元ownerを破棄しても有効。getterはlvalue限定のreadonly親または固定値を返す。配列の生存を呼出側に残さない。snapshot中だけ入力の不変を要求する。

GlobalContactStateが認定するのは、この離散時刻の有限三角形の接触と接続である。AliasKinematicsStateは与えた全行の一致を認定する。いずれも物理写像の正しさ、状態間経路の正則性、ANCF写像の向きの資格を生成しない。後続consumerは元の写像provider資格も要求する。現在法線と参照法線の内積による反転判定は禁止する。180度の適正剛体回転を拒否しない。写像資格のない入力は接触検査できても、循環・荷重・状態commitへ進めない。この区別をclosedなqualificationとして保持し、boolean一つへ潰さない。

contextはOwnedGeometryの全世代・政策・時刻とDOF配置版を保持する。context/hash一致だけで入力内容一致を代用せず、保持する元所有者または全内容を照合する。現在座標変更は接触/運動学資格を無効にする。参照、接続、役割、alias、TE変更はtopologyも無効にする。時刻/写像/DOF変更は運動学資格を無効にする。

## 2. 正規接続とTE

point aliasはA1の明示所有を使う。edge bindingもCanonicalまたはAliasOf{EdgeId}とし、別型IDを変換しない。alias先は同表のCanonicalに限り、chain/cycleを拒否する。正規辺の端点は正規PointIdの昇順。alias辺は同じ二端点でなければならず、元端点順と正規方向の符号を保持する。同じ端点を持つ別EdgeIdを暗黙合併しない。

各正規辺の整数頂点境界を $$ D_1=(-1,+1)^T $$ とする。保存済み幾何正境界incidenceを $$ D_2^G $$ とし、誘導核の時計回りsourceには $$ D_2^I=-D_2^G $$ を用いる。alias辺の出現は同じ正規端点順へ変換してから加算し、全facetで $$ D_1D_2^G=D_1D_2^I=0 $$ を整数で確認する。valenceが2を超えること自体では拒否せず、全ての面対について3節の許容集合を調べる。元PhysicalFaceIdと固定quad子、BandIdは残す。

TEは入力で明示した有向正規辺集合である。全TE辺はBound部分複体の境界に属し、Bound側のfacet出現が厳密に一つでなければならない。同じ物理quadの対角、二つのBound面の共有継ぎ目は不可。面IDの最小値や符号和の非零を使って一つ選ばない。付着Wake側の出現はBound出現数に数えない。重複TE、未知辺、零長辺を拒否する。

TE入力向きは付着wake前縁の**誘導sourceの出現向き**を指定する。正規辺基底に対する符号を $$ s_e\in\{-1,1\} $$ とする。wakeの幾何正境界符号とは逆である。Bound物理面の一定係数への整数散布を使い、

$$
T_{e,p}=-s_e(D_{2,b}^I)_{e,p},\qquad
D_{2,b}^I\Gamma+s_e(T\Gamma)=0
$$

を行ごとに確認する。左右帯の係数を平均しない。Tは疎な符号付き参照として所有し、循環値は本処理で計算しない。

TE節点の放出所有キーは型付きtuple `(CanonicalTEPointId, ReleaseEpoch)` とする。全uint64値を許し、hashや乗算で一整数へ詰めない。同一epochの全帯・全辺出現は一つの放出ownerを参照する。別epochは別owner。付着点のidentityはAttached、放出後はMaterial{キー}という別variantとし、現在同座標だから同一履歴としない。ここで作るのはimmutable放出計画である。同じ計画の再計算は同一順・同一内容。実際の一度だけのcommit、流体速度、移流履歴、J1時刻表はIM03-0の後流trial/commit子作業が担当する。

## 3. 全面対の接触

候補から自己facet対を除く。親面が同じ、頂点/辺が共通という理由で面対全体を除かない。固定quadの二子も検査する。epsilonは接触許可に使わない。

許容集合は明示aliasを解決した共通正規頂点と共通正規辺の実像だけである。三角形の閉集合交差を $$ K $$ とする。共通辺があれば $$ K $$ 全体がその閉線分内であることを要求する。辺がなく共通頂点だけなら、空集合またはその一頂点だけの交差を許す。交差線分の両端が別々の共通頂点に一致しても線分内部は許されない。二つの共通頂点があっても同じ正規EdgeIdがなければ共有辺とはしない。正面積の重なり・包含・一致面はID/親面にかかわらず拒否する。別IDで同座標の点/辺も無断接触となる。

製品の狭域処理は正規化しない整数外積/内積を使う。三角形Aの法線 $$ C_A=(a_1-a_0)\times(a_2-a_0) $$ とB頂点の符号距離分子 $$ h_i=C_A\cdot(b_i-a_0) $$ を計算する。全て厳密に同符号なら空。両法線が平行で平面が異なれば空。

非平行平面では、各三角形を他方平面で切る。零頂点を直接含め、厳密な異符号辺の切断点を

$$ P=\frac{h_u x_v-h_v x_u}{h_u-h_v} $$

という正分母の同次座標で保持する。重複をexactに除き、各切断は点または線分。$$ C_A\times C_B $$ の絶対最大非零成分の軸（同値はxyz順）で有理端点を比較し、二閉区間の共通端点を選ぶ。端点は必ず元の切断点から選び、浮動小数へ射影しない。

同一平面では法線の絶対最大軸を落とす。二次元の閉三角形内にある両方の元頂点と、9辺対の交差点を列挙する。平行共線辺の重なりは含まれる元端点で表す。最大24候補をexactに重複除去する。その凸包がKである。許容線分への包含なら全候補の線分所属、許容点だけなら全候補が同一点かつ許容点であることを確認すればよい。外積零・内積による閉線分所属を使い、許容差は0。

## 4. 整数容量と候補漏れの防止

A1と同じ $$ 2^{-2148} $$ 格子へ元high/lowを別々にdecodeする。各成分絶対値は $$ 2^{3173} $$ 未満。差分と符号加算余裕を含む基礎上界 $$ b=3176 $$ を使う。二進多項式のdegree d、展開項数mの絶対値上界を $$ 2^{db+\lceil\log_2m\rceil} $$ と評価する。

| 一時式 | d | m上限 | 必要bit上界（符号を別に1bit） |
|---|---:|---:|---:|
| 法線成分 | 2 | 2 | 6353 |
| 平面値h | 3 | 6 | 9531 |
| 切断点分子/分母 | 4/3 | 12/12 | 12708/9532 |
| 非平行切断点の座標比較 | 7 | 288 | 22241 |
| 共通線分への所属 | 5 | 48 | 15886 |
| 二次元切断点比較 | 5 | 32 | 15885 |

従って既存checked固定32768bit整数で足りる。正規化/除算不要の同次演算に限定し、未約分有理式を一般的に繰り返し乗算しない。容量超過を切捨てて成功にしない。degree/mが表を超えるコード変更は再導出する。

広域は全facetの閉AABBを使うbalanced BVHとする。箱端点は極値を与える元点ordinal6個を保持し、比較時だけexact座標を復号する。浮動小数AABBによる除外は禁止。構築は最大exact箱幅の軸、箱中心分子min+max、FacetIdの順で決定的median分割する。heapsortの比較/swapを計数し、opaqueな確保をしない。参照treeの葉順を保持して現在箱をrefitする。0/1 facetは空対集合。

自己node対は(L,L),(L,R),(R,R)、別node対は葉数の大きい側を分割する（同数はnode ordinal）。閉箱が交わる時だけ降り、各unordered葉対を一度だけ生成する。木高 $$ h=\lceil\log_2\max(F,1)\rceil $$、node数 $$ \max(2F-1,0) $$、深さ優先pair stackは $$ 1+4h $$ 枠で足りる。各分割で深さ和が増加し最大2h、追加保留は高々2枠だからである。stack順を固定して最初の失敗対も決定的にする。漏れがない理由は、交差点があれば両葉と全祖先の閉箱に含まれるためである。最悪候補数 $$ F(F-1)/2 $$ は残るが、全候補/全有理点を同時確保しない。

### 正のコア半径を持つ閉輪郭の点照会

`ContourContactIndex` は全体接触認定を終えた現在配置と参照木を所有し、元と同じ閉箱のrefitを一度保持する。利用する照会は自由コロケーション、または明示した非後縁の正規束縛辺の厳密な中点とする。中点の場合、消費側はその束縛辺の誘導を除く。照会は既存の元定義から $2^{-2148}$ 格子に展開し、所有する幾何と照合する。箱の端点に一致する場合も枝を訪問する。葉では面符号を厳密に比較し、零の場合は非零の法線成分に対応する投影面で3辺の向きを比較する。面符号は3次、投影した辺の向きは2次のままである。

束縛面の内部と四辺形の分割対角線、正規支持と一致する束縛辺中点を許す。接触する他の閉三角形は `GeometryInvalid` とし、照会IDと接触した最小FacetIdを返す。支持辺と同じ正規IDを持つ後流辺が中点に接触する場合や、適格な束縛辺を一つも特定できない場合は `Unsupported` とする。接触資格は点の包絡・最近実体・遠近判定・P1の積分精度から独立しており、それらを必要とする呼出しは既存の局所パケットを使う。

準備は既存のrefit上限を使い、`Node` 配列と認定済み幾何の所有記録を確保前に計上する。点照会は固定深さ優先stackを使用し、広域と狭域の訪問、照会元の厳密展開、辺IDの確認を一つのstep上限から消費する。照会ごとの確保は行わず、共有する木は読取り専用とする。個別の整数作業領域には既存の排他leaseと取消を適用する。上限内で全訪問を終えた場合のみ、接触資格と候補数を返す。

## 5. aliasの速度と全微分

global DOF配置は強いDofIdの昇順表、各成分の物理単位・基底・版を所有する。同じ長さだけの別配置を等しいとしない。全入力PointIdのxyz各行に明示時間項 $$ g_t $$ [m/s]、全 $$ g_q $$ のCSR疎行を指定する。列は共通DofId、係数の単位はm/各DOF単位。行内列は厳密昇順、重複列は禁止、厳密零は除く。係数・global速度・時間項は元1/2有限binary64語を保持する。符号付き零は値の比較で同一とし原本は保存する。

共通global速度vも所有し、各点で $$ v_a=g_{a,q}v+g_{a,t} $$ をexactな積和として**導出する**。速度結果は下記固定格子の符号付き整数で所有し、1/2語へ表現を制限しない。例えば $$ 1+2^{-200}+2^{-400} $$ も成功する。外部の速度を併記する場合は同じExactVelocity形式によるassertionとし、導出値と厳密一致しなければ拒否する。暗黙の許容差や任意速度への置換はない。CSR列のない成分は0。aliasは正規点と全列、時間項、導出速度、存在する速度assertionが厳密一致しなければ拒否する。正規ownerの導出値を共有し、外部丸め値は資格に使わない。位置と参照map ID/版もA1の検査に加え再照合する。一方向JVPだけを一致資格に使わない。各点が同一の物理写像を表すことは対応providerの資格に属し、任意の呼出側データを比較しただけでその資格を作らない。

運動学exact積和は二つのbinary64語の積を含むため格子 $$ 2^{-4296} $$ とする。各元量の整数幅は3174以下、積は6348bit以下、N項和に $$ \lceil\log_2\max(N,1)\rceil $$ と符号余裕を足す。Nはchecked uint64なので合計幅は6413bit以内。演算scratchは既存checked32768bit、持続値は8192bit容量に分ける。持続表現はsign、最下位非零uint32 limbの位置、長さ、packed limb offsetとlittle-endian非零帯配列（内部零limbは保持）とする。零はsign=0/長さ0、非零は両端limb非零の一意形。全32bit limbを256個以上使わず、actual sizeof(header)+4×保持limb数と整列を計上する。aliasは正規ownerを参照し値を複製しない。

ExactVelocityは認定snapshotから切り離さない。下流がTwoWordを必要とする時は、N4の明示投影で元exact値から中心と包含半径を作り、消費側の既存精度予算で認定する。二語だけを再び速度資格と呼ばず、予算不成立はPrecisionUnresolved。投影API/費用は当該consumerの開始契約で閉じる。単位の異なる列を数値だけで合併しない。参照写像の関数全体の同一性を、この一点の検査から主張しない。

## 6. 資源・失敗・並列

公開入口は既存AeroResultとcore FailureCodeを使い例外を漏らさない。無効ID/接続/TE/CSR/lease再入はInvalidInput、世代不一致はStaleGeneration、非有限はNonfinite、資源不足はResourceExceeded、キャンセルはCancelled、整数容量を超える要求はPrecisionUnsupportedへ対応する。具体operation、facet/point/edge、元detailsと親/試行contextを保持する。接触拒否は両FacetIdと許容集合を示すGeometryInvalid。失敗時は部分資格を返さず、親状態/Plan所有領域を変更しない。

全call計画はchecked countを先に行う。topology各表/lookup/CSR/元語/DOF/TE owner/T/所有controlblock、BVH node/並べ替えordinal、stack、最大24同次候補とexact scratch、診断、結果payloadを項別に数える。各配列はactual sizeof(T)×capacityとalignment paddingを計上し、所有共有済みpayloadは二重課金しない。sort比較上界は各長nに $$ 2n(\lceil\log_2\max(n,1)\rceil+1) $$ sift段、各段最大2比較/1swapとして計上する。

候補数Kは無確保BVH traversalのcount passで求める。node数Nとしてunordered自己を含むnode対の上界 $$ V=N(N+1)/2 $$ をchecked計算する。各node対のbox比較は最大6座標比較、stack pushは最大3、popは1。復号は各比較最大4元語なので、count passは $$ V(6c_{cmp}+24c_{decode}+3c_{push}+c_{pop}+c_{visit}) $$ を先にadmitする。各cは固定limb primitiveの静的opcode計数値であり隠れた任意定数ではない。各stepを消費し、その後K×狭域最大opcode数と運動学CSR積和数をcheckedで加え、全領域予約後に第二passを行う。未admitの候補を一つだけ先に評価しない。狭域opcode表は3節の固定展開から生成し、integer add/compare/multiply/shiftを固定limb走査数で課金する。事前counting型と実行型は同じ制御順を通し、実際stepが予約を超えたら失敗する。integer/候補scratchはWorkspace leaseで確保し、busyはInvalidInput。関数内heap確保やNに依存したstack配列を持たない。

継続利用する所有payload、同時生存scratch、全call累積workを分け、Planの規約に従う。Planの成功済みallocateは累積課金であり、一時領域を解放しても残量は戻らない。count pass用領域と結果領域を別々に確保する場合も両方を課金し、同じ予約領域を再利用する場合だけ二重計上しない。budget不足/overflow/cancelを全境界で検査し、catchして成功へ縮退しない。並列化は独立面対に限り、専用leaseを使う。失敗は最小の正規facet対へ決定的に縮約する。集約資格/TE順は1/2/4 threadでbyte同一を要求する。R-03の製品規模性能認定は本設計witnessの対象ではない。

## 7. 独立検証と開始判定

設計oracleはFractionによる6重心座標の非負実行可能多面体を使う。等式 $$ A\lambda-B\mu=0,\ \sum\lambda=\sum\mu=1 $$ の独立行rank rを求め、6列からr列を選ぶ全basic feasible解の像を列挙する。製品の平面切断/投影clipを再実装しない。literal交差頂点、許容集合判定、整数境界/TEのliteral期待値を固定する。許容差0。

開始判定の必要証拠は次の各行と対応する検証記録である。製品の合格を設計oracleの成功で代用しない。

- 横断/共面/平行/包含/一致/点/線分、共有頂点を越す交差、同端点別EdgeId、固定quad二子。各三角形の全6×6置換。閉AABB候補の全対参照との包含。
- binary64下位語による非零gap、$$ 2^{\pm900} $$ 尺度、巨大共通原点、正規方向反転と適正180度回転。alias速度/J/時間項の一成分不一致を個別拒否。
- D1D2、wake誘導符号を含むliteral T、逆向きalias辺、内部Bound継ぎ目/quad対角のTE拒否、三面共有辺での接触政策とTE政策の区別。同一頂点/epochの唯一owner、別epoch、uint64最大、Attached/Material区別。
- 開境界AABB、共有頂点の面対全skip、端点だけの許容点判定、誤ったT符号、alias基底未変換、内部TEの面ID選択、J一列無視、epoch無視の各変異が負対照を破ること。
- 製品実装時にはsourceword容量端、全資源equal/一単位不足、各確保故障、busy/cancel/非有限/別context/移動元、失敗後親不変、計画実step、1/2/4並列一致を別途実行する。これらは設計段階で実装済みとしない。

独立レビューで式・所有・資格境界・資源計画とoracleの独立性を確認し、新仕様がmergeされれば本契約の三処理を実装開始可能とする。写像provider資格、既知流れ、quad共有16変数、循環、荷重、後流更新、実務受入は既存IM03-0/IM01-0/IM06-0の各子作業が担当する。本契約だけでは全体開始審査を完了しない。
