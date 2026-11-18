# PG全原辞書・full Vの生成と資源契約

2026-11-17。Issue #108、IM01-3。[PG係数契約](pg-coefficient-contract.md) §3と§4の構成mapを内部処理へ具体化する。SH-D03/05/07、N4-02/08を適用する。方式、列空間、精度、物理範囲は変えない。原辞書・anchor補完・30必須場と並進のexact構成までを対象とし、活動rank、独立化、条件付け、物理N、全trace、多角形とPreparedPatchは親Issue17の後続に残す。本書の有限recipe・bit・生存領域を独立審査してから数値bodyを開始する。

## 所有する列と構成map

入力は実装済みPgGraph、明示limits、ResourcePlanである。結果はgraphの実内容snapshotとGraphContentIdentity、全背景セルのgとR/P/C、列由来、anchorのtrace/持上げ対応、必須場の構成mapを所有する。GraphContentIdentityは完全な物理参照identityの一部であり、Q/o/L_ref/厚さや活動域の資格を生成しない。callerの自己申告hashだけから内容を認定しない。

原列はr→d→paired→c。元の零列もそのID・由来を残し、独立列番号やCartesian展開番号と混同しない。scalar列iとCartesian成分alphaは別型とし、後者の展開はcheckedな3i+alpha。pairedのjは持上げの成分であり、Cartesian展開alphaではない。全成分はCellUVのcanonicalなexact多項式とし、係数微分と分母を含む物理場の微分を区別する。

O=2n、C=O²、r=(n+2)²+1、d=(2n+1)²+3、p=96n、c=(n+1)²とする。原列数S=r+d+p+cは122/246/530。anchor Rは高々n+3+48n、Pは高々2n+4+48nなので、返却列数上限は $$ F=S+99n+7=228,451,933 $$ 。実際にはt方向全区間のtraceがexactに非零の列だけを、R群→P群、各群内元列順で追加する。異なる元列から来る同じtraceは統合しない。計画用のcapacity上限と実際の論理列数を別に返し、limits不足時に列をclipしない。

s持上げ族はs=0で値と全一次jetが零である。実際、沿s微分はそれぞれ-f_s I0、-∫f_t k、I0であり、横微分もs=0では零、w/aの線形変換後も零なのでpのtraceも零となる。基本r/dのs=0 traceは先頭s-basisだけ、参照fは高々1列、法線分子は高々3列。t族は各R/Pで高々48n列であり、上記Fを得る。実装は全元traceをexactに検査し、予期しないs族非零を無視する短絡にはしない。

anchorは必ずc_s=0/u=-1の原traceを使い、t側の各片側区分を保持して全sセルへ延長する。pairedのR/Pがともに非零なら別の2列になる。Pは元のg(s,t)を分母とする。Cの延長は既存hat空間へ構成mapを持ち、新しいC列は加えない。T_a L_a=IとC包含を係数で確認し、V=S+L_a(T_aS)の由来と射影mapを保持する。clamped辞書だけを返さず、rankや独立性を証明したとは報告しない。

必須30列はX/Y/f、jごとのs全域曲率L_k(2s-1)の8場、法線X/Y/Z。t族も原辞書には全て含めるが、この30列の順へ混ぜない。定数並進は別mapで検査する。構成mapは原列の座標を保持し、anchor追加列の座標は零。係数を再構成してexact零を検査する。局所独立化のB/Sとは別の処理である。

## 有限な生成recipe

以下のA/M/S/D/I/E/Kは既存exact多項式の加減、積、scalar倍、明示軸微分、下端付き積分、点評価、cloneを表す。native定数factoryやmakeも省かず計画へ数える。 anchorのTはu=-1の係数制限であり、vも消す点評価Eではない。Tは各u指数の偶奇で係数を符号付けし、v指数へ集約するmakeを1回行う。負号の有理演算は高々入力項数回、borrowed入力列は253項以下のstack列、負係数のownerは同じ有限範囲へ保持する。t族は軸と元fの指数を入れ替え、最後に成分順とw/a変換を元へ戻す。以下は係数を再選択する新方式ではなく、採用直接積分式の実行順である。

| 段階 | 回数・有界な保持・呼出しの上限 |
|---|---|
| graph/cell | graphを1回cloneし、C個のcell多項式を構成。各cellでfの両chart微分からf_X/f_Y、U=1+f_X²、V=1+f_Y²、f_X f_Y、gを構成する。3M+4A+2D+4S以下を数え、1/wの有理演算は別leafとする。x/y/f/f_X/f_Y/U/V/f_Xf_Y/gの9多項式を保持 |
| univariate basis | 両軸・各O区間で、開いた二次C1/C0 B-splineをCox-de Boorの2段で構成する。分母零の枝は除算せずexact零。再帰node数は高々2O(6n+8)、各nodeはlocal affine係数を直接作り、2M+A、4個のnative分数と2個のaffine make以下。degree0のmake上限2O(3n+7)。一次hatは2O(n+1)個のaffine make |
| local powers | 各cellのs/t affineの0〜4乗を逐次生成。各軸でKを1回、定数makeを1回、Mを3回。一般の巨大powerループを使わない |
| curvature | L0=1、L1=u、L_(k+1)=((2k+1)u Lk-k L_(k-1))/(k+1)、k<7。各段1M+1A+3S以下。8個の係数を全族で再利用。t軸版も同じ有限漸化式を第2軸へ適用し、近似点補間をしない |
| moments | 各軸/origin/kでm=0..4。origin区間でsigma=c+h u、h=1/(4n)としてh sigma^m Lkを作り、下端-1からIを1回、上端+1をEで求める。前の区間は零、後ろは完了積分の定数。momentのactive計算は80O回以下、各1M+1S+1I+1E、最大次数12。全参照表は80O²個、各高々13項。m=5は不要で、f_t項のi+1<=4を使用する |
| basic scalar columns | 各cellでrのtensor積、dのtensor積とg倍、3法線分子、cのtensor積。Mは高々(n+2)²+2(2n+1)²+2+(n+1)²、Sは2、f/gのKを各1回 |
| paired J/v | 各vector seed/全cellでJ/J_t/J_stを元fの高々15/10/10項から合成する。35M+35S+35Aと係数j倍10回以下。vの3成分に3M+2S+4A以下。momentの借用に隠れた再積分はない |
| paired derivatives/p | 6D+9S+3A以下でv_X/v_Y jet、4M+4Aでt_X/t_Y、6M+4Aでpを作る。上記と合わせ、1 seed/cellは48M+46S+50A+6D以下。原scalar pairは3成分なのでvector seed数はp/3=32n、呼出し回数は(32n)C。native係数jと4n/wの有理演算は21回以下を別に数える |
| anchor | r+p本のR、d+p本のPを、各O個のt区間についてc_s=0のu=-1へ係数制限する。高々(r+d+2p)O回。各非零traceをCセルへell倍で延長するのでMは高々(F-S)C。重複も別列のまま保存する |
| required maps | X/YをGreville係数で、f/法線を補完列で、全域曲率を各originのlocal Legendreへ次数7のexact三角変換で展開する。30場と3並進の再構成は高々33F本のscalar倍/加算を各cellの各R/P/Cへ適用。C包含mapは高々c²項。これは局所RREFではない |

native/makeの照合は、実装の直接affine経路についても行う。geometryの各cellはcenterのnative2回とs/tのaffine make2回、s/t power cacheは各1 native constant/make・1K・3M。basisのdegree0 constant makeにはnative1回、hatは1rowあたりnative4回/affine make2回を数える。momentの1seed（軸/origin/k）はcoordinate/端点のnative4回とpower用native1回、affine/constant make2回・1K・3Mを持ち、5個のactive Iとendpoint Eに加えてcompleted constant make5回を持つ。下端付き積分をEに置き換えず、この追加storageも累積payloadへ数える。anchorのTは非空入力ごとにnative -1を1回、負号有理積は高々入力項数回、makeを1回持つ。

多項式の加算ループの先頭zero、空の入力、構造的に零のjetなどで呼出しを省けるが、計画値を越えて呼出しを追加しない。vector seedの作業用J/jetは次seedへ持ち越さず、最終v/pは出力へmoveする。moment carryはorigin区間で確定した上端値を後続へ渡すため、セルごとに積分定数を零に戻さない。

高次R/P全体をChartSTから一括substituteせず、最初のdegree4のfと有界なaffine/momentからlocal係数を作る。既存exactの指数展開上限を回避するためにadmissionを迂回することは許さない。既存call自体のB/U/byte forecastが親計画を越えないことを照合する。

## 次数・項・bitの上限

M(q)=(q+1)(q+2)/2を総次数qの項数上限とする。univariateならq+1。全てchecked64bitで計算し、a+bのIndex上限も格納前に確認する。

| 係数 | 総次数上限 |
|---|---|
| f、f_X/f_Y、g | 4、3、6 |
| I_m、v各成分 | m+k+1<=12、12 |
| t_X/t_Y、p各成分 | 14、20 |
| 原r/d/c | 4、10、2 |
| anchor R/P | 13、21 |

Aのraw上限はm1+m2、Mはm1*m2、S/Dはm、Iは2m、Eは保存する項数mと各軸powerを明示して既存plannerと照合する。canonical項数はraw上限とM(q)の小さい方。supportがunivariateやtensorで分かる場合はそれぞれq+1、(a+1)(b+1)も使える。fのaffine代入は既存graphと同じ最大raw129。次数を越える微分はaxis/order/tagの検査後、scaleやfactorialを作る前に零とする。

入力各分子/分母のbit上限をb、実際のf項数をmとする。bit計画では共通整数D0を構築せず、存在だけを使う。

$$
D_0=16!\,|\operatorname{num}(w)|\prod_{x\in\{x_0,w,a,f_1,\ldots,f_m\}}\operatorname{den}(x),\qquad
\operatorname{bits}(D_0)\le45+(m+4)b.
$$

native knot/scale/Legendreの各除数は16!の因子から作れる。Legendreの逆leading係数もk<=7の既知有理数に限定し、任意値の逆数を追加しない。wの逆は上のnum(w)で扱う。各計画nodeに、係数分母がD0^deltaを割るdeltaと、係数絶対値L1和が2^h以下であるhを持つ。入力/非整数native leafはdelta=1、整数leafは0。hは入力bitまたは既知nativeの絶対値から求め、零は別flagとする。

Aはdelta=max(delta1,delta2)、h=max(h1,h2)+1。Mはdelta/hとも和。native整数微分因子はdelta不変でhへbit因子を加える。次数dを下端±1で積分するIはdelta+1、h+1以下。Eの±1制限はdelta/hを増やさない。native scalar倍は同じ積規約で追う。したがって各返却係数の整数bitは $$ V\le h+\delta\,[45+(m+4)b]+1 $$ 。raw係数とGCD前の一時bitは、このcanonical値上限と別に、既存exactの有理四則・U(V+1)を適用してBを求める。Bを次段のcanonical Vと同一視して再帰的に膨らませない。

実plannerは上のD0を一律に累乗せず、native素因数2/3/5/7/11/13の指数と、入力由来の残因子の共通積の次数を別々に追う。加減で指数の最大、積で和、次数12以下の積分で既知除数の最小公倍数を使う。64bitを越える入力整数を未検査の既知素数へ分類しない。残因子はbit上限を共通積へ加算するため、入力の因数分解を必須としない。これはD0の上界を狭める整数帳簿であり、係数値の演算を近似しない。

X/Yのmapではx0/w/aとGreville係数の分母2nをすべて含め、x0=0でも他の座標係数を零と扱わない。再構成の途中和はtensorの全(n+2)²列または全8O区間曲率からの絶対和を数え、消去後の小さい値で途中和の上限を置換しない。

必須場の再構成にも同じ規約を適用する。map係数/和を省略して『raw係数と同じbit』とはしない。数学V、algebraic work B、guard込みplanned work、観測最大値を別に返す。元入力が極端でchecked上限を越える場合はResourceExceededであり、小項を落とさない。

## 一つの生成callのメモリ

所有payloadはraw/正規項/各storage、列/セル/由来/offset、mapと作業表の全allocateをResourcePlanへ累積計上する。前段の一時配列が解放されても累積課金をresetしない。planは各recipe nodeの実際のsizeofとraw/canonical上限から $$ P_{call}=P_{packet}+P_{tables}+\sum_j P_{exact,j} $$ を求め、最初の確保前にremaining_bytesへ比較する。既存graph clone/cell callのdata storageと累積payloadも加える。既存RawTermのpure sizeof accessorを使い、dummy予約をしない。

外部整数heapは累積payloadと別に、同時生存する値のproducer別capacityを数える。正規化がsizeだけを縮めるので、Vが小さくてもcapacityをVから作り直したとは仮定しない。各nodeで残るnum/denにはそのproducerが必要としたBからC(W)を付け、明示copyでcapacityを縮めた場合だけcopy sourceのsizeに対応する値へ更新する。入力は実capacityからH_inを得る。

fresh出力へ行う有理加減・積・除算について、交差積／GCD前の分子和は入力と出力のcanonical上限Vの最大を使った $$ 2\max V+1 $$ bitで包む。このstorage用上限と、既存exact入口が要求する $$ B=U(V+1) $$ およびword丸めguard込みの受入上限は別である。全leafで後者の検査を維持し、親forecastは実leafのexternal estimate以上を確保する。copy後に残るcapacityもC(W)で保持し、fresh値のL1には多項式の全項の絶対和とraw途中和を含める。

全出力項の保持集合の上限は、以下を別に数える。係数が零であることを確定していない段階のcapacityを計画外へ出さない。 n=4の密な3FC表現なら179136個のpoly slotとなる。実装のsparse record capacityは高々C(F+p)=84288個であり、下記原/anchor係数の粗い最大は9175872個になる。これは全fieldの総次数上限を一律に埋めたcapacity上限であり、実非零数や実確保数と呼ばない。零は既存exactの空表現を使い、各零のterm/storageをallocateしない。r/d/cの既知局所supportを保ち、gはcellごとに一度所有して列ごとに複製しない。非零係数と由来へのoffset/indexを所有する通常の疎格納を許し、列順や再構成を失わない。

- 原辞書： $$ C[(r-1)9+15+66d+(91+231)p+4c] $$ 個の有理係数。
- anchor追加： $$ C[26(n+3+48n)+42(2n+4+48n)] $$ 個。 traceはvだけの多項式なので、Rは高々13項、Pは21項、ellのs一次2項を掛ける。二変数の全M(13)/M(21)で埋めた値を実非零数としない。
- cell幾何cache：9poly、 $$ C[2+2+15+10+10+4(28)]=151C $$ 個。
- cellのs/t power 0〜4を保持するcache：各軸1+2+3+4+5項、合計30C個。
- univariate basis： $$ 2O[3(n+2)+3(2n+1)+2(n+1)] $$ 個。active再帰3段の保持は高々12n+10個を別に数える。
- moment表：80O² poly、高々1040O²個。endpoint有理数はactive Iの5個を別に数える。
- 必須/並進mapとC包含：高々33F+c²個。descriptorや±1だけで表すmapでも、そのnative storageはowned payloadへ数える。

vector seed作業の保持は、J/J_t/J_stの3poly、v_X/v_Yの6poly、t_X/t_Yの2poly、式のactive左右2polyとnative係数の2有理数。v/pは上記出力側のcapacityにも配賦するため、保守的な重複を許す。新しい同時生存変数を増やす場合はこの集合を更新する。native scalar factory・map構築・恒等式検査のleafは、呼出し前後の保持集合を同じ計画へ加える。

nested polyは既存契約の16+6U slotとG/D/K/resizeを使用する。親の保持capacityとactive leafのscratchを足した各時点の最大を外部forecastとし、全出力項へ最大scratch Bを一律掛けたり、各callのcapをresetして蓄積値を隠したりしない。既存leaf admissionの推定より親が小さくなる場合も拒否する。allocator管理領域/stack/OSの実RSS上界ではなくEstimateAdmissionの推定である。RequireProvenBoundは未証明Unsupportedを保つ。forecastの値をFailure.resourcesのobservedへ入れない。

native定数の字数、output容量、raw/canonical項、次数、V/B/guard、Plan/external cap、checked overflowを明示検査する。定数mapと巨大orderの零短絡でも入力のtag/内容/limits検査を省かない。全call preflightと各leafの実call/allocateの対応をprivate表・試験へ残し、根拠のない固定slot数へ置換しない。計画が保守的であることと、R-02の実graph/RSS受入は別である。 本体の最初のn=4実経路で工程別時間、実非零項/次数/bit、Plan累積payloadと同時生存/OS peakを早期に観測する。過大なら実測した保持や再生成の原因に対して表現を修正し、まだ測っていない高速性や48GiB資格を宣言しない。

## 公開境界・独立検証・開始条件

公開opaque所有型のfactoryはcore::ResultとOperation::PrepareModel、guard_entryで例外を値へ変換する。未知tag/shapeはInvalidInput、内容不一致はStaleGeneration、checked上限/確保失敗はResourceExceeded、取消はCancelled。内側の既知Failure詳細を保持する。readonly borrowはconst&だけを許しconst&&を削除、同じworkspaceの再入は拒否、別workerの作業を共有しない。失敗は元graph/旧辞書を変更せず、途中packetを成功型にしない。

製品の工程計測は係数構成・anchor補完・必須場/恒等式検査を分ける。全callの事前plan、metadata確保とgraph cloneはこの3工程の外側であり、3工程の和を全call時間と呼ばない。計測が無効なときはclockを読まない。同じ結果のcanonical係数/由来/mapがEnabled/Disabledで一致することを検査する。反復評価で辞書を再生成しない。短い準備工程の記録を実務性能や質量/条件付けの認定へ読み替えない。

独立参照は手書きliteral、別の有理式、版付きcanonical exact fixture。元R&D general_paired_polynomials.pyとpaired_exact_polynomials.pyの採用hashは保持するが、前者の既定clampedや後者の二次面限定を主経路にしない。原source/recordのpassed値を製品合格へコピーしない。全15曲面単項式、全n/全cell、両族全kと全3成分、carry、w/a変換、director恒等式と残る膜成分、R/P別anchorとC包含、全30map/並進、所有/取消/資源を許容差0で検査する。

制御byte同一の隔離control後、t族除去、carry零戻し、director符号/4n、当該cell左端への誤anchor、R/Pの一方省略、必須場map欠落などの意味ある変異を検出する。rawを保存し公開派生は匿名化/LF化後にhash固定する。必要対象だけD/R・staticを行い、旧代数/graph全suiteを変更理由なしで再実行しない。採用仕様と本書の式/順序/資源表、独立oracleが一致したことを別担当が確認した時点で数値bodyへ進む。
