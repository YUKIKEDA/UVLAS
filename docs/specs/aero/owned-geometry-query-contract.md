# 所有geometryとexact queryの実装境界

2026-11-17。Issue #95、IM03-2。[局所幾何契約](local-geometry-contract.md)2〜5/8節と[厳密述語](exact-triangle-predicates-contract.md)を具体化する。変更するのは内部の型、所有、接続の解決、位置式の整数入力経路であり、核・半径・接触・精度の政策は変更しない。

## 返すものと、返さない資格

内部型は `src/aero/geometry_snapshot.hpp` に置く。`OwnedGeometry` は参照/現在座標、全入力ID/接続、固定分割と散布、core/context、解決済みlookupを所有する。`ResolvedQuery` は元query ID、生成元実体、元の物質座標、位置式と正規支持を所有する。両者は移動専用で、getterは固定サイズの値をコピーして返し、内部配列の借用を公開しない。移動元は空と判定できる。

geometry payloadは `shared_ptr<const GeometryData>` で所有共有できる。queryは元のgeometry ownerを破棄した後も有効であり、queryごとに全点配列を複製しない。利用者が渡した配列・container・Plan・workspaceへの借用を残さない。snapshot完了までは入力の生存と不変を呼出側に要求する。

両型は `GlobalGeometryNotChecked` と、速度/写像微分が未検査である状態を持つ。正規支持は代数的な生成元の識別であり、`VerifiedIncidence`、`Detached`、自己除外maskではない。法線・sqrt・面積/長さ/距離包絡、ε=0自己/非incident接触政策、局所適格packetはA2が担当する。

## ID、入力配列、固定順

PointId、EdgeId、FacetId、PhysicalFaceId、BandId、QueryId、ReferenceMapIdは別の `core::Identifier<Tag>` とする。全uint64値を保持し、0を欠落にしない。幾何/精度/数値/モデル世代は既存の強い型を再利用し、参照写像・接続・核・半径の版も別tagで保持する。

| 入力 | 固定する内容 |
|---|---|
| points | 同長のPointId、参照/現在の各 `positions[point,xyz]`。各配置はBinary64または同形high/low、packed xyz、SI m |
| point binding | 各点の参照写像ID/版と、CanonicalまたはAliasOf{正規PointId}。alias先は同じテーブルのCanonicalに限り、未解決のchain/cycleを受けない |
| edges | EdgeIdと順序付き端点PointId。端点を正規所有者へ解決し、正規方向を所有者ID昇順にする。端点の原IDも保持する |
| physical faces/bands | 重複しないPhysicalFaceId/BandIdテーブルと面→帯の参照。点やFacetIdへ型変換しない |
| triangle sources | 安定FacetId、親PhysicalFaceId、頂点3PointId、正向き境界の3EdgeId。親面IDの反復はFacetIdの重複ではない |
| quad sources | 親PhysicalFaceId、前左/前右/後右/後左の4PointId、境界4EdgeId、内部対角EdgeId、生成する2FacetId。子IDを配列indexや切り詰めhashで生成しない |
| core/context | 共通SI半径の元語、核政策ID/版・半径版、モデル指紋/版、参照/現在幾何版、写像/接続/精度/数値版、有限な評価時刻。親/試行/最後の確定IDはoptional |

facetの3辺は、与えた3EdgeIdで解決し、対応する頂点対と照合する。同じ端点対の別EdgeIdをlookupで合併しない。各面の整数境界を構築して `D1 D2=0` を厳密検査する。四辺形は必ず(0,2,1)/(0,3,2)とし、対角辺の出現符号は逆になる。独立三角形の全3P1基底はそのまま保持し、四辺形には(L,L,R,R)から(L,R,L)/(L,R,R)への固定散布を付ける。μ値の平均・相殺は行わない。

aliasの値一致は、参照/現在それぞれの両語の厳密和で判定する。別のhigh/low分解や±0でも真値が等しければ一致する。写像ID/版は明示値を照合する。物理値の一致と、元語を保存するsnapshot内容identityを混同しない。速度/Jの一致はこの検査から推論しない。

## queryの型と正規支持

入力は `variant<FreePoint,SourceVertex,SourceEdgePoint,SourceFacePoint>`。FreePointだけが独立実座標を持つ。残りは実体IDと物質座標だけを受け、位置pを外から上書きできない。

| 入力 | exactな位置式と正規支持 |
|---|---|
| FreePoint{QueryId,words[xyz]} | 元の実座標を所有。物質支持なし |
| SourceVertex{QueryId,PointId} | 正規点の現在位置。支持はそのPointId |
| SourceEdgePoint{QueryId,EdgeId,t} | 正規端点方向の $$ p=(1-t)x_a+tx_b $$ 。 $$ 0\le t\le1 $$ を厳密検査。端点ならPoint、他はEdge |
| SourceFacePoint{QueryId,FacetId,u,v} | $$ p=(1-u-v)x_0+ux_1+vx_2 $$ 。 $$ u,v,1-u-v\ge0 $$ を厳密検査。正重みの個数からPoint/Edge/Faceへ正規化 |

支持がEdgeなら、その生成元facetが明示したEdgeIdを使う。元QueryId、生成元のEdgeId/FacetId、元物質語、元の位置式は支持正規化後も残す。任意の評価対象facetとの接触/適格性はA2が判定する。ここではqueryの生成元参照がsnapshot内に存在することを確認する。

`GeometryContext`の世代/時刻/coreと、factoryが実際の所有payloadから作った `SnapshotIdentity` をqueryへ結び付ける。identityは元語・typed ID・接続・写像・core・精度政策の正規byte列のSHA-256で、構造体paddingやアドレスをhashしない。callerのモデル指紋とは別の、局所内容の識別である。違う内容を同じ世代番号で差し替えて再利用しない。入力語表現が異なる同値aliasは許すが、その元表現を含むsnapshotの内容identityは異なり得る。呼出しごとのresource cap/取消/診断は内容hashへ混ぜず、各入口でその回のlimitsへ再admitする。不足はResourceExceeded、未証明strict要求はUnsupportedとし、内容/世代不一致のStaleGenerationへ読み替えない。

## 二語より長い位置を整数coreへ渡す

元語を $$ \delta=2^{-1074} $$ 単位の整数へ直接復号する。位置の両語和をX、物質座標の両語和をU/V/Tとすると、元座標は $$ x=\delta X $$ 、重みは $$ \beta=\delta B $$ 。暗黙の重みは整数で $$ B_0=2^{1074}-U-V $$ とする。範囲検査後に $$ P=\sum B_iX_i $$ を作れば、位置は既存coreと同じ $$ p=2^{-2148}P $$ で表せる。途中にも二語への変換を入れない。

| 段階 | magnitudeの厳密な保守上限 |
|---|---|
| 有限な元語2個の座標和X | $$ |X|<2^{2099} $$ |
| 範囲確認後の重みB | $$ 0\le B\le2^{1074} $$ |
| 一つの積BX | $$ |BX|<2^{3173} $$ |
| 3項を順に足す途中値 | 絶対和 $$ <3\,2^{3173}<2^{3175} $$ |
| 完成した位置P | 凸和により $$ |P|\le\max|X_i|2^{1074}<2^{3173} $$ |
| 元頂点のcore単位へのlift | Xを1074bit左shift、 $$ |X2^{1074}|<2^{3173} $$ |

完成値の凸和boundと途中の絶対和boundを分ける。したがって既存coreの差成分bound3176、式別最大31770bit、checked32768bitを維持できる。分母を持つ一般有理位置や、入れ子の任意深さの位置式へこの証明を拡張しない。

A0には値・単位・boundをprivate factoryで封じたexact入力経路を追加する。元word wrapperとquery経路が同じcoreへ入るが、leaseは入口で一度取得し、内部coreから再取得しない。一般の整数spanを公開し、尺度や範囲検査を迂回させない。

独立回帰には $$ x_1.y=2^{-1074},\ t=1/2+2^{-1074} $$ を含める。この辺上位置は $$ p.y=2^{-1075}+2^{-2148} $$ で、Binary64二語の和にできない。丸めた点を自己一致で再認定せず、正確な生成辺上という関係を既存述語で確認する。

## 処理I/F、作業域、資源

内部の論理I/Fは次の3つとする。`AeroResult<T>`は `expected<T,OwnedGeometryFailure>`。失敗型は共通 `core::Failure cause`、`optional<variant<PointId,EdgeId,FacetId,PhysicalFaceId,BandId,QueryId>> problem`、所有する `GeometryContext` を持つ。core共有型は変更しない。未解決shape等で局所contextを構築できない場合も、入口でcopyしたcontextと既知causeを保持する。いずれも外へ例外を出さない。

```cpp
prepare_owned_geometry(input, limits, workspace, plan, stop, details) -> AeroResult<OwnedGeometry>
resolve_query(geometry, definition, expected_context, limits, workspace, plan, stop, details)
    -> AeroResult<ResolvedQuery>
classify_resolved_query(geometry, query, target_facet, limits, workspace, result_plan, stop, details)
    -> AeroResult<ExactTriangleResult> // 数学的述語だけ。A2の適格packetではない
```

workspaceは実際のA0整数域と、次の7個の固定整数を所有する。重み3個、1語の復号域、座標両語和、積、位置和である。queryは元式をcompactに所有し、全queryについて32768bit位置配列を常時複製しない。exact位置はlease内のworkspaceで評価し、そのままcoreへ渡す。必要な述語結果は既存A0と同じ独立所有とする。

`ResolvedQuery`はgeometryへの所有共有参照と、Planで1個確保する `QueryStorage` を持つ。元座標/物質語は製品dataなので、固定サイズでも単なるmetadata overheadへ隠さない。snapshotの可変長payloadは先行factoryで計上済みであり、queryがその計上を複製したことにはしない。

ResourcePlanがdefault構築不能なstrong IDを確保するために、公開型へ意味のない0のdefault constructorを足してはならない。内部のtyped ownerでuint64配列を保持し、明示したIdの値getterで取り出す。入力配列と解決済み配列の両方を保存する場合は両方を計上する。

 $$ N_p,N_e,N_f,N_b,N_t,N_q $$ を点・辺・物理面・帯・独立三角形・四辺形の数、 $$ N_\triangle=N_t+2N_q $$ とする。位置の語数を配置ごとに $$ w_r,w_c\in\{1,2\} $$ とすると、座標payloadは $$ 24N_p(w_r+w_c) $$ byte。その他は次の全配列をcheckedに加算してから確保する。

| 所有payload | 個数と要素サイズ |
|---|---|
| point ID/正規owner index/写像ID/版 | $$ N_p S_{point} $$ |
| edge ID/原端点ID/正規端点index | $$ N_e S_{edge} $$ |
| physical face→band、band ID | $$ N_f S_{face}+8N_b $$ |
| 原triangle/quad接続と安定子ID | $$ N_t S_{tri}+N_q S_{quad} $$ |
| 展開facet、±incidence、生成元とscatter | $$ N_\triangle S_{facet} $$ |
| ID→ordinalの保持lookup | $$ 16(N_p+N_e+N_f+N_b+N_\triangle) $$ |

storageの固定フィールドは、pointがID/正規owner ordinal/写像ID/写像版の4個の64bit値と、元のCanonical/Alias指定を保持する16byteのvariant、edgeがID/原端点ID2個/正規端点ordinal2個、faceがID/元band ID、原triangleがFacetId/親面ID/点ID3個/辺ID3個、quadが親面ID/点ID4個/辺ID4個/対角EdgeId/子FacetId2個である。従ってpoint/edge/face/tri/quadのpayload要素はそれぞれ48/40/16/64/96byteとなる。

展開facetはFacetId、親面ordinal、点ordinal3個、辺ordinal3個、incidence符号3個と、`variant<IndependentTriangle{source_ordinal},QuadHalf{source_ordinal,half}>` の生成元を持つ。独立三角形に無意味なhalf欄を持たせない。 $$ S_{facet} $$ はこのstorageの実 `sizeof`、lookupはuint64 IDとint64 ordinalの16byteとし、その他のpaddingも `sizeof` へ含める。

`QueryStorage`のフィールドはraw uint64 query ID、所有する定義variant、正規支持variant。内部定義は `variant<FreeSingleWords,FreeDoubleWords,VertexId,EdgeIdAndParameter,FacetIdAndTwoParameters>` とし、座標3語/座標のhigh/low計6語を別alternativeへ直接保持する。公開FreePointの型は変更しない。各parameterは1語/2語のvariantである。支持は `variant<NoMaterialSupport,PointId,EdgeId,FacetId>`。raw storageでは各IDの数値を保持し、getterで正しいdomainへ明示構築する。query payloadは $$ S_Q=sizeof(QueryStorage) $$ の1個。固定query metadata/ownerは別として、可変長配列を追加しない。

```cpp
struct QueryArithmeticScratch {
    std::array<ExactInteger,3> weights;
    ExactInteger word, coordinate, product, sum;
};
```

workspaceのPlan要素はA0の `PredicateScratch` 1個と上記scratch1個だけで、payloadは $$ W_0+sizeof(QueryArithmeticScratch) $$ 。7整数の同一型配列なので追加分は $$ 7sizeof(ExactInteger) $$ 、現MSVC型では28728byteである。A0既存値W0=156008byteなら合計184736byte。値をhardcodeして計上せず各構成のsizeofからchecked計算する。ID lookupは保持済み配列をin-place sortし、別の可変長検索scratchを作らない。paddingをbyte出力へ流用しない。固定owner/control-block/context overheadは別のprocess overheadで、可変長配列をそこへ隠さない。

重み積の外部Boost確保forecastはA0と同じ7L limb、28672byteを維持する。根拠は同じ固定32768bit/et_offの `cpp_int/multiply.hpp`（SHA-256 `25413c62b71c24e9a4a42ea210ee206701c2a97b976960a2f9a39299905ab64c`）で、fixed setupの5sまたはsz+5s、s≤L、sz≤2Lによる。範囲検査済み重みは最大1075bitで、この積の小さい側は64bit limbで17、32bit limbで34であり、cutoff40未満でもある。位置を完成してから述語coreへ進み、一つのtop-level積の呼出しが完了して外部scratchが解放される前に別の積/coreを入れ子にしない。core内部のKaratsuba再帰は既存setupの共有bufferを使う。このscopeを超える同時確保やallocator overheadのproven boundを新たに主張しない。

数学bit容量、Planの累積payload、外部forecast、実heap/RSSを分け、proven bound不在のstrict要求はUnsupportedとする。snapshot/queryの要素数、接続/alias/述語の段階数にchecked上限と取消を設ける。根拠のない無制限値や飽和カウンタを作らない。

## 失敗順と最終開始審査

lease→shape/checked寸法→全外部語の有限性→ID/接続→context→半径/物質座標の実範囲→両配置の厳密非退化の順に検査する。同段はテーブル種別固定順、正規ID、入力ordinalの順で最初の原因を返す。alias先未解決/重複ID/不正接続/負半径/物質範囲はInvalidInput、非有限はNonfinite、版不一致はStaleGeneration、alias実位置不一致/厳密退化はGeometryInvalid。資源不足、取消、bad_alloc、未証明strict要求は既存A0と同じ分類とする。

固定順を段ごとに適用する。座標の有限性は参照/現在を合わせた論理点表から、非有限を持つ最小の(PointId,入力ordinal)を選ぶ。その後にtime、radiusの非有限を検査する。重複IDはpoint、edge、physical face、band、展開facetの順、接続はalias point、edge、physical face、独立triangle、quad展開子の順である。各表内はID/元ordinal順とし、quadは子FacetIdで選ぶためhalf0が常に先とは限らない。最後の非退化段は単一の展開facet表をFacetId/元ordinal順に検査し、各facetで参照→現在を検査する。所有配列・元ordinal・half・固定splitをこの順へ並べ替えず、既存lookupから検査対象を選ぶ。

点数を $$ N $$ とする。有限性の走査前に、設定済取消を先に確認し、checkedな $$ 6N $$ scalar分の残step予算を確認する。この確認自体では課金せず、各座標scalarの実訪問で1stepを数える。TwoWordもhigh/lowの検査を合わせて1scalarである。whole-pass予算拒否または走査途中の取消では、その時点で制御失敗を返し、未完了のcanonical不備を推測してproblemへ埋めない。完了した有限性走査だけから最小の非有限点を選ぶ。snapshotコピーの進捗とsemantic検査順は分け、推定stepや未訪問語をobserved量へ入れない。

新IDを既存 `core::ProblemId` の別domainへcastしない。決めた `OwnedGeometryFailure` に共通Failureを原形で所有し、補助problem/contextを値で保持する。内側の既知cause/resourcesを保持し、欠落だけを補完する。

実装開始前に、この型/配列順、Sの全フィールドと資源式、private exact入口のbound/lease、独立oracleをrootが審査する。oracleはcpp_rationalで元の両語と位置式を直接計算し、全query/零重み/全6頂点置換/alias同値別分解/±0/uint64端/固定split/整数境界/両配置を許容差0で比較する。false alias、split、incidence符号、weight-low、内容/世代を壊す隔離変異を選ぶ。A0refactorを含む場合は既存A0の影響3testと新規対象を、固定sourceからD/R各一度にまとめる。包絡・物理核・大域認定の完了は主張しない。

状態照合と内容hashには親AcceptedStateId、TrialId、直近AcceptedStateIdを含める。component、problem、resourcesなどの報告付加情報は照合から除外し、元情報の保存とは分ける。
