# PG多角形と体積・全外面の共通分割

Issue #117、親 #17。SH-D07、[主シェル仕様](main-shell-adopted-specification.md)2/5節、[PG係数契約](pg-coefficient-contract.md)2.3/5〜7節を、実装開始用に具体化する。採用した空間・厚さ・層・精度を変えない。読者はexact代数とPGグラフを使用する実装者である。

本処理が保証するのは、表現された領域のexactな交差・共通細分・活動根拠である。物理求積の精度、法線管の正則性、rank、条件付け、PreparedPatch、物理M/Kを認定しない。旧仕様はbyte保持し、本書を非矩形領域処理の追加規約として読む。

## 1. 入力の物理領域と補助線を分ける

無次元chartの単位正方形内にある、正面積の単純多角形Dを入力する。頂点は循環順、終点の反復なしで3点以上。向きは両方を受け、符号を記録して内部では反時計回りへ正規化する。連続同一点、自己接触、非隣接辺の交差・重複、隣接辺の折返しをGeometryInvalidとする。同一直線上の前進する中間点は許し、元頂点・辺の由来を残す。穴や複数の外周を1本の反復頂点列で表現しない。Dと1つの背景cellの共通部分は複数成分を持ち得る。

座標は既存ExactRational。有限十進とBinary64の変換入口は既存exact代数を使い、二進語の値を別の十進値へsnapしない。不確かさ区間だけの頂点から内外を選ばない。原domainが単位正方形外へ出る入力はInvalidInputとする。

全厚Hは既存Thicknessの頂点・三角形・正の頂点値による区分一次場、層は正規化厚さ座標ζの昇順区間である。各厚さ三角形は非退化、3つのindexは相異とし、D内の正面積部分をちょうど1枚で被覆する。D外の部分は積分せず、D内の重複・未被覆を拒否する。異なる三角形の実共有原子segment両端で、両側のaffine barycentric式から得たHをexactに照合する。端点一致はそのsegment上のaffine制限の同一性を意味し、面全体のH係数一致を要求しない。T-junctionや部分共辺も同じ原子分割へ揃える。不連続厚さを同一patch内の平均値へ変えない。厚さ段差は既存の別patch/接続処理へ渡す。層は隙間・重複なしで[−1,1]を覆い、界面両側を別traceとして保持する。

名付き境界は元domain上の順序付き線分列であり、各segmentの閉包がDの境界に含まれることを検査する。名前・元segment番号・向きを保持する。異なる名前が同じ物理辺を指すことは許すが、幾何辺を重複生成しない。同じ名前の列内での折返し・正長重複は拒否する。

| 入力型 | 固定する内容 |
|---|---|
| PartitionRequest | const PgGraph、domainのExactPoint2列、ThicknessMesh、LayerIntervals、NamedBoundary列、RefinementSegments。同期呼出し中だけ借用 |
| ExactPoint2 | chart(s,t)の2つのExactRational。物理xyzやCellUVと混同しない |
| ThicknessMesh | chart頂点、各3Indexの三角形、同長の正のH値。全三角形の向きを正規化し、元順を保持 |
| LayerIntervals | 各LayerIdとexact下/上ζ。材質自身を複製せず、消費側が同じLayerIdへ対応づける |
| RefinementSegments | 元区切りID/版とchart両端。材料/荷重等の既知区切りを追加する線分列。零長を拒否し、D外の部分も支持線は細分に使用するが物理外周を作らない |
| PartitionLimits | 全件数/bit/操作/payload、maximum_text_bytesと外部確保policy、stop_token。0を無制限としない |

Hの値はSI変換後のexactな長さであり、本処理はHを物理法線や体積測度へ変換しない。零・負のH、層間隔の零・負はInvalidInput。Hの同一性も結果の参照内容へ含める。

## 2. 元直線の共通細分から全交差を作る

方式は有理直線による単位正方形の共通細分である。凹多角形を凸と仮定したclipは行わない。次の直線を集める。

- domain全辺、厚さ三角形の全辺、名付き境界とRefinementSegmentsの全segmentの支持直線。
- 全入力segment端点を通るs一定・t一定の直線。境界の途中にある名付きsegmentの端点も失わないためである。
- 各方向の内部half-cell境界j/(2n)、j=1,…,2n−1。
- 単位正方形の4境界。

直線は整数同次座標(A,B,C)、式As+Bt+C=0で所有する。点の有理座標を(Ps Qt,Pt Qs,Qs Qt)へ変換し、2点の外積から作る。3係数のgcdを除き、最初の非零係数を正にしてcanonical化する。重複直線は1本にし、全由来を入力順で保持する。零長segmentの支持線を作らない。線は(A,B,C)の整数辞書式順とする。

全非平行線対の交点を同次外積から直接作る。正方形内の交点だけをcanonicalな有理点として所有し、(s,t)のexact辞書式順でVertexIndexを与える。交点を既存clip点から反復補間して分子分母を増大させない。各線対から同じ頂点への対応も残す。

正方形から始め、標準線順で各凸片を直線の両閉半平面へ分ける。交差点は、その片の元support lineと新しい直線の交点表から取る。線と一致する辺は二重の正面積片を生まず、接触点/線分は接触記録へ渡す。境界点を2片が共有しても面積の重複にはならない。片の頂点列は反時計回り、先頭は最小VertexIndex。頂点の重複と正長を持たない辺を除き、最後に全交点で辺を原子segmentへ分割する。面積はexact shoelaceで計算する。

各正面積凸片の、非共線な連続3頂点の重心を内点として選ぶ。元Dへの非零winding判定、厚さ三角形への符号判定、背景cellの厳密内外判定で所属を定める。補助直線で全元辺が切断されているため、内点の分類が片の内部全域に一定である。境界上に内点が来た場合は成功にせず、構成不変条件の失敗とする。各片にはcellがちょうど1つ、D内の片には厚さ三角形がちょうど1つ必要である。セル/厚さ所属を面積和だけで推測しない。

D内の片をcellごとにまとめる。原子辺を共有する片だけを連結し、点接触だけでは内部分を連結しない。同じcellの成分番号は構成片の最小番号順。ゼロ面積の接触は、元辺とcellのexactな閉集合交差の点/線分として別に所有する。活動フラグは正面積の有無からfactoryが導出し、callerから受け取らない。全4/16/64背景cellを欠落なく返し、空交差も根拠付きで保持する。

## 3. 体積・辺・上下・層は同じ片を参照する

凸片を先頭頂点からのfanで正向き三角形へ分割する。collinearなfan三角形は零測度として識別し、正面積三角形だけを積分domainへ返す。頂点列を曖昧なepsilonで単純化しない。三角形順は片順、fan位置順である。

| 出力record | 配列/向きと意味 |
|---|---|
| PlanarPiece | cell、成分、厚さ三角形index、CCW VertexIndex列、exact面積。無次元chart量 |
| IntegrationTriangle | pieceと3VertexIndex、正のchart二倍面積、元厚さ三角形。反復で再分割しない |
| AtomicEdge | 小さいVertexIndex→大きいVertexIndexを正向きとし、左右piece index、元domain/名付き境界由来の列を所有。隣接の不在は明示optional |
| VolumeCell | IntegrationTriangleとLayerId、ζの下/上。体積chartは三角形×閉区間、内部は重複なし |
| HorizontalFace | 同じIntegrationTriangleとζ、lower/upper、LayerId。ζ=−1/+1は外面、内部ζは両層のtraceを別recordにする |
| LateralFace | 物理外周AtomicEdge×層区間。domainのCCW向きの差(Δs,Δt)と、未正規化chart conormal(Δt,−Δs,0)を与える。物理法線/物理線長は後続 |
| EdgeTrace | 原子segment、元向き、名付き境界、隣接cell/厚さ区分/明示片側。補助辺に物理外周の資格を付けない |

左右pieceがともにD内なら、その辺は計算上の内部辺である。grid、厚さ、支持線の延長、fan対角線は、元domainの境界由来と外側不在の両条件が揃わない限りLateralFaceを持たない。box境界に一致するD辺も原domain由来から外周を認定する。名付き境界は物理外周集合の参照であり、同じ辺を二度体積境界へ足さない。

標準三角形の座標を(a,b)、a,b≥0、a+b≤1とし、chart点はv0+a(v1−v0)+b(v2−v0)。ζは層の下/上を用いるaffine写像である。後続は同じpieceのHとgraphから実法線厚さHζ/2、物理X0、Jacobianと正測度を作る。chart面積を物理面積として消費しない。上下の符号はζ増加方向に対して上が正、下が負。内部層界面の位置mapを共有し、材料/応力traceを平均しない。

本処理は求積点や重み、近似したN/Gを返さない。後続の求積はこの全domainを消費し、幾何/材料場に別の区切りが追加される場合は、同じfactoryへその明示線分を追加して再準備する。利用者が未分割の独立面求積を差し込む経路を作らない。

## 4. 所有・同一性・失敗を閉じる

入口は `prepare_integration_partition(const PartitionRequest&, const PartitionLimits&, ResourcePlan&) -> Result<IntegrationPartition>`。内部処理名であり公開ABIの固定ではない。結果はmove-only、元graphのsnapshot、元領域/厚さ/層/境界/追加区切り、全由来、全cell交差、片/頂点/辺/trace、資源報告を所有する。readonly viewはconst lvalueだけ。移動元の数値queryはInvalidInput。Workspaceの借用やcallerのarrayを保持しない。結果は専用factory以外で構築できない。

QualifiedActivityMapは同じ結果の全cell/交差証拠からprivate factoryで作る。旧UnitSquare mapを別領域へ読み替えず、graphと領域の両実内容を照合する。SHAだけで一致を決めない。canonical SHA-256は既存graphのstreamにdomain、Hのtriangles/値、層、名付き境界、追加区切り、直線/分割版を長さ付きで合成する。名付き境界のUTF-8名、RefinementSegmentsの区切りID/版とその全由来もcanonical streamへ含める。元入力順は意味のある由来として保持する。同値領域の異なる頂点開始点から同一digestを得ることは保証せず、同じ入力/版/環境の決定性を保証する。

検査順は、タグ/件数/限界/元graph有効性→入力shape/有理値/層/H→予測資源→domain/境界の単純性→直線/交点→片/所属/厚さ被覆→辺/trace/全cell被覆→所有結果の公開。資源拒否は未訪問の幾何不備を推測しない。全段の前と有限ループ内で取消を確認する。

InvalidInputはshape/index/未知tag/範囲/非正値、GeometryInvalidは自己交差・重複・厚さ欠落/不一致・境界不適格。StaleGenerationは元graph/参照内容不一致。ResourceExceededはbit/操作/byte/quota/checked overflow/bad_alloc、Cancelledは取消。外部の非有限は元変換入口でNonfinite。不足した幾何根拠や依存判定器の未解決はPG契約に従いRankUnresolvedとし、問題cellと活動根拠未認定を補う。正常な有限有理入力をtolで未解決へ落とさない。内側の詳細を上書きせず、例外は既存guard_entryで値へ変える。失敗後は入力/旧結果/親状態不変、成功済みResourcePlan累積計上を返金しない。

## 5. 有限件数と整数bitから呼出し全体を計画する

以下は保守上界の導出でありRSSの実測ではない。入力segment数S、元点数N、厚さ三角形数T、名付き境界segment数Nb、層数L、全入力有理scalar数Iを数える。Sはdomain辺・3T・Nb・RefinementSegments数の合計。元点Nは重複を含む入力の総数。追加の明示分割線分もS/Nへ含める。r=S+2N+4n−2、M=r+4を重複除去前の直線上限とする。

$$
C=1+r(r+1)/2,\quad V=CM,\quad P=M(M-1)/2,\quad U=2M^2,\quad K=4n^2.
$$

Cは正方形をr本で切る正面積片数、Vはclip中の片頂点出現数、Pは全直線対/候補交点数。各直線上の交点は他の直線の本数以下なので、全原子辺はM²以下、両側の辺出現/最終片頂点出現はU以下である。原子化passは候補交点を全て照合してよいが、その候補数VPを実出力数としない。各片のfan三角形もU以下、体積cellはUL、水平faceは2UL、側面はUL以下。全接触記録はK(P+M²)、由来対応はU(S+Nb)以下。全積/和をuint64でchecked評価し、Index/size_tへの変換前にも検査する。

入力有理数の分子/分母最大bitをbとし、b0=max(b,ceil(log2(2n))+2)、h=2b0、d=2h+1、k=2d+1と置く。入力点の同次成分はh、支持線はd、交点同次成分はk bit以内。交点の有理座標も分子/分母ともk以内。元線対から交点を作る理由は、この上限を繰返しclip数へ依存させないためである。

内点の3頂点重心はc=3k+3。shoelaceの各差は4k+1、全片/元領域の面積と面積和はq=U(4k+2)+2 bitで囲う。厚さ値は元三角形の有理barycentric式から評価する。pointのbitがjなら、orientation分母/分子の上限は2j+6b+5、元三角形のorientationは8b+5、各重みは2j+14b+10、H値3項の和は6j+45b+32である。共有原子辺の両端j=kで両側Hを照合し、affine trace全体の一致を証明する。H係数を浮動小数点でfitしない。

最大作業bitは、これらの式の各乗除/加減に[exact資源契約](pg-exact-algebra-resource-contract.md)の分子/分母式を適用した最大値Bである。例えば面積比較の交差積は2q、H比較は2(6k+45b+32)、point-line判定はd+2c+2まで含める。返却値のbitだけを中間上限にしない。既存のlimb丸め/guardを加えたplanned work bitもquotaと比較し、予測超過時に小片や小項を捨てない。

### 操作数は幾何bodyと有理演算を別に数える

NumericOpcodeは有理add/subtract/multiply/divide/compare/clone、整数gcd/exact-divide/符号正規化の各1演算である。初版は式をfresh非alias出力で逐次実行し、1opcodeにつき所有scalar storageを高々1個作る。操作上限を秒数の予測と呼ばない。bitは上記Bで別に制限する。

| body | opcode数の上界と内訳 |
|---|---|
| 2D orientation | 差4、積2、差1の7 |
| point-on-segment | orientation7、bounding-box比較4の11 |
| segment対の単純性 | orientation4回28、端点包含4回16、区間比較8の52 |
| 同次cross/正規化 | 積6、差3、gcd2、exact-divide3、符号3の17 |
| 線対交点の2座標構成 | 同次cross17、2つの有理正規化各gcd/divide/divideの6、計23 |
| point-line分類 | 積2、和2、符号比較1の5 |
| windingの1辺 | orientation7、y比較2、符号1の10 |
| 三角形内判定 | orientation3回21、符号3の24 |
| 1点のH評価 | orientation4回28、重み除算3、H積3、和2の36 |
| Hの両側trace | 両端×両側×36、等値比較2の146 |

同次crossの入力整数を有理p/1として処理しても、この数に内部gcdを再びopcodeとして足さない。nested整数scratchは外部forecastへ数える。整数原始演算を直接使う場合も上界を狭めず、正規化の生存slotを確認する。

全callのnumeric上限は次のRで計画する。各項は上表のbody数、同次入力点3積、定数/入力clone、centroid6、shoelace各辺4演算と面積/被覆集約、境界の両端包含22に対応する。

$$
\begin{aligned}
R={}&I+18+3N+3M+17S+23P+52S^2+5rV+7CM+6C\\
 &+C(10S+24T+4K)+4U+4N+4C+3K+3L+8\\
 &+11VP+146U+22U(S+Nb)+4K(P+U)+R_{sort}.
\end{aligned}
$$

直線・頂点sortは安定merge sort。件数xの比較上限をD(x)=x ceil(log2(max(1,x)))とし、Rsort=3D(M)+2D(P)+2P ceil(log2(max(1,P)))とする。後段のindexだけのsortは数値opcodeに含めず、同じDでtopology側へ追加する。canonical有理scalarを比較するときの内部乗算bitもBへ含める。

TopologyVisitはloopの1反復、recordの書込/比較、binary-searchの1比較を各1回とする。上限は各実arrayの全走査回数を保持する計画表から合算する。固定したloop形はS²、P、rV、CM、C(S+T+K)、VP、U²、K(U+P)、U(S+Nb+T)、U(L+1)である。配列構築は実record数、sortはD(件数)を加える。隠れた再帰や予定外の再走査を認めない。bodyの追加はRとbit/slot計画の同時更新を必須にする。実行はopcode/visitを呼ぶ前に残量を検査し、予定上界を超える場合は成功を公開しない。利用者はnumeric/topologyの別quotaを指定し、片方の未使用量を他方へ流用しない。

### 累積payloadと同時scratchを混同しない

所有payloadは実allocateごとのsizeof×容量を全call累積で計上する。計画配列は、元入力snapshot I scalar/全index、線M×3 scalarと由来、交点P×2 scalarと線対、片ping-pong 2V index、最終片U index、辺U/出現U、三角形U、体積UL、水平面2UL、側面UL、全cell K、接触K(P+M²)、由来U(S+Nb)、sort/queue各2max(U,P,M,C)である。各arrayはchecked上限で1回確保し、無計画なvector増設をしない。名付き境界名等の全UTF-8 byte合計をTextPayloadとしてmaximum_text_bytesとchecked比較し、所有byte配列1本と各offset/countを確保前に課金する。元文書/呼出しのstringを借用しない。コピー1byteをTopologyVisit 1回として全量を操作計画へ加える。recordの可変列はoffset/countと共有の所有arrayへ平坦化する。実sizeofはprivate純粋layout accessorで得る。

これに既存graph cloneの累積payload、NumericOpcode全R個分のRationalStorageを加える。snapshot/交点/線scalarのstorageがRと重複する保守計上は許すが、実Planには実allocateだけを課金する。scalar handleとscalar heap storageを両方数える。結果を破棄してもPlanは返金しない。Bad_allocやchecked overflowはResourceExceeded。

整数外部forecastは旧W/C/G/D/K/Rを継承する（この段落のK/Rは旧資源式の記号であり上の件数K/操作Rと区別する）。入力の実capacityをHin、保持有理slotをJ=I+18+3M+2P+C+K+6、逐次の幾何body scratchを16slotとして、Hin+limb_bytes[(2J+32)C+max(G,D,K)+R]を全callで1回受け入れる。6は面積集約2、centroid2、Hの両端比較2である。直線/交点/面積以外のrecordはindexで保持し、物理厚さは元入力を参照する所有indexとする。幾何bodyは中間を逐次解放し、16slotを超える式をそのまま入れない。

保持する値が途中で大capacityを持ったまま小さくなっても、Bから得る旧Cはそのfresh演算経路のcapacityを囲う。Hinは実capacityであり入力bitsだけから作らない。EstimateAdmissionは固定Boost/et_off/fresh非alias経路の推定受入であり、allocator overhead/stack/RSSのhard capではない。RequireProvenBoundは未証明なので確保前にUnsupported。各exact関数のcapをresetして保持配列を予算外へ出さない。

## 6. 独立照合と製品実装の開始条件

全比較の許容差は0。設計oracleはFractionによる共通細分と、別の垂直断面積分/手計算literalへ面積・一次momentを照合する。限定設計検査の合格を製品の交差判定器の完成へ読み替えない。製品では実返却頂点/所属/面/traceから以下を検査する。

| 系列 | 固定する照合と負例 |
|---|---|
| 全n=1/2/4 | 全4/16/64cell、空/正面積/辺点接触、標準順と活動漏れ/追加を拒否 |
| 凹形/細片 | 1cell内で離れた複数成分を持つU形、非格子有理交点、正面積2^−120、同一直線上の境界点。面積と両一次momentを独立照合 |
| 厚さ/層 | cellを横切る厚さ三角形、複数層、共有traceのH一致、未被覆/重複/不一致を拒否。全体積cellと面の共通頂点/片indexを照合 |
| 境界 | 逆向き入力、名付き部分辺、boxと一致、内部補助線。内部辺を外周へ出す変異、面だけ未分割、片重複/欠落を検出 |
| 幾何不適格 | bow-tie、自己接触、折返し、零辺、範囲外、退化厚さ三角形、名付き境界のdomain内部通過 |
| 所有/資源 | 入力破棄/変更、clone/move、同shape別参照、全quota同値/1不足、件数積overflow、返却/Boost確保失敗、取消、独立Plan並列とbyte決定性、旧入力/結果不変 |

実装開始には、具体型のshape/所有/失敗、各演算列のbit伝播・累積storage・同時slot、限定oracleの独立性、全系列の負例と変異計画を別担当がレビューする。実装Issueはこの契約の全能力をDoDへ接続し、製品の交差判定器を独立試験で閉じた後に、多角形rankと物理求積の後続を開始する。方式候補の探索や任意の切断配置の実務精度認定へ話を広げない。
