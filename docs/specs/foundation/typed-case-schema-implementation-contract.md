# 全schema validatorの所有型・資源・開始契約

2026-11-18。Issue [#118](https://github.com/YUKIKEDA/UVLAS/issues/118)、IO4-01、既存IM06-0の「全schema」を実装する担当者向けの追加契約である。[全schema契約](case-schema-contract.md)CS-05の内部表現・比較・資源・診断を具体化する。入力の全keyと物理的意味は[FD-01〜12](case-schema-fields.md)、初期摂動は[完全物理場追補](case-initial-physical-state-contract.md)が正本である。旧仕様・examples・証拠を編集しない。本書は実装開始の設計審査対象であり、製品validatorの合格報告ではない。

初回は所有型から失敗までを通読し、実装時は検査表と製品試験表を参照する。式評価、SI変換、snapshot、importer、prepare、run/resumeとCLI公開は既存IM06-0の依存段階である。schema成功だけでCLI `validate` の成功を返さない。

## 所有する型を固定する

入口はCS-05の `validate_case_schema(const InputDocument&, SchemaLimits, ResourcePlan&, SchemaDiagnostic*) noexcept -> core::Result<SchemaCase>` とする。diagnosticのnullは補助診断不要を表す。nullでも同じ検査・失敗分類を保つ。操作は既存 `Operation::PrepareModel`、親Model/Trialは存在せずoptional空。入力・plan・diagnosticの重なり、同一plan/slotの同時使用と再入を禁止する。異なる所有者の独立呼出しは並列実行できる。

`SchemaCase`は非公開factoryだけが構築できるmove-only不変型である。元InputDocumentを深く複製した所有snapshot、root payload、型別arena、ID索引、解決参照、出典表と資源内訳を持つ。元文書や呼出しworkspaceを破棄しても結果は存続する。snapshotは元bytes・SourceIdentity・node/member順・全decoded text・元decimal tokenを保持し、再parse/再serializeで作り直さない。成功前に複製へ全検査を終える。`InputDocument`内部の公開借用は同期処理内だけで使用する。

| 元の型表記 | 成功payloadの具体表現 |
|---|---|
| Bool、Tag、Index、Count、Id<T> | bool、閉じた型別enum、int64、uint64、namespace Tを持つuint64。Count→添字は範囲照合後のみ |
| Number / Q<d> | `ExactNumber{source_node, kind, sign, digits_slice, exponent10}`。Qは別の `Quantity<Dimension>` 包装。整数も元snapshotの出典とdecimal正規形を持ち、doubleを持たない |
| Rational / ExactQuantity | `ExactRational{numerator,denominator,source_node}`、または閉じた `variant<ExactNumber,ExactRational>`。分母の符号は正。約分を成功条件にしない |
| Text、Digest | 所有text arena内のoffset/byte_count、32byte配列。TextにNUL終端を要求せず、意味別の空/NUL条件を保持 |
| Vec<T,n>、Mat<T,r,c> | `array<T,n>`、`array<array<T,c>,r>`。各形状を照合してから書込み。混合次元の行列は成分ごとのDimensionを保持 |
| List<T> | 型別 `PlannedBuffer<T>` の `Slice<T>{offset,count}`。別型のsliceを交換できない。可変長の再帰vectorやmapをpayload内へ隠さない |
| 単一object | 辞書の型名に対応するnamed record。表の各keyを同名fieldへ一対一に写す。任意keyはoptional、noneも実variant |
| variant object | 辞書のタグ集合ごとの閉じたvariantと、各枝のnamed record。タグとpayloadの不一致を表現できない |
| Ref<T> / Target | `ResolvedRef<T>{user_id,typed_row,source_node}`、TargetはShell/Beam/Rigidのvariant。rigid行はMassのrigid枝だけ |
| FileRef、Source、出典 | role enum・path text・digest・SourceId、Source registry行への参照。全scalar/recordに元node ordinalを保持。元member pathは出典表から再構成 |

辞書の全object行にnamed recordを作る。例えば `TimeSeries<N,Dims>` はConstant/Polynomial/Pchip/Expression、`FullShellField<D>` はGlobalPolynomial/HalfCells、`ExprNode` はConstant/Variable/Add/Subtract/Multiply/Divide/Negate/IntegerPower/Sin/Cos/Exp/Log/Sqrtを持つ。`ReferenceField`は出力次元Dと成分数Nを、`TimeSeries`は成分ごとの次元配列を型引数または不変descriptorで持つ。一般JSON node、string→任意値辞書、未検証nodeを返す公開経路は設けない。内部descriptorを使う場合も全fieldの型・shape・variantと後述のdispatchは同じであり、未検証の汎用bagを成功型にしない。

beam stiffness/viscosityやsupportの混合行列は、FD-04/05の一般化ひずみ・座標から各(i,j)次元を決める。観測はFD-09の全行を用い、first_pk=9、reaction=6、stress他=6、vector=3、mode=1の長さと成分順を別descriptorへ固定する。interface=bothでもnを倍にしない。層入力のMandelと観測テンソル成分を同じenumにしない。

全17 namespaceはCS-04どおりで、Layer/Loadは部材内部も全域登録する。順序付きの `IdEntry<T>{id,row,source_node}` を別索引へコピーし、非再帰heapsortでid、同値なら元node順に整列する。探索は二分探索、同namespaceの重複は元走査で後に現れる定義を報告する。ID=0、UINT64_MAX、先頭zeroを許し、別namespaceの同値を許す。Sourceだけは同じid/path/hash/roleの反復宣言を同一行へ写し、別定義を拒否する。元payloadの入力順をsortで変更しない。

## 十進値と有理数を丸めず比較する

Numberの正規形は符号、先頭/末尾zeroを除いたdigit列c、十進指数eとする。零はsign=0、c=0、e=0。小数点以下桁数の減算と末尾zero除去をcheckedで行い、元tokenは別途不変に保持する。native整数の負端は符号付きabsでoverflowさせず、unsigned magnitudeとして変換する。負の零は値零、元字面はsnapshotに残る。

非零の絶対値はまず `digits(c)+e` をsigned checkedで比較し、同値なら左右のdigit列を右側virtual zeroで補って辞書順比較する。巨大な `10^e` を生成しない。例えば1e−400>0、1.000000000000000000000000000001>1、1e0=1.0をexactに判定する。順序逆転は負符号で行う。指数字句を読む前に最大桁数を検査し、超過をResourceExceededにする。

予算合計・時刻幅・刻みの整数商など算術が必要な場合は、符号付き10進digit配列で加減乗算・長除法を行う。shiftで生じる全桁数を先にchecked計算し、確保・workに余裕がある場合だけゼロ領域を具体化する。商を丸めず、(end−start)/stepの剰余が零かを検査する。有理数の比較は符号と分子×相手分母の交差積で行う。DecimalとRationalの混在は有限decimalのexact分子/分母へ写して同じ検査を使う。不要なgcdを要求しない。資源不足は未比較であり、不適合や等値を推測しない。

## 全辞書を同じ検査段階へ写す

CS-04の段階順を維持する。段階2ではrootからmember/arrayの元順に全treeのkey・kind・shape・variant・構造quotaを検査し、その全体が成功してから段階3のexact値/範囲/順序/排他へ進む。前方の値不正より後方の未知keyを先に報告する隣接段階の負例を置く。欠落keyの同順位は辞書のkey列順とする。ID登録・解決後に相互条件を評価し、同段階の最初の元nodeを失敗位置にする。既知variantへ他枝keyを混ぜる場合はInvalidInput、未対応tag/schemaはUnsupported。数値検査を先に行って未知keyを見逃さない。

| 辞書 | 必ず実装するschema検査 | 後続へ渡す未認定量 |
|---|---|---|
| FD-01 | root9必須/3任意、model8、単位enum、全null/欠落/空の区別、analysisとの全相互条件 | 単位倍率・角度のSI変換 |
| FD-02 | 全場/系列/軸/点/領域の型、次数と成分、重複powers、正length/error非負、時刻/イベント/period、使用座標/unused零 | 補間係数/C1/C2、法線射影、SO(3)、実境界所属 |
| FD-03 | 全幾何variant、degree/総次数≤4、重複項、width正、n=1/2/4、厚さ配列と全index、境界名一意、continuity相手/境界index | polygon正則性/被覆、全trace、外部source/importer |
| FD-04 | 全材料/断面/部材/層、正値とν範囲、層ζ端点/順序/無隙間被覆、全field/混合次元、2節点、owner記録 | SPD/PSD、自然計量、反り/断面適格性 |
| FD-05 | 全Joint/Support/Mass/Ledger/Load、mapping項/shape/index、ID、回転を持つtarget条件、replace自己/重複/循環、測度/次元 | 取付/拘束整合、仕事・接線、質量二重計上の物理判定 |
| FD-06 | 全面/attachment/sharing、chord/span厳密増加とindex、Flow/Wakeのnone適用、全正値/政策enum | 向き/incident、dual-quarter実幾何、gust直交/正規化 |
| FD-07 | 全analysis/branch/continuation/search、区間/刻みexact整数商、隣接/イベント、seed/coverage参照、coupled flutter必須true | 同一平衡枝、モード資格、探索/後流残差 |
| FD-08+追補 | initial/restart排他、全entity被覆、全frame/planes/wakeのshapeとID、birth時刻、下記完全場検査 | 初期基底/活動空間、姿勢/拘束/非貫通/最終状態 |
| FD-09 | 全観測・基底・比較・予算、layer owner、全成分shape、全12誤差群を各1回、代表/covers一意・自己/連鎖/循環拒否、配分のexact和≤親 | 参照の真正性、適用物理の実誤差、実受入予算認定 |
| FD-10 | 全政策/尺度/seed/cluster、正値・shape・相互条件、mixed初回全primary/拘束政策、dimension指数範囲 | 実layout/nev/ncv/seed展開、rank/抽出資格 |
| FD-11 | 全資源/出力/保存、非空/NULなしpath、参照ID、論理quotaのホストとのmin | RSS/24・48GiB受入、atomic保存、安全なpath |
| FD-12 | 全DAGのoperand<自index、outputs形状/範囲、全node到達、depth/edge/power、domain形状/順序/event/変数、全次元伝播 | 式値/微分/区間包絡、定義域正則性・物理被覆 |

参照の循環禁止は、入力が表す材料/座標継承・source依存、Mass.replace、誤差群代表へ適用する。現在のcaseに継承keyがなければ架空の辺を作らない。物理Joint、相互Continuity/AeroSharingの閉路を一般DAG検査で消さない。正当な共有の相互参照と禁止継承を別edge enumへ分類する。禁止グラフは元node順の三色非再帰DFSで全頂点を検査する。

DAGの深さはconstant/variableを1、他を1+最大operand深さとする。binaryの同じoperandはedge2本と数える。outputsから逆走し、未到達nodeを拒否する。次元は(L,M,T,angle)のint16の4要素で、加減は一致、乗除は和差、冪はexponent倍、negateは保持、sqrtは全偶数を半分、sin/cosはangle入力→無次元、exp/logは無次元→無次元。全中間の指数をcheckedで検査する。integer_power exponentはint32、abs(INT32_MIN)もunsignedで計算する。出力次元はconsumerの全成分と一致させる。domainの型・下上限・eventはschemaで検査するが、log(−1)など値の成立はprepareが区間包絡から拒否し、schema成功を評価資格にしない。

完全初期場は全Shell/Beam/rigid Massを各1回含み、R/P/Cと変位/速度を必須とする。VectorTermのpowers相異、R/P/C総次数13/21/2、全ExactQuantity、halfcellの重複/2成分indexを検査する。内部PGを直接参照できる場合は各index<2nも検査し、外部から実幾何が定まる場合は後続検査義務を保持する。正面積の全活動halfcell被覆・Cの区分双一次への所属・全共有後の独立係数はprepareで確認する。beam/rigidは両端/全frameとbody角速度を保持する。PointCheckはshellだけside_s/side_tを必須、s/t∈[0,1]とζ∈[−1,1]、beamのunused2成分exact零、全tolerance≥0を検査する。片側の向き・knot適合、点と場の物理整合と最終q/vはprepareであり、schemaで平均・補間・補完しない。同配列内の重複点はtargetと全exact座標・片側tagの同値で検査する。旧MaterialPointの他consumerをInitialPointへ置き換えない。Sourceのpath/hash/role一致やNumberの値同値から元SourceIdentityの同一性を推定しない。

## 全呼出しの資源を先に数える

SchemaLimitsの全項目は明示uint64で、0は零上限である。CS-05/FD-11のquota全てに加え、ホスト専用に `maximum_scratch_bytes, maximum_work, maximum_exact_digits, maximum_exponent_magnitude` を持つ。caseへ新keyは足さない。入力宣言で上限を広げず、既存quotaはホストとcaseのmin、追加4項目はホスト値を使う。元文書の7quotaも再照合する。巨大な宣言Countを使って配列を確保せず、実配列長とchecked照合する。

| 計数 | 単位と課金対象 |
|---|---|
| entities | CS-04のIDを定義する各record、Layer/Loadを含む。Source反復は一つの定義、宣言走査自体はworkへ課金 |
| references/scalars | 全Ref出現。Number/Rationalの分子分母/Bool/Tag/ID/Index/Countを各scalar、Numberはdigit桁数と別計数 |
| polynomial_terms/table_samples | 使用箇所を問わず全term、全tableの一つの時刻と対応value rowを1sample。coefficients列は各rowを1sample |
| expression_nodes/edges/depth/power | 全式合計node/operand出現数、最大深さ、最大絶対冪。pointsは後続評価点上限として保持し、schemaで評価した点は0 |
| typed/text/path | 成功所有snapshotの全可変payload+型arena+索引+出典+正規digits、複製text、最長JSON PointerのUTF-8 byte数。末尾NULは所有する場合だけ1byte加算 |

走査前に、文書所有dataの実node/member/edge/text件数を内部helperで取得し、最初のstackはN+1枠としてchecked計画する（深さを知るための未計上stackを作らない）。Nをホスト上限へ照合してから確保し、case側のquotaを読んだ後はminで再照合する。計画用の走査もworkを消費する。全サイズを数えた後、payload Pと同時scratch Sを計画し `P≤maximum_typed_bytes`、`S≤maximum_scratch_bytes`、およびResourcePlanの残容量を照合する。Pは元source/treeの深いコピーも含む。各PlannedBufferの `count*sizeof(T)` をABIごとに積み上げ、alignmentを必要とするbyte arenaはchecked align_upを加える。索引/DFS/重複項sort/decimal比較・算術/出典path構築の同時領域をSに含む。実装は再利用する固定scratch buffer群を一度確保し、検査中にvector成長やBoostの隠れた確保を行わない。桁算術は計画済みbase10 byte配列を使う。

容量式を次表で固定する。Bはsource bytes、Nは全node、Mは全member、Aはarray edge、Tは元text bytes、Dはcontainer深さ、Ntは型tのarena要素数（named recordのほかList内のNumber/Index/Vec/Rational/row等の全element型を含む）、RはRef出現数、EはID定義数、Gは全式node数、Jは全禁止graph edge数、Vはその頂点数とする。Lは全nodeの最大JSON Pointer bytes、Cは正規digits合計、Hは下記exact算術の最大temporary桁数である。いずれも実入力から数え、宣言上限を実件数としない。

| buffer | 固定容量とpayload bytes |
|---|---|
| snapshot | `sizeof(InputDocumentData)+B+N*sizeof(InputNode)+A*sizeof(ArrayEdge)+M*sizeof(ObjectMember)+T`。既存src/io/input_document_detail.hppの型を使い、構築専用keysは0容量で複製する。元ordinal/edge/member順を保つ |
| 型付き成功領域 | `Σt Nt*sizeof(Record_t)+C+U`。Uはsnapshotにない新しい所有text bytes。snapshot内textや同じarenaの共有sliceは実容量を一度だけ数え、重複保持する配列は別に計上する |
| 索引・解決・出典 | `Σnamespace Et*sizeof(IdEntry_t)+R*sizeof(ResolvedReference)+N*sizeof(SourceLocation)`。SourceLocationはparent ordinal、member/array ordinal、field enumを固定長で持つ |
| 段階走査 | `(N+1)*sizeof(ScanFrame)+N*sizeof(NodeTypeDescriptor)`。descriptorはphase2→3の渡し用scratch。元doc走査ではrootからDFSし、未用意のdescriptorを読まない |
| 重複項/点/ID sort | `max(M,E,K)*sizeof(SortEntry)`。Kは一つの重複禁止リストの最大件数。SortEntryは元行/ordinalの固定長対。heapsortは追加再帰stackなし |
| 禁止graph・式 | `V*sizeof(Color)+(V+1)*sizeof(DfsFrame)+J*sizeof(GraphEdge)+G*(sizeof(Dimension)+sizeof(uint64)+sizeof(bool))`。式のreachability stackはさらに`G*sizeof(Index)` |
| exact算術 | `6*(H+1)`byte。6本のbase10 digit buffer（2 operand、2交差積、商、剰余）を最大長H+1で一度確保。virtual-zero比較は追加配列0 |
| path | `(D+1)*sizeof(PathFrame)+(L+1)`byteのscratch。失敗結果の所有pathは別の`L`byte以下の専用bufferを計画する |

Pは上3行の総和、Sは後5行の総和とする。失敗pathは成功型と同時に返さないが、構築途中payloadと共存し得るため累積要求・同時要求へさらにLを予約する。P/Sをmax_typed/scratchと別に照合し、全体を既存ResourcePlan残容量へ照合する。上のcapacityはupper bound、実際のallocate countも記録する。配列のsizeof/alignは実ABIで求め、推定sizeofの固定値を両OSへ流用しない。

Hは行うexact演算の全候補から、値を計算する前に決める。digit長a,bの積/交差積はa+b、桁位置を揃えた加減はW+1、k個の非負和はW+digits(k)+1、長除法の剰余は除数桁数+1、商は被除数桁数以下を要求する。Wはoperandの最大上端桁位置から最小下端桁位置を引いたchecked幅である。decimal→有理化では分子/分母へ付くzeroを指数から先に数え、交差積のa/bへ含む。比較だけの巨大指数はWを具体化しない。最大要求Hがmaximum_exact_digitsを超えれば確保前に資源失敗とする。6本を越える同時temporaryを使う実装変更は開始条件の変更であり、未計上の追加確保を許さない。

累積ResourcePlanはscratchの破棄・再利用で返金しないため、計画には「累積確保 P+各scratch buffer確保の総和」と「最大同時生存 P+S」の両方を記録する。再確保するscratchがあればその全回数を累積へ加える。先に呼出しが使ったallocated/reserved bytesは消さない。全要求は確保前に既存allocateで予約する。内部型のconstructor/destructor/moveは非確保・非再入・非例外にする。

workは処理前に1単位を課金する。単位はnode/member/配列要素/edgeの訪問、比較するtext/digit1byte、桁の加減乗除の一つのdigit演算、sortの比較/交換、dimension1成分、ゼロ初期化/コピー1byteである。反復と全呼出しの合計をchecked加算し、上限に達した次の作用を行わずResourceExceededを返す。固定descriptor比較の走査も課金し、exact算術や二次の重複照合を無課金にしない。workは時間/RSSの上界とは呼ばない。処理中の経過予算をcoreの物理反復数へ偽装せず補助診断に保持する。

allocator metadata/array cookie/プロセス基礎メモリは既存ResourcePlanのpayload保証外である。今回のvalidatorは外部parser/動的多倍長を呼ばない。snapshot複製にもPlannedBufferを使い、parserを再実行しない。RSS上界を本契約から主張せず、上位の外部確保admissionと区別する。計画値を実測量へ入れない。checked加算/乗算/align/符号変換が表現不能ならResourceExceededで、飽和値や0を必要量として報告しない。

## 確保失敗でも診断を所有する

SchemaDiagnosticはreason enum、SourceIdentity、optional元node ordinal、期待type/field/variant enum、actual kind、optional typed ID、quota/limit/required/workと、所有JSON Pointerを持つ。CS-05の「所有path」を保ちつつ、確保前に失敗できるようpathは `Complete | Unavailable` の状態を必須とする。Unavailableは理由 `NotBuilt | PathQuota | Allocation` を持ち、短く切ったpathを完全と表示しない。不存在node/byte位置/IDを創作しない。欠落keyは親nodeと欠落field enum、循環は起点nodeと相手IDを保持する。

呼出し先頭でslotをresetする。固定部分の生成/複製/移動は確保不要。通常失敗では元順の親stackからJSON Pointerを作り、keyの~と/を~0/~1へ逃がす。path長を先に数え、quotaとworkを検査し、既存ResourcePlanに計上して正確な長さだけ確保する。診断pathは失敗返却時も存続し、借用やwhat文字列を含まない。診断path自体がquota/確保不足なら元の失敗code/reasonを保ちpathをUnavailableとする。元の原因が診断不足だけの場合はResourceExceeded/DiagnosticPathである。slot省略時もpath quotaは計画で検査し、成功/失敗が変わらないようにする。可変pathを持つ診断の所有コピーには明示ResourcePlanを使うcloneを用い、暗黙コピーで未課金確保しない。

| 原因 | code / reason |
|---|---|
| key/kind/shape/値/ID/参照/DAG/次元不適合 | InvalidInput / 個別enum |
| 未対応schema/tag | Unsupported / UnsupportedSchemaまたはUnsupportedTag |
| 注入された非有限不正document | Nonfinite / NonfiniteDocument |
| quota/work/scratch/overflow/確保失敗 | ResourceExceeded / 対象quota、OverflowまたはAllocation |
| その他の標準・非標準例外 | 既存guard_entryのInvalidInput / UnexpectedException |

最外側をguard_entryで囲み、内側Failureは詳細を保持する。例外変換後に診断がなければ固定fallbackだけを設定する。成功時slotは空。失敗時は部分SchemaCaseを返さず、入力、以前のSchemaCase、Model/Accepted/Trialを変更しない。自身の失敗予約だけを戻し、過去に成功したResourcePlan計上は残す。ファイルを開かず、GIL/callback/外部コードを呼ばない。

## 製品の全系列試験へ接続する

[旧設計例](case-schema-examples.json)144fragmentと5root、[初期場設計例](case-initial-physical-state-examples.json)4正常/22負例を、元hashとJSON Pointerで固定する。fragmentは実validatorの同じ内部型入口へ独立registry/consumer contextを渡す試験adapterを使い、全rootを無理に合成しない。公開APIに未検証fragment入口は増やさない。期待は手書き型/shape/値/参照で、製品から再生成しない。

旧 `Perturbation.physical` fragmentと旧PointValue経路だけはsupersededの拒否例である。旧rootにその旧型がある場合は元rootを負例として残し、全初期場追補を適用した新rootを別fixtureとして検査する。黙って旧入力を修復しない。4正常のschema成功は活動空間/SO(3)/非貫通の成功を意味しない。22負例は個別にschema拒否・prepare拒否・精度等の未解決へ分類し、prepare専用反例がschemaを通ることも正しい境界として試験する。

| 製品gate | 独立期待と負対照 |
|---|---|
| 全辞書 | 全行・全tag・全shape・全enum成員に正常と未知/欠落/別kind/他枝key/長さ±1。全root条件、全namespace/禁止辺/正常共有閉路 |
| exact | 0/負零/UINT64端/1e−400/1+1e−30/等字面、Rational符号/分母0、時刻剰余、12群予算和。float経由・pow生成の変異でassert失敗 |
| 式と初期状態 | 全opの型/次元/到達/depth、R/P/C全成分、両場・全halfcell・全frame・PointCheck両側。未使用node/次数/欠落/誤targetの変異 |
| 資源 | 全quotaちょうど/1不足、P/S/累積/全work、checked overflow、巨大指数、各確保点bad_alloc、非標準例外。失敗後の入力/既存結果bytesと累積計上 |
| 所有と診断 | 元doc破棄/移動、rvalue getter拒否、別source/同値別字面、全失敗reason/欠落位置/path不能、slot nullと有り、独立並列の決定性 |

全製品試験はIO4-01をタグに持つ。実装の変更をcompile/linkできる隔離変異と正常controlを走らせた場合だけ検出と数える。本設計の限定oracleは辞書inventoryとexact字句比較・段階分類の独立照合に限り、全製品validator/式評価/物理準備を実行したとはしない。

開始審査では本書とFD全行の対応、所有コピーの計画、全例の適用段階、診断fallback、exact/workの境界を独立レビューする。未解決のtyped variantや資源経路があれば、そのvalidator本体を開始可能としない。設計採用、実装適合、数値認定、実務受入を別に記録する。実装後に全系列gateが通るまでSchemaCase能力の完成を主張しない。
