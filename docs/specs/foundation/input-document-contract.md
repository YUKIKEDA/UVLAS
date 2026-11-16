# 入力文書kernelの内部契約

2026-11-16。Issue #22の最初の依存閉包。対象はUTF-8入力を所有する文書へ移す処理である。本書は[IO4-01/02](io-execution-adoption-specification.md)の内部具体化であり、物理schema、Model、prepare、runの成功条件を変更しない。独立した技術・日本語レビューを終えた契約候補であり、採用判断はPRレビューに委ねる。

## 正本と今回の境界

IO4のケース形式・意味、[N4の型・所有・失敗](numerical-foundation-adoption-specification.md)、[失敗契約](failure-certificate-contract.md)、[寸法・資源契約](dimensions-resource-contract.md)を継承する。物理値の意味は構造・空力・J1の採用仕様と追加規約が持つ。本書はこれらと矛盾する例外を作らない。[IM-00終了審査](../../plans/im00-completion-review.md)の再利用範囲内で、`io`の文書読取り部品を閉じる。

成功は「このbyte列が本書の文書条件とroot envelope条件を満たし、値を失わず所有できた」である。全nested schemaの検査、参照先ファイルの解決、SI化、物理Model構築、幾何/rank/材料/精度の認定は行わない。結果型は`InputDocument`とし、`ValidatedCase`、`PreparedModel`、`QualifiedTrial`へ暗黙変換できない。`valid_case`や`ready_to_run`の真偽値を返さない。CLIの`validate`正常終了やPythonの物理case構築成功へ、この結果だけを接続してはならない。

nlohmann/json 3.12.0を内部の構文解析に使い、独自のJSON parser、別依存、parser型の公開を追加しない。採用済みのnlohmannを使うことと、case全体が検査済みであることは別である。以下の型名は内部契約の名前であり、公開ABIを固定しない。

## 入力と所有する結果

入口は同期処理 `read_input_document(bytes, source_id, limits, resource_plan, diagnostic)`。`bytes`は長さ付きの`span<const byte>`であり、末尾NULを要求しない。`source_id`は呼出側が割り当てる文書入力専用uint64識別子で、物理ID、ModelFingerprint、配列indexと別の強い型にする。0も有効である。同一IDの文書を同内容と見なさず、内容digestも必ず照合する。

この入口はファイル名を開かず、URL解決も行わない。ファイルadapterとPython adapterは読取り済みの同じbyte列を渡す。相対pathや含まれる文字列をこの入口が自動で解決しない。参照解決の成果物と元入力snapshotの接続は別の依存閉包で定める。

成功結果は次を所有する。所有者の寿命は入力buffer、parser、ResourcePlanの寿命から独立する。

| 所有値 | shape・意味 |
|---|---|
| SourceSnapshot | 元bytesを長さnで正確に保持。書換え不可。UTF-8正規化、改行変換、末尾NUL追加をhash対象へ混入しない |
| SourceIdentity | source_id、n、元bytesのSHA-256。digestは32byte。hash対象はファイルpathや再serializeしたJSONではない |
| Value tree | rootを含む全値node、配列の入力順、object memberの入力順、復号済みkey/string、整数値または十進tokenを所有 |
| Envelope state | schemaが`uvlas.case/1`で必須root keyが揃い、root未知keyがないことだけを表す状態。物理case認定を含まない |

nodeはNull、Boolean、String、SignedInteger、UnsignedInteger、DecimalNumber、Array、Objectの区別を持つ。Arrayは子node ordinalの列、Objectは復号keyと子node ordinalの列を持つ。node ordinalは文書内の0始まりint64であり、物理IDではない。nullと欠落、空objectと空array、整数と小数を同じ状態にしない。rootはobjectとする。

利用者へnlohmann型、parserのiterator、第三者例外を返さない。借用getterは`const &`に限定し、一時所有者からの呼出しを禁止する。treeの移動・破棄を越えるviewを許さず、Pythonへviewを返す将来のadapterは所有者を保持する。文書は不変であり、同じinstanceを書き換えてhashだけを再利用するAPIを設けない。

## JSONの受入条件を上書き前に検査する

JSONの構文・string escape・UTF-8検査はnlohmannの厳格モードで行う。全inputを消費し、comment、trailing comma、trailing tokenを拒否する。空bytes、先頭のUTF-8 BOM、UTF-16/32、壊れたUTF-8、不正surrogate escape、string内の生NULを拒否する。JSON stringの`\u0000`は正しい復号値として保持するが、将来のpath型への変換では別に拒否する。JSON外側の空白はJSON標準のspace/tab/CR/LFのみとする。

object開始ごとにそのobject固有の復号key集合を持ち、次の値で上書きする前に重複を拒否する。`"schema"`と`"\u0073chema"`は重複である。大文字小文字やUnicode正規化は同一視しない。親子objectのkey集合は別々に照合する。SAX等で値を所有treeへ構築し、構文不正をskipした部分treeを成功として返さない。

rootの許可keyは必須の`schema,units,model,flow,analysis,observations,numerics,resources,output`と、任意の`initial,restart,provenance`だけとする。schema値はstringの`uvlas.case/1`に厳密一致させる。未知root key、欠けた必須key、schemaの型不一致を拒否する。未知schema文字列はUnsupportedであり、別のschemaへ読み替えない。

root以外の値はこのkernelでは所有するだけである。たとえばunitsの内部key、modelの8集合や各集合のshape、analysisのvariant、材料・幾何・観測の未知keyを検査済みとはしない。後続の全schema validatorは、ここで保持したtreeの全nested key/shape/型を検査してから別の結果型を構築する。未検査のnodeを検査済みnodeに混ぜたpartial-success packetを作らない。

## 数値tokenとIDを丸めずに渡す

数値の構文判定はparserに任せる。`true`、`"1"`、`null`を数値として扱わない。整数eventは値をint64またはuint64で厳密保持し、doubleを経由させない。DecimalNumberはnlohmann SAXの`number_float`に渡される元token文字列をそのまま所有する。指数や小数点のない大整数でもparserがfloat eventへ送った場合はDecimalNumberに保持し、native integerへの切詰めを行わない。

native整数の元の字面はSourceSnapshotに保存される。treeの値が0であることから元tokenが`0`だったと推定しない。DecimalNumberのtokenは復号・書式変更しない。`1.0`と`1e0`の字面を保持し、nodeを整数1へ暗黙変換しない。tokenによってparserが生成した近似doubleは、物理量の認定値として公開しない。

十進tokenの値は、符号をs、数字列の非負整数をc、小数点以下の桁数をf、指数部をeとして $$ x=s\,c\,10^{e-f} $$ で定義する。この式は独立oracleと後続の型変換の意味であり、任意精度数を全tokenに対して製品内で構築せよという要求ではない。元字面の精度をこのkernelがbinary64へ消してよい理由にはしない。

`NaN`/`Infinity`はJSON tokenとしてInvalidInputである。引用された同文字列はStringとして保持する。parserが数値overflowにより有限のfloat eventを生成できないtokenはNonfiniteで拒否する。underflowでparserの近似値が0になっても元十進tokenを保持し、元値が0だったとは主張しない。物理量へのbinary64/TwoWord変換・表現不能・丸め予算の認定は、対象量の契約を持つ後続処理が行う。ここに暗黙の`as_double`は設けない。

後続が使うID文字列とindexの厳格変換は、この部品の独立関数として閉じる。IDはASCII digitの非空文字列で、符号・空白・小数点・指数・区切り文字を認めない。先頭zeroは認め、数値として同じIDへ変換する。元字面はsnapshotに残す。累積 $$ v_{j+1}=10v_j+d_j $$ の前に $$ v_j\le\lfloor(U-d_j)/10\rfloor $$ を検査し、 $$ U=2^{64}-1 $$ を越えればInvalidInputとする。字面とphysical IDの一意性検査は別である。

indexはnative整数nodeだけを受け、0以上int64上限以下を許す。DecimalNumber、String、Booleanを拒否し、`1.0`を1にしない。parserの整数eventの`-0`は値0として許す。indexを長さ/byte数へ使う際はchecked arithmeticとResourcePlanへ接続する。ID文字列をindexとして使う逆変換は提供しない。

## 有限の資源を、実際の計上単位で要求する

`InputLimits`は明示値を必須とし、暗黙の無制限や0=無制限を設けない。max_source_bytes、max_depth、max_nodes、max_members、max_decoded_string_bytes、max_decimal_token_bytes、max_total_text_bytesをuint64で持つ。0はその対象を許さない上限であり、必須rootを保持できない設定はResourceExceededになる。設定値のsize_t/int64境界を縮小前に検査する。

- source byte数はUTF-8の元bytes長。入力を複製する前に上限と照合する。
- depthはroot objectを1とし、containerへ入るごとに1増える。scalar rootはdepth0でもroot型不正である。parserがさらに深いcontainerへ進む前のeventで止める。
- nodesはroot、各memberの値、各array要素を各1と数える。member keyはnode数に重ねない。membersは全objectのmember数の合計とする。
- decoded_string_bytesは復号した一つのkeyまたはStringのUTF-8 byte数。total_text_bytesは所有する全復号key/StringとDecimalNumber tokenのbyte数の総和とし、元source長は別に数える。
- decimal_token_bytesはfloat eventで保持する一つのtokenのASCII byte数。native整数は固定64bit保持で、元字面はsource上限に含む。長い整数がfloat eventへ送られる経路にも同じtoken上限を適用する。

quotaの加算/乗算はcheckedで行い、超えた値を先に生成しない。node/member/string/token/source等の製品所有payloadをResourcePlanで確保前に計上する。vector容量・object key集合・parser stack等の同時に存在する領域を内訳へ含める。DOMを別に全保持して二重領域を隠さない。

`InputLimits`は7quotaに加え、`ExternalAllocationPolicy{mode, max_bytes}`を必須で持つ。modeはEstimateAdmissionまたはRequireProvenBoundの排他的な値であり、未知タグはInvalidInput、max_bytes=0も有効な零予算とする。これらは既存backendの外部確保と同じ受入区別で、新しいallocator方式の採用ではない。診断slotとResourcePlanを含め、呼出し引数の重なり/再入/同時共有を禁止する。

この入口自身が構文解析前に`plan_parser_memory(source_bytes, max_depth, quota)`を実行する。計画のconstructorは専用factoryだけに公開し、呼出側が自己申告の上界を構築しない。返す内部計画は、固定版とモデル版、同時領域の内訳、checkedなestimate_bytes、証明済み上界がある場合だけproven_upper_bytesを持つ。SourceSnapshotと製品treeはResourcePlan側へ、lexerのtoken_buffer/token_string、parserのstates、構文エラー生成時の一時文字列は外部計画側へ数え、同じ領域を両側で消したり二重保証したりしない。モデルの式と使用するコンパイラ/標準ライブラリ/nlohmann版・根拠sourceを実装時の契約補足と記録に固定する。任意倍率だけのestimateや、実測peakを普遍的上界へ読み替えた計画は採用しない。

EstimateAdmissionはestimate_bytesがmax_bytes以下の場合だけ進む。これは明示された推定受入であり、実確保やRSSのhard capではない。RequireProvenBoundは全外部同時領域のproven_upper_bytesがあるときだけ進む。証明がなければUnsupportedをparser呼出し前に返し、証明済み上界がmax_bytesを越えればResourceExceededとする。checked計算のoverflowもResourceExceededである。未知の必要量を0として通さない。両modeで実割当失敗はResourceExceededへ変換する。比較したplanned bytesは計画値として補助診断へ保持し、推定bytesをcore::Failureのobserved/実測量へ入れない。実際の必要/実測量が不明ならresources欄を省略する。

parser内部のtoken bufferとstackはSAX event前に確保し得る。従って7quotaのcallback検査を、外部確保の事前証明として使わない。成功文書は要求policy、外部計画/受入方式と明示保証範囲を所有し、推定受入を上界保証と表示しない。このkernelの製品payloadと外部parser領域の受入はここで閉じ、プロセス基礎メモリや後段の物理prepare領域まで保証したとはしない。

入力bytesを消費する同期処理中に、呼出側はbytesを書き換えない。処理に利用者callbackや可変Workspace共有を入れない。異なる入力・独立した資源計画は並列処理可能とし、同じResourcePlanの同時共有を許さない。

## 例外と診断を値で返す

返値は既存の`core::Result<InputDocument>`で、入口全体を`guard_entry(Operation::PrepareModel, InvalidInput, details, action)`で守る。`diagnostic`は呼出側が所有する`std::optional<InputDiagnostic>&`で必須引数とする。文書専用原因を確保不要の補助enumで返し、coreのFailure型・Result alias・guard_entryを変更しない。core::Failureの処理/識別/内側の既知詳細は保ち、Model/親/試行の不在をoptionalで表す。未知のIDを0で埋めない。

補助診断は、reason、source_id、取得済みならsource digest、判明した場合だけbyte_offset、event_ordinal、quota種別/上限を持つ。byte_offsetは元bytes先頭を0とする。parserが報告する1始まり位置は検査して換算する。EOFはnを指す。正確な位置を取得できない重複/キー/値変換では、架空のbyte位置を作らずevent_ordinalを使う。event_ordinalはparser callbackの文書内0始まり順で、未実行eventを数えない。

InputDiagnosticは確保なしで生成/複製/移動でき、parser message/what()/一時stringへの借用参照を持たない。呼出し先頭でdiagnosticをresetし、通常の失敗ではreason/判明位置を所有値として設定してからcore::Failureを返す。guard_entryが例外を変換した後、診断が空ならResourceExceededにはAllocation、その他にはUnexpectedExceptionの固定診断を確保なしで設定する。成功ではdiagnosticを空にする。位置やdigestが未取得ならそのoptionalを空にし、計算したことにしない。人向けmessageは固定reasonから別のformatterが生成する。JSON Pointerの追加確保を失敗構築に要求しない。

| 原因 | core codeと補助reason |
|---|---|
| JSON/UTF-8不正、BOM、重複key、root型/必須/未知key、型不一致 | InvalidInputと個別reason |
| schema文字列が未対応 | Unsupported / UnsupportedSchema |
| parserの有限数値範囲を超えたtoken | Nonfinite / NumberRange |
| quota超過、寸法表現範囲、bad_alloc | ResourceExceededとquota/Allocationの区別 |
| nlohmannの通常parse/type/rangeエラー | 上記の該当code。第三者enum/文字列を公開しない |
| その他の標準/非標準例外 | InvalidInput / UnexpectedException。空の成功へ置換しない |

最初にsource長と資源設定を検査し、snapshotとhashを取得してから構文/treeを作る。parse/重複/構造quotaの先行する失敗で止め、構文全体が通った後にschema型/対応版、未知root key、必須root keyの順に検査する。同種root不備が複数ある場合は固定key表の順を用いる。成功時も失敗時も原bytes、既存文書、Model/親状態を変更しない。失敗した処理自身の予約解放と累積payload計上はIM00-2の規約を維持する。

呼出側は返ったFailureと補助診断のコピーを一組として所有できる。同じ診断slotを次の呼出しで再利用すればそのslotだけが置換されるため、過去の失敗を保存するときは所有コピーを残す。結果もFailureもslotを借用しない。並列呼出しで同じ診断slotを共有しない。診断だけを見てcoreの成功/失敗を逆転させない。CLI exit/status保存はこのkernelの返値ではない。将来adapterは入力/未対応/Nonfiniteを2、資源を5へ対応させる。status保存に失敗した場合は6として別原因を保持する。未収束物理状態や架空のAcceptedStateを、この文書読取りから返さない。

## 独立literal corpusと初回実装gate

期待値は製品parserで再生成しない。短いUTF-8/JSON literalと明示treeを試験内に固定し、SHAは既存の独立hashlib/標準既知値へ照合する。数値についてはinteger/decimalのkind、値、元token/snapshotを厳密比較し、物理誤差の許容差を使わない。

| 対象 | 全成員を実行する対照 |
|---|---|
| 構文と文字 | 空、正常改行、日本語/emoji、正しいsurrogate対、BOM、UTF-16/32、不正UTF-8、単独surrogate、生NUL、escaped NUL、comment/trailing comma/token |
| 重複と順序 | 同じkey、escapeで同じkey、同じkeyを別objectに置く成功、親子同名、順序を入れ替えたobject、array順、unknown/missing root、全任意keyの有無 |
| 整数とdecimal | int64両端、uint64両端とその次、2^53の前後、-0/0、1.0/1e0、bool/string型不一致、1e400、1e-400、minsubnormal近傍、元token保持 |
| IDとindex | 0/最大、先頭zero、空/符号/空白/指数/小数/nonASCII digit、あふれ、index最大/負/float、IDとindexの混用拒否 |
| 資源 | 各quotaちょうど/一つ超過、depth/nodes/member/decodedUTF-8長/decimal長/総text、checked加算、0上限、推定受入の予算ちょうど/不足、証明の有無とproven予算ちょうど/不足、bad_allocと他例外注入 |
| 所有と診断 | 元buffer破棄/書換え後のdoc、結果移動/一時getter拒否、後続失敗後の既存doc不変、別source_id/同ID別内容hash、欠落位置optional、例外what寿命、独立並列入力 |

同じ木を別の空白/escape/改行で表した場合、意味値が同じでもsource digestは異なる。元snapshotからdigestを再計算し、再serializeしたJSONのdigestで代用していないことを検査する。複数不正を含むliteralで固定検査順も確認する。

意味ある隔離変異は、重複key検査省略、IDのdouble経由、元token破棄、元bytesでなく再serialize内容をhash、string長をUTF-8 byteでなく文字数で計上、sourceをborrowする変更から選ぶ。controlが通り、同じ入力をcompile/linkでき、実行assertが失敗したものだけを検出と数える。

初回実装開始条件は、本書の型境界・数値token・quota・例外/位置の契約、固定版parser外部計画の内訳/forecast式/保証可能範囲と、nlohmann 3.12.0のevent/例外位置/一時確保の実ソース対応を独立レビューで確認することである。実装は文書kernel、ID/index変換、対応する独立試験まで。parser資源の成立性に具体的欠陥があれば実装前に共有し、上限を黙って無制限へ変えない。

全nested schema validatorと参照snapshot解決、SI変換、CLI/Pythonのcase操作、runtime配布はIssue #22の後続契約で閉じる。採用済みのNastran/Gmsh・物理prepare・run/resume・結果/checkpoint・fresh環境の必須範囲を削らない。公開先/対応環境/署名はREL-1が持ち、本書で仮決定しない。ここまでに製品目的・精度・方式を再選定するowner判断はない。

## 読取り確認した正本とparserソース

本候補を作成したbaseは`8b6f2ae08cf5f55cbbd63b4fa740605ab4019adb`。既存の仕様・証拠は変更していない。以下は内容照合のためのSHA-256であり、凍結R&Dのsource-mapを書き換える指定ではない。

| ファイル | SHA-256 |
|---|---|
| `io-execution-adoption-specification.md` | `df926d9f6726a62af313e534e04c6a9b9aa4b838f7fa252567f3a50519039123` |
| `failure-certificate-contract.md` | `9555b700a595e8b80beec39147cc76ee25f27b9f31d0c070cb8fe5f184a99171` |
| `dimensions-resource-contract.md` | `dcd245ece27e7eec1bf446727813fc54f665f33d80222d679ed1208e5b419f89` |
| nlohmann `detail/input/parser.hpp` | `b4a762967b174de0c43896b465e4a7c8c92c3a8636f669debdec836050f7912a` |
| nlohmann `detail/input/lexer.hpp` | `a8b73e55242731375c088133f3a355af7b1bb82d6b9531678cb20560d3733c1d` |
| nlohmann `detail/input/json_sax.hpp` | `0ed6fdb0ef91f163cca312b7ffd6e0842f12e1b8340204ab2d0c05277e0d1f18` |

parserソースはロックしたvcpkg依存の`include/nlohmann/`以下で読んだ。`abi_macros.hpp`の3.12.0定義も照合した。`parser::sax_parse`のstrict末尾検査と`number_float`前の有限検査、`lexer::scan_number`の整数からfloatへのfallback、`lexer::get_string`のlocale小数点復元、`json_sax::number_float`のraw token引数を確認した。parserの`states`とlexerのtoken領域がSAXの外で確保されることも通読した。この候補作成ではparserの製品実行や資源上限の実測を済ませたとは主張しない。