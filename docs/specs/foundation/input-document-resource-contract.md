# 所有入力文書の型・parser資源の補足契約

2026-11-17。Issue #93、IM06-2。[入力文書契約](input-document-contract.md)を内部実装へ移す補足である。IO4の入力意味、全nested schema、物理資格を変更しない。型・資源式を独立審査し、先行する依存補修のMSVC対象検証成功を受けて、文書kernelの実装開始を承認済みである。Linuxの実行・標準lib識別が未確認の部分を成功とはしない。

## 所有型と文書内の参照

`SourceId`はIO内の専用tagを持つuint64識別子であり、0も有効とする。coreの物理ID・世代・Indexを流用しない。`SourceDigest`は32byteのSHA-256値であり、既存Boost.Hash2を使う。hash対象は保持した元bytesの全長だけで、path、末尾NUL、再serializeしたJSONを加えない。

`InputDocument`はmove-onlyであり、元bytes、source識別とdigest、全node、子edge、object member、復号textとdecimal tokenをそれぞれ固定長の`PlannedBuffer`で所有する。ResourcePlanと呼出し元bufferを破棄しても結果は生存する。移動元は文書不在であり、正常な空文書や有効caseとは扱わない。lvalueの文書だけから借用viewを取得できる。

nodeはNull、Boolean、String、SignedInteger、UnsignedInteger、DecimalNumber、Array、Objectの排他的variantとする。ordinalは0始まりint64だが、文書が生成する専用型にし、公開数値constructorを持たせない。利用者が任意ordinalや別文書のordinalを文書へ入力して未検査参照を作るAPIを設けない。公開node viewは生成元の不変storageとordinalを一組で借用する。数値ordinalの取得は報告用であって、他文書の同数値との同一性の証明ではない。

Array/Objectの直系子は、parentごとのhead/countとedge/memberのnextによる不変列で保持する。nested SAXのpreorder node区間を直系子の連続区間とみなさない。arrayの要素順、objectのmember順を入力どおりにたどれる。object memberは復号keyのtext参照と値nodeを持ち、keyとStringのescaped NULも長さ付きで保持する。DecimalNumberはparserのfloat eventに渡されたASCII tokenを所有し、近似doubleを公開値へ使わない。

`InputDocument`が持つ成功資格はroot envelopeだけである。`ValidatedCase`やModelへ暗黙変換せず、nested値のshape/SI/材料/参照先を検査済みと呼ばない。既存の厳格scalar入口へ渡す場合も、この文書のnative整数kindまたはStringを先に確認する。

## 製品payloadを先に計上する

source長をnとする。全7quotaとsize_t/int64への縮小可能性を検査し、SourceSnapshotを作る前に専用factoryのparser計画を受け入れる。snapshot/hashの後、構文検査を行う。7quotaは、source長、container深さ、node総数、member総数、一つの復号key/String長、一つのdecimal token長、所有text総長である。

初版はSAX中に容量を伸ばさない。保守容量は次のように設定する。

| 固定配列 | 容量の上限 |
|---|---|
| 元bytes | n |
| nodeとarray edge | min(max_nodes,n) |
| object member | m=min(max_members,n) |
| 重複key用の固定hash索引 | checkedな2m+1slot。親objectと復号key全bytesで一致を確定 |
| 復号key/String/decimalの共有text領域 | min(max_total_text_bytes,n) |
| containerの組立stack | min(max_depth,n) |

各値が少なくとも一つのsource byteを消費し、key/Stringの復号は元字面より長くならず、保持tokenも元の非重複部分であることから、nは安全な保守容量となる。hash索引の2m+1は空slotを少なくとも半数残すopen addressing用で、hash衝突を重複と扱わない。配列間で同じquotaを共有する場合も、それぞれの実sizeof×容量をchecked算術で課金する。input sizeだけを掛けた任意倍率でpayloadを隠さない。配列、control metadata、重複key索引の具体sizeofは実装と記録へ固定する。

callbackは書込み前に対応quotaと配列容量を照合する。key重複は同じ親objectの復号keyだけで照合し、次のmemberを追加する前に拒否する。hash等で索引を作る場合も、最終的には全UTF-8 bytesの一致で決める。親子や別objectの同名keyは重複ではない。

ResourcePlanは成功した確保を累積計上する。後続構文失敗や重複失敗で、この累積を過去の値へ戻さない。失敗した確保自身の予約解放は既存ResourcePlanに任せる。元入力と既存文書を書き換えない。

## parserが先に確保する領域

固定nlohmann 3.12.0のstrict SAXを使い、DOMは別に作らない。lexerの`token_string`、`token_buffer`とparserの`states`はcallback前にも確保される。decoded文字数quotaを、それらの事前容量証明へ読み替えない。

`lexer::get()/unget()/reset()`からtoken_stringの有効要素数はn以下、復号と数値scanからtoken_bufferはn以下となる。`states`はcontainer開始callbackが成功した後にpushするため、要素数は $$ d=\min(max\_depth,n) $$ 以下である。保守容量はparserが実際に消費する元source長から計画し、構文不正の場合も同じ外部予算を受け入れてからparserを構築する。

厳格な全消費には補助検査が必要である。固定lexerはraw NULをend-of-inputに分類するため、snapshotの全bytesからraw zeroを構文不正として拒否してからparserへ渡す。`\u0000`というescapeを拒否せず、その復号値を長さ付きStringに保つ。先頭BOMもparserの自動skipへ任せず拒否する。これは独自JSON parserではなく、既契約の不許可bytesを検査する入口である。

### OS別capacityと再確保の同居

Vchar(n)、Vbool(d)、Schar(l)は、固定版のcontainerごとに旧領域と新領域の同居まで加えたpayload forecastとする。allocator管理領域、全プロセスのRSS、stackの保証を表す量ではない。

| 実装 | sourceにある成長式 | forecastへ含める領域 |
|---|---|---|
| MSVC vector | max(requested,old+floor(old/2)) | char vectorとvector<bool>のword vectorの旧・新capacity |
| MSVC string | max(requested OR alloc_mask,old+floor(old/2)) | NULを含む旧・新char領域。Debug proxyは別slot |
| libstdc++14 vector | size+max(size,追加個数) | 旧・新capacity。boolはword単位へ切上げ |
| libstdc++14 string | requestedがoldの2倍未満なら2old | NULを含む旧・新char領域 |
| libstdc++14例外COW文字列 | 長さ+NUL+_Rep、必要に応じ4096byte単位への増分 | std::runtime_errorに保持するmessageの実要求領域 |

モデルはMSVC x64とlibstdc++14 x64へ限定する。MSVCのproxy量を、iterator debugが有効ならsizeof(container proxy)、無効なら0としてpと書く。各領域について、実size以下のold capacityと最大2倍までのnew capacity、stringのmaskとNULを含む次の保守式を使う。small-stringのinline部分も式上は外部領域として過大に計上し、必要量を小さくしない。

| 領域 | MSVC x64 | libstdc++14 x64 |
|---|---|---|
| Vchar(n) | 3n+p | 3n |
| Vbool(d) | 3ceil(d/32)×4+2p | 3ceil(d/64)×8 |
| Schar(l) | 3l+32+p | 3l+2 |
| R(l)：例外messageのcopy | l+1 | l+1+24+4096 |

Scharの係数3は旧capacityと新capacityの同居、MSVCの32は旧・新の各mask/NUL分である。内訳を次に固定する。

| 容量要素 | 個別の保守量とsource |
|---|---|
| MSVC旧char領域 | l+15+1。旧capacityの保守枠、alloc_mask15、NUL1 |
| MSVC新char領域 | 2l+15+1。1.5倍成長を2l以下で包み、mask15とNUL1を別に加算 |
| MSVC string proxy | iterator debug有効時に1個。xstringの_Alloc_proxy |
| MSVC char vector proxy | iterator debug有効時に1個。vectorの通常_Mypair側 |
| MSVC bool vector proxy | bit-containerの_Vb_val側1個と、内部_Myvecのvector<word>側1個。vector2833付近と通常vector672付近 |
| GNU旧・新char領域 | l+1と2l+1。二つのNULを含め3l+2 |
| GNU COW message | l+1+sizeof(_Rep)。固定x64の_Repは24byte、page増分は4096byte以下 |

SSOだけで終わる文字列には外部char領域が存在しないが、同じ保守式へ含め、過小に見積もらない。proxyはchar capacityやSSOとは別の確保であり、上表のpへ一度ずつ加える。再確保を実行しない小さい値へも同じ保守式を使う。未知構成を既知のOS名だけで同じ式へ通さない。model版と標準lib版・構成・根拠hashを計画に保持する。

### 構文エラー生成のlive slot

制御byteを`<U+XXXX>`へ展開する`get_token_string()`の返却長は最大8nである。失敗時には通常token領域が残ったまま、次の呼出しDAGの文字列が存在し得る。引数の評価順へ依存して同居領域を消さない。

| 経路 | 個別に数える文字列とcopy |
|---|---|
| lexer/parse_error | SAX引数用escaped token、exception_message内の別escaped token、last-read concat、error_msg、while-parsing concat、expected concat、context文字列 |
| exceptionの名称 | class名文字列、idの十進文字列、name concat |
| parse_error位置 | line/column各十進文字列、position concat |
| 例外完成 | 空diagnostics文字列、最終w concat、runtime_errorが所有するmessage copy |
| 数値overflow | SAX引数用escaped token、別escaped token、number-overflow concat、例外名称、空diagnostics、最終wとmessage copy |

各経路のliteral長は固定sourceに現れる文字列から得る。最大はlexer error77、token名21、context16、例外class12、id3桁、line/columnは各20桁である。 $$ k=8n $$ と置くと、last-read concatはk+92、exception_messageはk+170、parse_errorの最終wはk+275、overflow messageはk+26、range最終wはk+60以下となる。先行する依存修正候補が追加する原文字列とcopy結果の二つも数え、構文経路を次式で計画する。Sは上表のScharである。

$$
\begin{aligned}
E_{syntax}(n)={}&4S(k)+2S(k+92)+2S(k+170)+2S(31)+2S(32)+S(16)\\
&+2S(34)+S(12)+S(3)+2S(58)+2S(20)+S(0)\\
&+2S(k+275)+R(k+275),\\
E_{range}(n)={}&4S(k)+2S(k+26)+2S(34)+S(12)+S(3)+S(0)\\
&+2S(k+60)+R(k+60).
\end{aligned}
$$

4S(k)は二つのescaped token生成で各原文字列と返却copy、2S(34)+S(12)+S(3)は例外名のconcat原/copy・class引数・id、2S(58)+2S(20)は位置concat原/copy・line/column、S(0)は空diagnosticsを表す。その他の係数2は各concatまたはexception_messageの原文字列と返却copyである。並行しない局所領域も和に残した保守forecastであり、実測peakや任意倍率ではない。計算の全加算・乗算・切上げをcheckedにする。

parser全体のforecastは $$ P(n,d)=Vchar(n)+Schar(n)+Vbool(d)+\max(E_{syntax}(n),E_{range}(n)) $$ とする。各V/Sには再確保時の旧領域を含む。MSVC Debugではstring/vectorのcontainer proxyも別途加算し、vector<bool>の外側とword vectorの二つのproxyを一つにしない。

専用factoryだけがPを生成する。EstimateAdmissionは計画値Pをcapと比較し、P<=capなら受け入れる。不足とchecked overflowはResourceExceeded。証明済み上界は初版では不在であり、RequireProvenBoundはparser呼出し前にUnsupportedを返す。Pをobserved量へ入れず、成功文書と補助診断へ計画値として所有する。

factoryは標準libのfamilyとrevision、ABI、pointer/size_t/word幅、CHAR_BIT、iterator debug構成を固定表へ照合する。既知の条件と一致しない場合はUnsupportedで止め、Windows/LinuxというOS名だけから既知式を選ばない。計画metadataはmodel版、parserの上流版と補修版、標準libの識別、構成、根拠sourceの役割とSHA-256を所有する。実際のCIで新しいrevisionを観測した場合も、単に数値を広げて既知扱いにせず、該当sourceと式を照合してから固定表へ追加する。

## 例外を観測できる境界

通常は`core::guard_entry`でbad_allocをResourceExceededへ翻訳する。第三者message/whatは同期処理中の分類にだけ使い、Failureや補助診断へ借用しない。error.idと取得可能なbyte位置から分類し、reason/source/digest/eventを所有値として残す。raw error.idを返却するAPIは追加しない。

MSVCの`std_exception.cpp::__std_exception_copy`はmessage長+NULをmallocで複製し、失敗してもthrowせず_Whatをnullのままにする。`vcruntime_exception.h::what()`はそのとき固定`Unknown exception`を返す。nlohmannが生成したparse_error/out_of_rangeだけに入口を限定し、通常生成messageの固定`[json.exception.<class>.<id>] ` prefixが失われたCRT fallbackを、通常id/byte分類より先にResourceExceeded/Allocationへ変換する。無関係な標準例外をこのfallbackへ分類しない。

private分類関数へ同じfallbackを与える対照は、実CRT mallocを失敗させた実行とは区別して記録する。C++ newの故障対照だけでCRT malloc全量を監査したとは書かない。libstdc++の例外COW文字列生成はbad_allocを伝播するので既存guardで翻訳する。

**例外到達の開始gate：** 固定MSVC Debugのvector/stringのdefault/move constructorは、std::allocatorの場合にnoexceptでありながらiterator proxyを確保する。nlohmannのstates/token_string/token_bufferとnamed string returnがこの経路を持つ。独立21行probeはcontrol0、bad_allocで外catch42へ届かずterminate77となった。

先行するIssue #98の依存修正候補はparser/lexer/string_concatの3header9箇所だけを対象にする。空値は同値のthrow可能な明示count constructorで作り、lexer自身の誤ったnoexceptを外す。get_token_string/exception_message/concatのnamed returnは、同内容を直接構築するcopy prvalueにして暗黙noexcept moveを避ける。JSON grammar、event意味、標準lib、ABI、global Debug設定を変更しない。

候補の独立probeではNRVO無効も検査した。正常、単独surrogate、数値overflow、80重container、長tokenと末尾不正の5literalについて、原依存と候補の結果・event数・通常error101/406は一致した。候補の200個のC++new地点に対する個別故障は全て外catch42へ届き、terminate0だった。これは当該JSON SAX経路の限定的実証で、全nlohmann機能やCRT malloc注入の実証ではない。

先行依存修正source `b11535dfb94f71b0e947e08e8d369a9df1a097cc` の正式port導入は専用prefixへ行われた。対象4CTestのMSVC Debug/Releaseは各4/4、SAXの実new故障はDebug200・Release121地点でcatchへ到達した。これは依存補修側の検証で、文書kernel自身の終了判定ではない。検出macroの抑制や故障注入対象からの黙った除外は行っていない。

## 固定順序と独立な対照

失敗順はsource長・quotaの表現範囲・policy/forecastの受入→snapshot/hash→raw NUL/BOMとstrict SAXの構文・重複・構造quota→存在するschemaの型/対応版→未知root key→必須root keyである。schema文字列不対応だけをUnsupportedにする。未知keyが複数ある場合は復号UTF-8 byteの辞書順、必須key欠落は入力文書契約の固定表順で選ぶ。schema自体の欠落を型不一致へ読み替えない。

診断slotを入口でresetし、source_id、取得済みdigest、取得可能なbyte位置またはevent ordinal、quotaと計画値を確保なしで保持する。parserの1始まりbyteを検査して0始まりへ変換し、EOFはnとする。正確なbyteがないcallback拒否はeventだけとし、未実行eventを数えない。内側FailureDetailsを保存し、不在の親や世代を0で生成しない。

全8kind、空/日本語/emoji/正負surrogate/BOM/UTF16/32/raw NUL/escaped NUL、全重複と順序、native整数両端/2^53前後/decimal fallback/overflow/underflow、全root key、全7quotaの等号/1超過/0、所有move/一時borrow拒否、同ID別bytesのhash、既存文書不変、独立並行呼出しを固定literalと明示treeで照合する。SHAは手書きbytesと独立hashlib既知値を照合し、製品parserから期待値を作らない。

source/hash、元token、重複、UTF-8 byte計上、所有の意味を壊す変異は、隔離したbyte同一controlと専用producerを確認してから行う。compile/link成功後のassert失敗だけを検出に数える。未変更scalar/numericsの試験を再実行せず、最終固定sourceで影響するD/Rを各1回へまとめる。実行OS、標準lib、parser sourcehashと未確認範囲を記録する。

## 根拠sourceと検証範囲

MSVCは19.51.36260.0、実header/CRT sourceは14.51.36231。GCCは一次sourceのreleases/gcc-14.2.0を取得して読み、Linuxの実行や宿主標準libとの完全一致を済ませたとはしない。実CIのcompiler/標準lib識別は別に記録する。source取得は環境や共有依存prefixの導入・変更ではない。

| source | SHA-256 |
|---|---|
| 原nlohmann parser.hpp | `b4a762967b174de0c43896b465e4a7c8c92c3a8636f669debdec836050f7912a` |
| 原nlohmann lexer.hpp | `a8b73e55242731375c088133f3a355af7b1bb82d6b9531678cb20560d3733c1d` |
| 原nlohmann exceptions.hpp | `5c9b21e0a13ac0bab8dea676b9804d739dfd2adc24f036f97ab243ef70e2aae7` |
| 原nlohmann string_concat.hpp | `fb30fb50c325a1b5610a62a9061a1966b4e850dd964635c1ad1c4c9da48f7daa` |
| MSVC vector | `b254d6fc250164a9f4e71ea66e2c1e9401cb2bf8b0056550fd92135518450b51` |
| MSVC xstring | `ab80a6b60196246d1c8653e6d55f78e3cbaaf4207860a37cf6996cdc833eea7a` |
| MSVC std_exception.cpp | `cc401650d7cd6fdefbfaec65b6cb83e504b2617a89bae3a42871ac9960eb09c3` |
| MSVC vcruntime_exception.h | `0dd3dcc67b655fbfc96708f131c9e7f86416d2db95fd33b1e7acc8489f2c32d4` |
| GCC stl_vector.h | `b896547442cc7fa6e6eea35f94a64e0b32431e1e3d8abee301bd5aa1f5ff8c0f` |
| GCC vector.tcc | `dc5dad860d9f9add7ccbbea1e67f315c6f0f74d6398ccdd6459391a976cd0dc2` |
| GCC stl_bvector.h | `a5e4a8743d833591f6809d5d2bdb9998f7f681c1248a5e2c49623676d29ab4ea` |
| GCC basic_string.tcc | `4f8dbb49a708622d994acc0b6ad3e69cf72a746732c4ff2362d078b4699d478d` |
| GCC cow-stdexcept.cc | `c40ec114b10ec3a67067d2e5f206cb31d4441afb713bb4d8789d7ec85e107f65` |
| GCC cow_string.h | `0cb671deeb0d99a38b6171290d056034ef6c2c804078a4fca30b4541f16f4d5b` |

私有候補と正式port導入後の3headerはhash一致を確認した。parser SHAは`8f219283e1b1a5ad3dc413acfa33c4ee3901c0366f60676c3a93b0d37aeda7b6`、lexerは`9ab10d7a6060e048842465836c9f11ec3b064125197f388e2e9792f7ca1211f2`、concatは`666b4e48add7381f4b995cd451d3923c1a6b261d633dd9c9ad6f80e54298020b`。将来のpatchでcomment等が変わる場合は新hashと対応を記録し、旧候補へ実行証拠を付替えない。

## GNU実環境の受入表への追補

2026-11-17、依存補修PR #101の[Linux Debug job](https://github.com/YUKIKEDA/UVLAS/actions/runs/44008552697/job/133789790410)が、実環境の識別JSONを出力した。採取対象はcommit `77f690414088ae2006442a4e289bb776d22cd27d`、packageは`g++-14`と`libstdc++-14-dev`の`14.2.0-4ubuntu2~24.04.1`である。`_GLIBCXX_RELEASE=14`、`__GLIBCXX__=20240908`、C++11 ABI=1、long/pointer/size_t各8byte、CHAR_BIT=8、C++23、debug系macro不在を確認した。この採取は環境識別であり、本入力文書のLinux試験成功を示さない。

実headerの`stl_vector.h`、`vector.tcc`、`basic_string.tcc`、`cow_string.h`は上表のGCC一次sourceと同じSHAだった。`stl_bvector.h`だけは`7abd919ef3579f6199fb87981f80239e25ec57bc0f80e16e0ca792d44bb14f29`と異なったため、[Ubuntu公式の同版deb](https://archive.ubuntu.com/ubuntu/pool/universe/g/gcc-14/libstdc++-14-dev_14.2.0-4ubuntu2~24.04.1_amd64.deb)を私有directoryへ取得した。依存の導入は行わず、該当headerだけを展開して実CIのSHAとの一致を確認した。deb全体のSHA-256は`6d223262b8cbe6d5150435056d58dffa8159d8582fc0df47fb45ea550cd37424`である。

差分は二つのconstructorだけだった。`_Bvector_impl`のdefault constructorは、requiresの条件へ`__glibcxx_type_trait_variable_templates`を加える。allocatorを受け取る`_Bvector_base` constructorは、`_M_impl(__a)`を`_M_impl(_Bit_alloc_type(__a))`へ変更する。SAXの`std::vector<bool>`は標準allocatorを使い、後者は同じword allocatorへの変換を明示する。word数への切上げ、保持layout、allocate/deallocate、`_M_check_len`の成長式に差分はない。このためVboolと他のforecast式・model revision 1を維持し、GNUの受入revisionへ20240908を追加する。返却するsource証拠のBitVector欄は、20240801には一次sourceのSHA、20240908には実Ubuntu版のSHAを保持する。未知revision/debug設定を許可する変更ではない。

Windows D/Rは追補前後に同じ製品source pathをMSVCで前処理した。追加したpreprocessor directiveの位置に四つの空行が増えたため、raw出力bytesは一致しなかった。空行だけを除去した後は、残る全行がbyte一致した。文字列内の空白やtokenを除去・変更した比較ではない。この結果からWindowsの前処理tokenが変わらないことを確認し、既存dbe8a03のWindows試験・変異証拠を保持して、数値試験を反復しない。rawと比較用表現の両hash・比較手順、および実CIの環境JSONは今回の検証記録へ別artifactとして残す。
