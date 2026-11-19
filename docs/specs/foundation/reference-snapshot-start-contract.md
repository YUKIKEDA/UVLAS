# 外部参照と所有snapshotの実装開始契約

ケースschemaの [参照解決](case-schema-contract.md) を、既存の [元入力文書](input-document-contract.md) と [SIケース](si-precision-conversion-contract.md) から独立して検査できる処理にする。Nastran・Gmsh・履歴・観測の消費処理が、同じ出典、元bytes、役割と依存関係を使う。

## 入力、所有権、受渡し

| 処理 | 論理I/Fと保証 |
|---|---|
| 参照を採る | `SchemaCase` の採用辞書のフィールド順と各配列の添字順で、元の `SourceId`、復号済みpath、宣言digest、role、元文書上のnodeと位置を取得する。JSONのキー出現順、実ファイルの列挙順、SI数量の並びを使わない |
| 参照を収集する | `ReferenceInput`、上限、資源計画、取消、任意の診断slotを受け取る。`ReferenceInput` は検査済み `SchemaCase` と、取得経路を表す `File` または `Memory` を組にして所有する。adapterは保持した根のハンドルに対して開き、利用者callbackを数値反復に残さない |
| 成功を公開する | `ReferencedCase` は親 `SchemaIdentity`、宣言元の位置、元と正規化後のpath、role、元bytesとdigest、依存辺を所有する。入力、読み手、Planを破棄しても内容は保持する。公開の借用viewは親結果の寿命まで有効 |
| 準備で消費する | `SiCase` と `ReferencedCase` の `SchemaIdentity` を比較してから、`(SourceId, role)` でsnapshotを得る。消費者は収集時に読んだ配列だけを読む。検査後にpathを開き直さない |

論理入口は `resolve_references(const ReferenceInput&, SnapshotLimits, ResourcePlan&, std::stop_token, SnapshotDiagnostic*) -> Result<ReferencedCase>` とする。外部参照の収集とSI変換は同じ `ReferenceInput` の `SchemaCase` から始める。SI変換後に参照を収集する場合も元の `ReferenceInput` を渡す。`SiCase::source_case()` 単独からファイル由来の情報を復元できない。CLIの入力検証は、参照整合のためにSI変換用の許容誤差を要求しない。

初期参照の採取対象は次の8役割とする。列挙器の生成・検証は辞書と同じ入力を使い、現行辞書の全 `FileRef` の型と実際の到達経路を独立に照合する。

| 起点 | role |
|---|---|
| Nastran Geometryのdeck、sidecar | `nastran-deck`、`nastran-sidecar` |
| Gmsh Geometryのmesh、sidecar | `gmsh-msh4.1-ascii`、`gmsh-sidecar` |
| ReferenceGeometryのsupplied | `reference-geometry` |
| BranchSeedのprepared | `preparation-history` |
| Initialのprepared_history | `initial-history` |
| ComparisonReferenceのfile | `observation-reference` |

ファイル内の依存宣言は、そのroleの採用形式が認める場合に限り、検査済みsnapshotと宣言位置から抽出する。抽出器はファイルを開かない。役割固有のパス規則がある場合は、抽出器が規則を適用し、収集処理へmanifest根に相対化した宣言を渡す。採用形式が未定のroleは、依存抽出を認定するまで実製品の全参照成功を返せない。単に抽出器が無いことを依存ゼロの証明にしない。

## パスと実体

`File` は取得処理だけが生成する。外部参照の根にはmanifestを含む直上のディレクトリを使う。取得処理はその `ReadRoot` からmanifestをファイル名一成分で開き、同じハンドルから元bytesと実体ID（Linuxはdevice/inode、Windowsはvolume/file ID）を取得する。読取りの前後でID、長さ、変更指標を照合する。その元bytesから検査済み `SchemaCase` を作り、同じ `ReadRoot`、実体ID、manifestハンドルと組にして解決器へ渡す。呼出し側が別の `SchemaCase`、根、IDを組み合わせる入口は設けない。`File` の寿命中は元ハンドルを保持し、解決中にファイルIDが別の実体に再利用されるのを防ぐ。パスを再度開いてIDを取得しない。安定したIDを取得できない環境は `Unsupported`、manifestの読取り失敗または取得中改変は `InvalidInput` とする。

`Memory` は呼出し側のbytesから作る。外部参照がある場合は呼出し側が基準となる `ReadRoot` を指定する。参照がゼロなら根を省略できる。`Memory` はmanifestの実体IDを持たない。`SourceIdentity` の呼出側ID、長さ、digestは両経路とも元bytesの識別に使う。外部ファイルとの実体照合には、開いたハンドルから得たIDだけを使う。

正規化ではUTF-8の元表記を保持したまま、相対パスの `/` と `\` を区切りとして扱う。`.` を除き、`..` では正規化済み成分を一つ戻す。空スタックから戻る場合は脱出、結果が空または末尾が区切りの場合はファイル指定の欠落とする。NUL、絶対パス、ドライブ相対、UNC、Windowsのストリーム指定を拒否する。余分な空成分は除く。大小文字、Unicode正規形、空白は意味を変える変換を行わない。Windowsで同一視される名前を黙って別名扱いせず、開いた実体のIDで確認する。

正規化後の解決では根のハンドルを使い、途中の成分をディレクトリ、最後の成分を通常ファイルとして扱う。途中も含めてリンク・reparse pointを拒否する。現在の作業ディレクトリ、文字列前方一致、事前のcanonical path文字列だけでは根への所属を認定しない。

- Linuxの実装候補は根のdirfdを使う `openat2` の `RESOLVE_BENEATH | RESOLVE_NO_SYMLINKS`。通常ファイルであることをハンドルの属性で確認してから内容を読む。検査中に親が移動できる状況で、成分ごとの `openat` だけを同等と見なさない。[Linuxの解決規約](https://man7.org/linux/man-pages/man2/openat2.2.html)
- Windowsの実装候補は `NtCreateFile` の `RootDirectory` に検査済みの親を渡し、一成分ずつ `FILE_OPEN_REPARSE_POINT` で開き、得たハンドルの属性でreparse pointを拒否する。根と途中ディレクトリは `FILE_LIST_DIRECTORY` を含む権限と、rename/deleteを許可しない共有条件で開いて保持する。属性確認だけのハンドルではrenameが通ることをNTFS上で確認した。検査例と再現手順は [Windowsの最小対照](../../verification/IM-06/reference-snapshot-windows-probe.md) を参照する。[相対openとファイルオプション](https://learn.microsoft.com/windows/win32/api/winternl/nf-winternl-ntcreatefile) [アクセス権と共有条件](https://learn.microsoft.com/windows/win32/api/fileapi/nf-fileapi-createfilew)

上記と同等の検査を実機で成立させられない環境では値として失敗する。配布対象の最低OS版はREL-1が決める。パス解決後は同じ最終ハンドルを読取りに使い、成功・失敗のどちらでも一時ハンドルを解放する。

開いた同一ハンドルの通常ファイルについて、読取り前の実体ID、長さ、変更指標を取得する。長さの上限とPlan残量を確かめてから所有先へ読む。一回の読取りが要求長より短くても、正の進捗があれば続ける。実長が不足した時点でEOFになった場合、または宣言した容量を読み終えた時点で次の1 byteが得られた場合は拒否する。読取り後の実体ID・長さ・変更指標にも差があれば取得中改変とする。固定した元bytesにSHA-256をかけ、宣言digestと一致する場合だけ次段へ渡す。digest照合後、子依存を開く前に今回の途中成分と最終ファイルのハンドルを解放する。DFSに残すのは所有した内容と出典である。この検査が保証するのは消費するbytesの一致と観測した変更の拒否である。全ファイルシステムで観測不能な書込みまで検出できるという保証は置かない。

## 重複と依存

同じSourceIdの再宣言は、元のpath・digest・roleが同じ場合に一つの保存先を使う。宣言の出典と依存辺は保持する。各参照先を開いて実体IDを得た直後、`File` ではそのIDをmanifest自身のIDと比較する。一致すればパス名や宣言digestにかかわらず、その宣言位置を指す `Cycle` を返す。`Memory` ではこの比較を行わない。次に、異なるSourceIdが同一の通常ファイルIDを指した場合は重複として拒否する。内容が同じ別ファイルは、両方の出典・宣言が正しければ許可する。

宣言の辺は **消費するsource → 依存するsource** とする。初出順に、未訪問・処理中・完了の3状態で反復DFSを行う。処理中点への辺はcycleとして拒否し、完了点への参照を共有する。完了順は依存先が先になる。各頂点と宣言辺を一度処理し、DFSは明示stackを使う。診断には最初の不一致または循環を指す元source、role、宣言位置、必要ならcycle列を保存する。

## 資源と失敗

上限は元宣言数、固有source数、辺数、正規化後の全path bytes、path成分深さ、グラフ深さ、単一ファイルbytes、全ファイルbytes、走査・読取り・hash work、診断bytesを別に指定する。0を無制限として扱わない。byte容量と要素数はchecked加算・乗算で求め、`Index` とPlan残量を検査してから確保する。既知のroot宣言はpathとIDを先に検査し、全部を環境から読み出す前に解決できる不備で止める。

固有source $$ N $$ 件、宣言辺 $$ E $$ 本、pathと出典のコピー $$ T $$ bytes、読み取る固有ファイルbytes合計 $$ L $$ に対し、所有payloadの計算は次の項を含む。

$$
P=\operatorname{bytes}(N\text{ source records})+
\operatorname{bytes}(E\text{ declaration records})+
\operatorname{bytes}(N+2\text{ adjacency offsets})+T+L.
$$

$$ P $$ は完成物の容量である。オフセットはmanifestを表す仮想頂点、固有sourceの各頂点、末尾の番兵を含む。成長時の旧新buffer、lookup、明示stack、hash状態、OSハンドル、診断を別に見積もる。子宣言は検査済みbytesを読んだ後に初めて分かるため、段階入場とする。辺・頂点・pathの追加で生じる確保はその追加前、ファイルのbytesは同一ハンドルから長さを得た時点で容量を検査する。一時的に旧新領域が同時に必要なら両方が入ることを確かめる。明示DFSの生存frame数はグラフ深さで制限する。解決中のWindowsの生存ハンドル数は根1・`File` のmanifest 1・途中成分数・最終1で制限する。読取りは所有するbytesへ直接行い、別の全ファイル複製を確保しない。宣言処理／比較、path文字走査、読取り・hashそれぞれのbyte数と最後の余剰読取り検査をworkに計上する。allocator・OS管理領域は別の費用として扱う。

最初に親と制限、root宣言を検査する。その後は宣言順でパス、実体（manifest自身と重複を含む）、容量、読取りと変更、digest、子宣言とDFSの循環の順に失敗を確定する。取消は開始時・I/Oの前後・公開直前で検査する。失敗は同じ原因を診断slotなしでも返す。診断の容量不足は診断の詳細だけを省き、元の拒否理由を置き換えない。成功済みのPlan課金と、解放する未公開payloadは区別する。

失敗のoperationは `PrepareModel` とする。path違反、読取り不能、重複、循環、取得中改変、digest不一致は `InvalidInput` とし、詳細の理由を分ける。上限と確保失敗は `ResourceExceeded`、根相対の安全な解決またはroleの検証方式を提供できない場合は `Unsupported`、取消は `Cancelled`、準備時に別の親世代へ接続した場合は `StaleGeneration` を返す。診断は宣言元の位置、拒否した成分、判明した実体・長さ、quotaと必要量を所有し、実測できない値は省く。

## 開始判定と検証

グラフ、宣言表、所有bytesの組立は、上記の論理I/Fから実装に着手できる。一般roleの全ファイル解決成功は、Windows/Linux adapterと対応roleの抽出器が次の検証を通過した後に公開できる。Importerの節点・材料・物理等価性の認定はIM06-0の各後続段階で行う。

| 検証 | 独立な期待値と観測 |
|---|---|
| 入力の結線 | 8役割を全部持つ正常例と参照ゼロの例。辞書から手で転記したSourceId・role・位置と厳密一致 |
| 所有と照合 | [元bytesとグラフの例](../../verification/IM-06/reference-snapshot-design-cases.json)。digestは独立SHA-256、元ファイルと根を破棄した後も全payload一致。異なるSchemaIdentity・roleの借用は拒否 |
| manifestの実体 | [固定例と照合手順](../../verification/IM-06/reference-snapshot-manifest-origin.md) のとおり、`File` の元名とハードリンクは `Cycle`、同じbytesの別実体は成功。`Memory` に同じbytesと根を渡した場合、ハードリンクは通常の外部参照として成功。呼出側ID・長さ・digestは共通とする。別の例で、上位と直上のディレクトリに同名ファイルを置き、直上を根に選ぶことを確認する |
| 同名・別名・循環 | 明示グラフ、hardlink、二つのpathの同一実体、同じ内容の別実体、root manifestへの辺を両OSで照合。親から子の順、完了順、循環の宣言位置は手書きの列と厳密一致 |
| パスと取得 | root内での`.`・`..`、root外、絶対・ドライブ相対・UNC・stream、symlink/junction、途中の親の改名、短いread、EOF前の改変、長さ不変の別bytes。制御したadapterまたは別プロセスで変更時点を固定し、時間の偶然に頼らず拒否理由を照合 |
| 資源と故障 | 各上限ちょうどと1不足、整数overflow、確保点ごとの故障、診断予算ゼロ、最後の確保中の取消。旧snapshotの全payloadと読み手のハンドル数を確認 |

この開始契約では数値精度・実務受入を認定しない。上記の実機対照、role parser、外部形状・履歴の物理的妥当性は製品実装の検証に残す。
