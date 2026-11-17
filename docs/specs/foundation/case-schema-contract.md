# 全ケースschemaと段階別検証の契約

2026-11-17。Issue [#110](https://github.com/YUKIKEDA/UVLAS/issues/110)。IO4-01の `uvlas.case/1` を、型検証の実装へ渡す追加契約である。構造・空力・J1の採用式を変更しない。旧仕様のバイト列は保持する。本書と[フィールド辞書](case-schema-fields.md)、[設計検査](case-schema-design-cases.md)を一組として読む。

## CS-01 対象と検証段階

`InputDocument` は構文、元のバイト列、値のkind、root envelopeを所有する。全schemaはその後の処理であり、物理kernelの完成を待たず定義する。根拠は[作業分解のIM06](../../plans/next-package-work-breakdown.md#im06入出力実行再開)と[入力文書契約](input-document-contract.md)である。

| 段階 | 入力→所有する成功型 | 成功が保証すること | ここでは保証しないこと |
|---|---|---|---|
| 構文 | 元bytes→`InputDocument` | 既存文書契約 | nested schema、物理値 |
| 全schema | 文書→`SchemaCase` | 本書の全key、kind、shape、variant、IDの種別と一意性、内部参照整合、未SI値の意味 | 外部参照ファイルの存在、SI丸め、物理適格性 |
| SI変換 | `SchemaCase`→`SiCase` | 次元、単位倍率、基底、値の変換誤差、解決した既定値 | 外部bytesの真正性、幾何/rank |
| 参照解決 | `SchemaCase`と`SiCase`の参照表→`ReferencedCase` | 同じsnapshotのbytes/hash/出自、役割・依存関係 | importerの物理等価性 |
| 準備 | 同じcaseのSI値・参照snapshot→`PreparedModel` | 材料、幾何、拘束、基底、rank、容量等の各採用契約 | 数値解の達成、実務受入 |
| 実行 | 準備版・初期/再開状態→既存の試行/確定結果 | 同じ離散問題の数値認定 | 未評価用途の妥当性 |

これらは別型である。暗黙変換、共通の `valid=true`、未検査nodeを含む部分成功を設けない。CLI `validate` の成功には構文、全schema、参照整合が必要であり、schema単体を接続しない。SI変換・snapshotの実装はIM06の後続依存段階、物理準備はIM01〜05との統合へ渡す。

## CS-02 字句型と閉じたobject

フィールド辞書の各行はobjectの**全許可key**を列挙する。`?` の付いたkeyだけ省略可能であり、その他は全て必須とする。記載されていないkeyは、階層やvariantにかかわらず拒否する。全型で `null` を拒否する。欠落、空配列、零、`false`、空文字列を代用しない。`none` が意味を持つ箇所はタグ付きvariantとして明示する。

| 記法 | JSONと意味 |
|---|---|
| `Text` | UTF-8 string。表示文には空文字を許す。識別名・列名には空文字を拒否。NULは表示文以外で拒否 |
| `Tag{a,b}` | 大文字小文字を区別する列挙string。未知値を拒否 |
| `Bool` | JSON booleanのみ |
| `Id<T>` | 10進stringの全uint64、0も有効。既存strict scalarを使う。Tは名前空間を識別し、別種IDで代用しない |
| `signed-native-integer` | JSON native signed/unsigned integerでint64へ正確に変換可能。負値可。bool/decimalを拒否 |
| `Index` / `Count` | native signed/unsigned integerのみ。前者は0〜INT64_MAX、後者は0〜UINT64_MAX。`1.0`、bool、stringを拒否 |
| `PositiveCount` | `Count`かつ1以上。配列長と掛算はchecked arithmeticで照合 |
| `Number` | native integerまたは`DecimalNumber`の元token。boolean/stringを拒否。binary64へ変換せず保持 |
| `Q<d>` | 次元dを持つ`Number`。root unitsに従う。dは辞書の型で決まり、利用者の自由な単位stringを各値へ重ねない |
| `Rational` | `{numerator: Text, denominator: Text}`。numeratorは`-?[0-9]+`、denominatorは`[0-9]+`かつ厳密正。既存PG exact scalarへ渡し、分母0/桁数上限超過を拒否 |
| `Vec<T,n>` / `Mat<T,r,c>` | 長さnのarray / r個のrow array（各c個）。常に行優先の論理順。平坦化や転置を推測しない |
| `List<T>` | Tのarray。空可否は使用箇所に記す |
| `Ref<T>` | 対象名前空間TのID。参照先は同じcase内に存在すること。配列indexをIDにしない |
| `Digest` | 小文字16進64文字のSHA-256。元bytesに対するもの |
| `FileRef<role>` | `{source_id: Id<Source>, path: Text, sha256: Digest, role: Tag{指定role}}`。役割を変えた使い回しを拒否 |

`Number`は元token、sign、digits、十進指数、出典nodeを所有する。十進数の値は既存文書契約の式を使う。例えば `1e-400` と `1.000000000000000000000000000001` を0や1へ先に丸めない。schemaの大小・等値・整数条件は元十進値の厳密比較で検査し、桁数/指数の資源限界超過を資源失敗として返す。指数の巨大さを扱うために巨大な10の冪を無条件に生成しない。

有理係数とJSON Numberは別型である。PGの厳密なrank根拠をdecimal→binary64→有理復元で作らない。元decimalを厳密有理化する経路はSI変換契約で定め、入力が明示した有理係数の原形も保持する。

## CS-03 所有・同一性・既定値

`SchemaCase`はmove-onlyで、元`InputDocument`と型付きpayload、参照表、単位系、未SIの値、出典位置を所有する。元文書の破棄後も独立に存続する。内部借用は所有者のconst lvalueからのみ取得し、rvalueからview/referenceを返さない。公開不変型へ未検査nodeのescape hatchを作らない。

型付き結果は元bytesの`SourceIdentity`とschema版を持つ。同値値を別字面で記述した文書は同じsource digestにはならない。schema結果の同一性は元digest、schema契約版、型付き内部表現版で識別し、物理Model fingerprintと同一視しない。内部参照に使うIDは利用者値であり、hashの代わりではない。

省略可能keyの全てに辞書で解決規則を付ける。物理値、材料、減衰、密度、支持剛性、基準量、精度、資源に暗黙の既定値を置かない。省略した `provenance` は空の来歴補足、`initial` と `restart` はanalysis条件で処理する。`resolved-case.json`は後続段階が元値、SI値、変換誤差、適用した構造上の既定値を区別して生成する。schema処理は元文書を書き換えない。

## CS-04 参照・IDと検査順

IDは `Geometry, Material, Shell, Beam, Joint, Support, Mass, AeroSurface, Layer, Load, Observation, Source, Frame, Branch, WakeRow, WakePoint, WakeBand` の各名前空間で一意である。Sourceは各FileRefの初出で文書全体のSource registryへ登録する（rootに別collectionは置かない）。複数FileRefから同じ定義を再参照できるが、同じSourceIdに異なるpath/hash/roleを宣言してはならない。同じ数字を別名前空間で使える。型付きtargetがshell/beamのどちらかを表す場合は、`kind`と対応するIDを両方持つ。名前空間をまたぐ曖昧検索をしない。出典が同じファイルであることと、物理entityが同じであることを混同しない。

検査の順序は次の一か所で定める。

1. 契約版・呼出し上限の正当性と文書identityを検査する。
2. rootから文書のmember順・array順で全key/kind/shape/variantを走査し、quotaを検査する。
3. 各値を未SI型へ所有コピーし、exactな字句条件、範囲、順序、排他を検査する。
4. 全IDを登録して重複を検査し、同じ順で型付き参照を解決する。循環禁止の関係は全グラフを検査する。
5. analysis/initial/restart、model/flow/observations/numerics間の条件付き必須と次元契約を検査する。
6. 全要素を検査してから一度だけ`SchemaCase`を返す。最初の失敗を返し、途中結果は公開しない。

重複weldや物理拘束の従属はschemaの単なる参照循環と異なる。参照座標/材料継承/参照ファイル依存は循環禁止、物理接続グラフは閉路が存在してもschemaだけで削除しない。固定・非線形拘束の大域整合はCN-D04のprepareで判定する。

パスはNUL/空文字を拒否し、`FileRef`のpathはmanifest directory基準である。schemaでは文字列と役割を検査する。存在、正規化、親ディレクトリへの脱出、symlink/junction、同じfileの重複・循環、取得中改変は参照snapshot契約で検査し、schema成功をその代替にしない。出力先pathは読み込む`FileRef`と別型である。

## CS-05 失敗・資源と例外

入口の論理I/Fは `validate_case_schema(const InputDocument&, SchemaLimits, ResourcePlan&, SchemaDiagnostic*) -> Result<SchemaCase>` とする。外部状態・ファイルを変更せず、GIL/callback/物理評価を必要としない。`SchemaDiagnostic`は入力source/digest、node index、key path、期待型、実kind、失敗理由、参照ID/対象種別、quotaと実必要量を所有する。存在しないbyte位置を創作しない。key pathは元のdecoded keyとarray indexから作り、JSON Pointerの`~0`/`~1`を用いる。診断slotは失敗後も借用に依存しない。

| 原因 | core分類 | 補助理由 |
|---|---|---|
| 未知/欠落key、kind、shape、variant、範囲、重複/参照先欠落 | `InvalidInput` | fieldと期待条件 |
| 未知schema版、未採用tag/表現 | `Unsupported` | 対応版/能力 |
| quota/checked寸法/予約容量の不足 | 既存の寸法/資源分類 | quota、必要量、上限 |
| 外部から非有限値を注入した不正document | `Nonfinite` | 入力位置。通常の文書入口では発生しない |
| allocation/想定外例外 | 既存guardの失敗値 | 診断を生成するため再度巨大確保しない |

`SchemaLimits`は元文書上限とは別に、typed payload bytes、全entity件数、全参照件数、全scalar件数、全多項式term件数、全table sample件数、式node/edge数・深さ・整数冪絶対値・評価点数、全text bytes、診断path bytesの最大値を全て明示する。無制限の既定値を設けない。元documentのquotaを後段の追加確保上限と呼ばない。

固定型metadataを除く成功payload、辞書、ID索引、参照配列、exact比較用の可変scratchは`ResourcePlan`で確保前に計上する。入力文書を保持する所有コピーも計上し、借用だから無費用とはしない。最大同時scratchと成功payloadのforecastを検査してから構築する。cumulative計上は既存契約どおり破棄時に返金しない。std容器等の外部確保を使うなら固定版の計画/推定受入と真の上界要求を分け、根拠なしにhard capと呼ばない。

quota失敗でも元文書、呼出し前のModel、Accepted/Trialは変更しない。`ResourcePlan`の過去の成功確保計上は巻き戻さない。並列呼出しは各workspace/diagnostic/planの所有規約に従う。

## CS-06 採用仕様への対応と開始条件

| 採用元 | 本契約で型へ落とす量 | 後続で数値/物理検査する量 |
|---|---|---|
| IO4-01 | root全9必須/3任意、model全8集合、全nested辞書 | SI、snapshot、prepare/実行 |
| SH-D01〜07 | PG係数、領域、厚さ、層、分割、アンカー、共有/折れ | 有理rank、幾何/trace/求積誤差 |
| MAT-D04〜12 | 全3D則、密度、軸、E0、層、出自 | SPD/自然計量/実材料適用性 |
| BM-D01〜09 | SE(3)参照、K6、質量/偏心/慣性、初期ひずみ、反り | Log分岐、断面/慣性適格性 |
| CN-D01〜09、SI-D01〜09 | 接着/取付/支持、質量台帳、荷重/減衰/回復 | 全trace/依存拘束/仕事/全接線 |
| UI/UW、AE-07 | 面/後縁/共有、核、流れ/補間、移流/窓 | 幾何/自己/誘導/荷重誤差 |
| J1 | static/transient/flutter、枝、区間、初期/再開、探索 | 同じ残差、BDF2履歴/微分、coverage |
| N4、IO4-05/06 | 精度/誤差/予算、観測/出力/再開要求 | 元残差、各consumer、実checkpoint |

実装開始には、辞書の全型が閉じていること、各variantの最小例、全種類に対する型/未知key/欠落/相互条件/参照の負例、巨大ID/十進低位保持、入力不変/所有/資源境界を独立に審査する。製品実装で試験するDoDは設計検査の表を使用する。例文のparse成功だけでschema適合を証明しない。

本書の設計審査、schema実装、SI/snapshot、物理prepare、全機能受入は別々に状態を記録する。設計が完了してもIssue #22全体、IM06-0、実務適用の完了を主張しない。
