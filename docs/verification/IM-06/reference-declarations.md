# 外部参照宣言の所有収集

検査済みの `SchemaCase` から、8種類の `FileRef` を所有する宣言列へ移した。内部入口は [reference_declarations.hpp](../../../src/io/reference_declarations.hpp) の `collect_case_file_references` である。元IDの数値と字面、UTF-8パス、宣言digestの32 bytesと字面、参照位置を保持する。元文書、schema、資源Planを破棄しても読める。

この結果は宣言の収集だけを表す。パスの正規化、根や実体の取得、bytesの読取り、実digestの計算、依存関係の閉包、形式ごとの解釈は行わない。`ReferenceInput` や `ReferencedCase` への変換入口も持たない。これらは既存IM06-0の参照snapshotの後続であり、この部分の成功を全参照の取得成功へ読み替えない。

## 順序と出典

採用契約の辞書フィールド順、次に配列の添字順で訪問する。元JSONのkey順やIDの大小で並べ替えない。同じSourceIdを複数の場所で使っても、各出現を残す。SourceIdの整合性検査は先行schemaの資格を使い、ここで再び別のregistryを作らない。

各行のnodeは元のFileRef objectであり、JSON Pointerもそのobjectを指す。rootのSourceIdentityと、外部参照先のSourceIdは異なる役割で保持する。key順を反転した入力では宣言の意味と順序が一致し、元bytesのdigestとnode ordinalは新しい文書に対応する。後者を無理に同じ値へ揃えない。

## 予算と失敗

10個のquotaはすべて明示値で、0は無制限ではない。元文書の全node数とcontainer深さ、全出現数、元パスの合計、所有文字数、所有payload、scratch、累積payload、work、診断Pointerのbytesを分ける。Countと未知enumのgetterは値失敗になる。

走査は再帰や動的なDFS stackを使わない。schemaが所有する親情報、辞書、平坦な配列indexを読む。最初の無確保走査で件数、文字数、最大Pointer、workを計数し、全callを認定してから二巡目で所有領域へコピーする。各演算の加算・乗算と、配列を確保するIndexへの変換を検査する。

所有payloadを `P`、一巡目のworkを `W1`、コピーする文字数を `T` とする。`P` はData 1個、出現数分のoptional Row、文字配列の実ABI `sizeof` の和である。work上界は `3*W1 + 2*T + 2*P` とする。二巡のmetadata、Pointerの再走査と書込み、文字コピー、配列の初期化、行とheaderの書込みを含む。workは訪問、比較したkeyのbytes、添字の桁、Pointerの走査、コピー・初期化したpayload bytesを処理前に課す。CPU命令数や実時間の計測ではない。

scratchはContext、計数値、行の一時値、Pointerの添字buffer等の固定値payloadの上界である。compiler frameや例外ABI、process stack・RSSの実測値ではない。新しい動的確保は所有Data、行配列、文字配列の3点で、空の配列は確保しない。診断Pointerの予備領域も累積quotaとPlan残量の事前認定に含め、診断slotの有無で通常入力の受理を変えない。

失敗は親schemaと既存の結果を変えず、成功した確保のPlan課金を返金しない。診断の固定情報は借用しない。完全なPointerを作れないときは元failureを残し、断片を完全なパスとして返さない。開始時の取消、集計予算や確保の失敗など、特定nodeの出典を認定しない場面ではnodeとPointerを付けない。同じPlanまたはslotの再入・並行利用は呼出側の禁止条件であり、runtime再入検出を認定した意味ではない。

## 独立に確認する内容

手書きの正常入力は、Gmshの3宣言、Nastranの2宣言、prepared seed、観測参照、初期historyを一つのtransient caseに置く。役割、ID、元パス、digest、Pointerの期待列は製品から生成しない。さらに公開InputDocumentを別途作り、Pointerからたどった元nodeと照合する。8種類の具象FileRef型を辞書から列挙し、手書きの集合と過不足を比較する。

正常試験には、0件、同じIDの再宣言、先頭zero、ID 0とuint64最大値、JSON key逆順、配列の入替、元ownerとPlanの破棄、move後の状態を含めた。資源試験は全quotaの境界、Planの等号と1不足、全3確保点の故障、診断の予算不足と確保失敗、開始時と最後の確保中の取消、独立2呼出しを対象とする。失敗後は旧結果の全payloadと、先行schemaの191 recordにわたる既存の独立visitorの出力を照合する。

確保監査が数えるのは置換したC++ newの全familyである。例外ABIを含む全C heapの確保を観測したとは主張しない。独立変異では無変異controlのcompile・link・試験成功を先に確認し、roleの省略とdigestの破壊が正常期待に反することを確認する。結果の件数、実行環境、source hashは併設manifestと生ログを正本とする。

Windowsの対象Debug/Releaseだけをローカルで検証する。外部snapshotのOS adapter、数値変換、構造・空力の計算には変更を加えておらず、その無関係なsuiteは繰り返さない。LinuxはCIへ渡し、未実行の成功を主張しない。数値認定、製品規模の性能、実務受入の判定は、この非数値の内部部品の完了からは得られない。
