# PR 128：外周辺と試験のoptional参照の訂正

2026-11-18。[レビュー指摘](https://github.com/YUKIKEDA/UVLAS/pull/128#discussion_r4609700784)が対象としたheadは 0bb20b52aca3c42e032babaaa3e04705d70e6336。この記録は、そのheadの分割実装に対する訂正を扱う。既存の製品検証manifestと凍結ログは変更しない。

## 指摘と修正

[Linux Debug job](https://github.com/YUKIKEDA/UVLAS/actions/runs/44233312914/job/134563402315)では、製品1箇所と試験6箇所にbugprone-unchecked-optional-accessが出た。Windowsで同じ診断が出なかったことを、Linuxでも検査済みとする根拠にはできなかった。

外周辺の位置traceを作る処理は、左右の隣接片を使用箇所で順に検査して値を取得する。どちらもなければ既存のGeometryInvalidを返す。生成側が左右いずれか一方の存在を保証していても、その保証だけで空のoptionalを参照しない。

試験側は、存在を通常のifで確かめ、取得した値をCatch2の比較へ渡す。fixtureのgraphもrequestを作る時点で確認する。比較する数値、許容差、独立参照、分割順序、資源上限は変更していない。再発条件と手順を既存のC++規約へ追記した。

## 対象検証

Windows 11、MSVC 19.51.36260、C++23、既存のout/build/issue125を使い、次の2 targetだけをDebug/Releaseで再構築した。

- uvlas_test_structure_pg_integration_partition
- uvlas_test_structure_pg_integration_partition_resources

ctest --test-dir out/build/issue125 -C Debug -R '^structure.pg_integration_partition(_resources)?$' --output-on-failureは2/2成功、3.63秒。Releaseの同じ選択は2/2成功、0.54秒。最初のDebug buildでは追加したローカル変数名が同じ試験内の既存名と重複して失敗したため、名前を訂正してから両構成を検証した。

正常側は7ケース・30,531 assertions、資源側は5ケース・168 assertions。存在確認のREQUIREを通常の失敗分岐へ移したため、旧記録よりそれぞれ2件と1件減っている。数値比較の検査は維持した。

clang-tidy 21.1.1は、Debugの重複を除いた既存専用DBで製品・正常試験・資源試験の3翻訳単位を検査し、exit 0、利用者コードの指摘0。変更4ファイルのclang-format 21.1.2とdiffの空白検査も通過した。

最初の修正commit（ef60053）時点では、修正後のLinux CIは未実行だった。レビューが求めるLinux静的jobの完走はpush後に追跡し、対象commitとjob結果を同じPRの返信へ残す。既存の数値認定・実務受入の範囲は広げない。

## 旧証拠との対応

旧headの必要なsourceとartifactをgit archiveで別directoryへ展開し、元のmanifestをそこで照合した。75ファイル・2 JSON、整合検査のエラー0を確認した。現在のtreeへ旧manifestを直接適用すると、今回訂正した4 sourceのhashは一致しない。これは元sourceに結び付く旧記録を保持した結果であり、旧hashを現在値へ差し替えない。[source bridge](im01-partition-product/optional-guard-bridge.json)に4 sourceのbefore/after、現在の全source hash、保持した旧artifactのhashと対象検証を記録した。

## Linux再検証と第2修正

ef60053の[Linux Debug job](https://github.com/YUKIKEDA/UVLAS/actions/runs/44249168147/job/134615450703)はclang-tidyで失敗した。旧7診断のうち5件は消えたが、正常試験が辺メンバーを再参照する2箇所が残った。Windowsの静的検査だけで完了とはしなかった。

第2修正のsourceは81e9f609f9f8cf1af04b03b0c371860b3349ad04。左右のoptional自体をif初期化でconstの局所値へコピーし、value_orへcore::Index{-1}を渡して範囲を比較する。-1は必ずpiece範囲外なので有効なpieceを捏造せず、存在する値の比較は維持する。対象2箇所のoperator*は除去した。

変更した正常試験だけを再検証し、Windows Debugは1/1成功、3.43秒、Releaseは1/1成功、0.44秒。正常側は引き続き7ケース・30,531 assertions。変更した正常TUだけのclang-tidy 21.1.1もexit 0。製品・資源試験のsourceは最初の修正から変えておらず、それらのWindows試験と静的検査は重複して走らせていない。第2修正後のLinux jobはpush後に完走まで追跡する。
