# Linuxレビューの実行確認とtidy訂正

指定run 44335012542のheadが0597204420ec324e70badfb00ced1b11f215d397であることをAPIで確認し、Linux Debug・Release・Sanitizerの3jobだけを追跡した。3構成で新しい2試験が実行され、Passedになった。表の値はその時点の履歴で、今回の訂正後headに対するLinux成功ではない。

| 構成 | 実CTest | normal / resources | job全体 |
|---|---|---|---|
| [Debug](https://github.com/YUKIKEDA/UVLAS/actions/runs/44335012542/job/134909296362) | 92/92成功 | 0.72s / 0.01s | 後続tidyの8診断で失敗 |
| [Release](https://github.com/YUKIKEDA/UVLAS/actions/runs/44335012542/job/134909296429) | 93/93成功 | 0.17s / 0.00s | 成功 |
| [Sanitizer](https://github.com/YUKIKEDA/UVLAS/actions/runs/44335012542/job/134909296387) | 92/92成功 | 3.23s / 0.06s | 成功 |

Debugの実ログで、owner optionalの未検査アクセス5件、Abort::code未初期化、試験入力optionalの未検査変更、stop_sourceのconst指摘を確認した。ownerを局所const参照へ保持し、存在を明示判定してから参照する。公開getterは不正状態をInvalidInputで返し、private buildは既存の値failure捕捉へ渡す。同類検索でcheck入口とrelease_keyにも同じ判定を適用し、後者のoptional抑制を取り除いた。Abortのdefault codeをInvalidInputで初期化した。公開資格の所有・成功結果・演算式・予算は変わらない。

試験はoptionalそのものを変更せず、ローカルarrayを変更してから格納する。const stop_sourceはMSVC14.51の非const request_stopでbuildに失敗したため、既存Aero試験と同じ理由コメントとmisc-const-correctnessの1行抑制を採用した。Linux libstdc++はこのメソッドをconstとして宣言する。この差を他のoptional警告の抑制へ広げていない。

対象Debug2targetの最終buildとCTest 2/2は成功した。通常clang-tidyを新規6TUに適用した初回は、製品4TUに診断がなく、testのMSVC const差だけでexit1になった。訂正したnormal test単独と、先行診断の状態が残ったresource test単独を最終確認し、双方exit0だった。先行する失敗出力も保持した。format/diff-check成功。検証出力は個人パスと行末空白を除いたLFの表示記録である。

数値方式・固定参照・全組合せは変わらないため、既存のLinux3構成の実行結果と隔離変異は履歴として保持し、新たな全体試験や変異を増やさない。Windows parserは変更しない。今回push後のCIは待機しない。旧仕様・証拠は書き換えず、この追補とsource SHA256だけを追加した。
