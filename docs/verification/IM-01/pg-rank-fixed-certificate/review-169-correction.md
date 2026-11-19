# PR #169：結果を全件合格後に一度だけ公開する

2026-11-19。指摘 `discussion_r4623898683`、対象head `a97961a3ce88da05c942e64c417570380ff07867`。対象はrunnerの結果公開であり、数値証明器や凍結12件の数値結果は変更しない。

## 欠陥と訂正

旧入口は1件が成功するたびに `--output` を上書きしていた。公開処理だけを切り分けたprobeで2件目の失敗を注入すると、既存出力がある場合もない場合も、有効な形式の1件JSONが残った。既存の完成結果のbyteは失われた。

現行の `run-verifier.py.txt` を修正した。12件をメモリへ集め、status/heap/rank/rows/candidatesを明示的に検査する。条件は `python -O` でも有効である。全件が成功してから隣接する一時fileへ一度書き込み、flush/fsyncして閉じた後、`os.replace` を1回だけ行う。書込や置換が失敗したときは自分のtempを片付ける。途中の進捗表示はstdoutへ継続する。

旧runnerと旧manifestは `review-169-original-runner.py.txt` と `review-169-original-manifest.json.txt` へbyte保全した。既存manifestは旧runnerの参照先だけをその原本コピーへ結び直し、hash値は変えていない。新runnerの検査結果を、旧数値実行の成功証拠へ転記していない。完成した結果を置き換えるrunnerの規則を、evidence規約へ追加した。

## 対象検証

`review-169-publication-probe.py.txt` は数値transportをstub化し、実際の現行runner入口とfilesystemで公開処理を試す。rankや容量の再認定には使わない。旧欠陥2例を先に再現した後、通常と `python -O` の各20系列を実行した。

- 既存出力あり/なしの両方で、12件成功時だけ隣接tempから一度置換し、tempは閉じられている。
- 2件目の失敗、status/heap/rank/rows/candidatesの不一致を拒否する。
- 書込途中、fsync、replaceの各故障を注入し、既存のbyteを保持する。元々出力がなければ部分結果を作らない。
- 全系列でtempのfileとhandleが残らない。

結果は通常版 `review-169-publication-result.json`、最適化版 `review-169-publication-optimized-result.json` に保存した。数値証明器・係数入力・全12結果のbyteは不変であり、重い全rank再実行は行っていない。独立レビュー担当が現行入口を通読し、従来の12件の順序、実行中の出力書込なし、close後の1回置換、最適化時も有効な条件を確認した。

```text
python -B review-169-publication-probe.py.txt --legacy review-169-original-runner.py.txt --runner run-verifier.py.txt --output <新規通常結果>
python -B -O review-169-publication-probe.py.txt --legacy review-169-original-runner.py.txt --runner run-verifier.py.txt --output <新規最適化結果>
```

PRとinline threadは利用者の確認まで開いたままにする。製品開始や全call資源契約について、親#140の未達条件は変わらない。
