# PR #204：派生原本の保護と添付来歴の事前検査

2026-11-20。レビュー対象headは `651123f74f3b8bd35a75aa711d51b24f93f763e4`。対象は2つのPython入口で、数値核、730候補の結果、元の検証記録は変更しない。

## 旧headで2件を再現した

派生器へbase/output、base/record、output/recordを同じファイルとして渡した。3組とも旧入口は終了code0で既存ファイルを書き換えた。base/outputは原本をn2 sourceへ、base/recordは原本をJSONへ、output/recordは生成sourceをJSONへ上書きした。

旧runnerへ空のdriver記録を渡すと、controls/first/allの3回を呼んだ後、空の来歴を添付した完成結果を公開した。この再現の数値transportは既存実行の捕捉値を返すstubであり、730列を再実行していない。

## 現行入口の修正

派生器は全3パスをresolveし、各組のパス一致と、既存ファイルのsamefileを確認する。相対表記とhardlinkで別名になった同じ実体も、全書込み前にValueErrorで拒否する。別々の3パスでは原本を保持して、従来と同じ派生SHA/byte数と記録を生成する。

runnerは生成記録がobjectであることと、6必須項目の型・値を、計算開始前に固定値へ照合する。

| 項目 | 期待 |
|---|---|
| schema | `uvlas.fixed-batch-n2-driver/1` |
| 原本SHA | `f46f8dbe506c1efa7e834c298b64dc1329c0d4b6b96ff1cb893ade44ce7cdb63` |
| 派生SHA | `3516c0253b60e06e73d14c9e3c9633767e11cbb9a8636160df0c6eaab33a32fe` |
| 生成byte数 | intの38,327 |
| 変更範囲 | `main n guard and diagnostic only` |
| kernel保持 | boolのTrue |

floatやboolを生成サイズのintとして扱わない。検査後の同じobjectを完成結果へ添付し、公開前の読み直しをなくした。この検査は指定された生成履歴の同定であり、実行ファイルとコンパイルsourceの暗号学的な結合を新たに認定するものではない。既存のbinary SHA記録と数値証明の範囲は変えない。

## 対象検証

通常と `python -O` の両方で、次を確認した。

- 旧headの3組の破壊と空来歴の公開を再現。
- 新派生器の3組×直接/相対/hardlinkの9衝突を拒否し、全既存ファイルのbyteを保持。
- 別パスで正常生成し、元原本と派生SHA/byte数を保持。
- driver正常、空object、配列、null、不正JSON、6項目の欠落/不一致/誤型という23種を、旧出力あり/なしで計46系列確認。
- 不適合記録は全transport呼出し0、旧結果byteまたは未作成状態を保持、temp0。
- 正常記録は3transportを通し、検査済みの6項目を完成結果へ保持。

失敗時の出力保持は実際の一時ファイルで調べ、数値transportだけをstubにした。数値C++、全730候補、古い共通probe、無関係な製品suiteは再実行していない。

## 凍結した旧証拠を保持する

元の2入口とmanifestをbyteコピーとして残した。旧manifestの引用は、その2つの同一byte原本へ結び直し、旧数値/結果/検証記録のhashを変えていない。今回の入口、再現probeと結果は別のreview manifestで固定する。

同じ種類の再発を防ぐ規則を `.cursor/rules/reference-data.mdc` へ追加した。3パス衝突の事前拒否と、結果へ添付する生成履歴の前段検査に範囲を限定する。
