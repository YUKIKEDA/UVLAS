# PG辞書・共通分割・rank試作の同居確認

2026-11-18、Issue #140。四次n=4の辞書、U形2層の製品共通分割、切断域のrank試作を一つのWindows x64 Releaseプロセス内で順に実行した。分割の製品オブジェクトは[PR #149](https://github.com/YUKIKEDA/UVLAS/pull/149)の固定版 `6049a768a458abe935c9d9bad7bf3a09a084d7f2` を先にリンクした。辞書入力・rank試作は[比較実験](../pg-rank-resource-study/verification.md)と同じ。

辞書と分割には同じ48 GiB `ResourcePlan` を使い、両結果をrank試作が終わるまで保持した。分割の全64背景セルについて、活動判定をU形の3矩形との正面積判定と照合し、18活動セルが一致した。rank試作は切断域5,634行でrank 498を厳密に認定し、単独試験と候補数と再構成の整数件数が一致した。試作終了後にも同じ要求を使って保持した分割を検証した。

| 観測 | 値 |
|---|---:|
| 共通分割前の累積Plan確保 | 4,513,254,640 B |
| 共通分割の累積確保予測 | 14,655,590,207 B |
| 共通分割の累積確保実績 | 30,029,119 B |
| 共通分割後の累積Plan確保 | 4,543,283,759 B |
| 共通分割後のPlan残量 | 46,996,323,793 B |
| 同一プロセスのピークRSS | 384,491,520 B |

分割のBoost外部予測は19,978,141,584 B、辞書側は11,747,038,224 B。これらの見積もり、累積Plan確保、RSSは異なる指標である。rank試作は作業用の標準コンテナと多倍長領域を使う。測ったRSSにはその領域を含む。rank試作の確保を製品Planで事前予約し、実行中に制限する検証は #140 の残条件である。

再現には、比較実験のexporterを `export-coefficients.cpp` に複製し `coresident-driver.diff.txt` をunidiff-zeroとして適用する。rank比較用C++ TUを `main=rank_comparison_main` でコンパイルして、変更したexporter、固定分割オブジェクト、製品structure/core、Windows計測のpsapiとリンクする。固定版の `pg_partition_test_support.hpp` をinclude pathへ入れる。通常のexporterと同じ引数の `4 quartic <新規出力>` で実行する。記録は[result.json](result.json)。
