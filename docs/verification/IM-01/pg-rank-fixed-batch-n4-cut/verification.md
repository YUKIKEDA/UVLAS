# PG固定バッチ：n4四次U切断の全1496列

2026-11-20。子Issue #206、親#140、IM01-7、SH-D03。`solver-rd` と `natural-japanese` を使用した。n1/n2全背景の全候補検査と、n4切断の最悪1列だけの持上げの間に残っていた穴を、同じ固定核と資源枠で閉じた。

## n4切断でも全候補が完了した

firstの1従属列は、全5634行・498座標を検査して成功した。固定保持1,695,539,200 B、監査1 ms、追加new0を確認してから全列へ進んだ。全1496列は、既存の8192bit・128桁・定義イベント20億・48 GiB、process60秒以内で完了した。

| 観測 | 全列結果 |
|---|---:|
| 全活動行 / 候補 / rank | 5,634 / 1,496 / 498 |
| 選択列 / 従属列 | 498 / 998 |
| 元式全行のexact零残差検査 | 8,428,464 |
| 独立certificateへのexact座標比較 | 745,008 |
| 定義イベント | 385,828,529 |
| 計画と実capacityが一致した固定保持B | 5,307,678,176 |
| 復元座標最大bit / 全行証明の事前中間幅上界 | 1,461 / 1,665 |
| 全batchの追加operator new試行 | 0 |
| 監査時間 / プロセスpeak B | 2,945 ms / 4,982,382,592 |

従属列の成功桁数は、1桁843列、2桁20列、4桁31列、8桁62列、64桁10列、128桁32列だった。基底の全順序、選択列の単位座標、従属列の未来basis拒否、元式の全行検査、全座標の独立比較は旧核のままである。capを観測値へ縮めていない。

## 派生は実験selectorだけ

旧C++原本の固定SHAを確認し、mainのn guard、診断文、cut領域のselectorだけを各1回置換する。逆置換して原本全byteへ戻ることも確認した。数値算術・因子準備・lift・全行証明・controls・監査・資源計画は同じbyteを使う。生成sourceのSHA/サイズは `driver.json` に固定した。

PR #204の修正を引き継ぎ、派生先/記録と原本/共通producerの全4パスを、書込み前に同一実体検査する。原本は一度だけ読み、その同じbyteから派生内容と来歴を作る。runnerはdriverの6必須項目を型と固定値へ照合してから計算する。

n4/quartic/cut/sparseを要求として凍結studyとcertificateを一意選定し、その後にheader/SHA・basisを確認する。各列の全5634行・498座標、選択bitmap、status、桁/幅/イベント/整数積、全体の時間/peakの型と範囲も検査する。全列が合格した場合だけ、指定された完成結果1個をatomic公開する。暗黙のraw出力はない。

## 検証と記録を分けた

公開resultは集計、選択ordinal、桁分布、最大幅、全列diagnosticsのcanonical JSONのSHA/byte数を保持する。実行時には全1496列を検査している。raw全文は保存しておらず、独立レビューが個別1496列やdigestを再計算したとは主張しない。同じbinaryへ元係数/certificate/`all`を渡してstdoutを捕捉すれば再生成できる。key順固定・空白なしUTF-8のcanonical化手順は旧 `compact` と同じである。

新caseの公開gateは、捕捉した集計と独立certificateのbasisから作ったschema stubで検査した。列別の桁/幅/イベントは合成値であり、数値oracleではない。通常と `python -O` の各28系列で、行/座標不足、誤選択、全背景行数との交換、診断欠落/誤型、空来歴/floatサイズ、timeout、置換失敗、平坦入力との交換を、旧結果あり/なしで拒否した。旧byteまたは未作成状態を保持し、temp0だった。別途、原本1 byte変更の拒否、派生SHA再現、原本/共通producerへの出力衝突拒否も確認した。

旧共通controlsは今回のbinaryから取得し、型と必須観測を同じvalidatorへ通した。旧n1/n2全列、48有限体系列、n4最悪2列、共通の広いprobe、製品suiteは再実行していない。

## 再現と限界

1. 凍結exporterとcertificate producerの原本を準備し、索引SHAと照合する。
2. `derive-driver.py.txt --base <旧cpp> --shared-producer <訂正済みn2派生器> --output <build先cpp> --record <driver.json>` で生成する。MSVC `/TP /std:c++latest /O2 /EHsc /MD /DNOMINMAX`、固定Boost include、旧source directoryのinclude、`/STACK:67108864`でビルドする。
3. `run-case.py.txt`へ `--inputs --certificates --batch --base-runner --capacity-producer --study --certificate-index --driver-record --output` を渡す。first→allの順で実行する。
4. `gate-probe.py.txt --root <repo> --output <別出力>`を通常と`-O`で実行する。これは公開条件のschema対照だけである。

Windows x64、MSVC19.51.36260 Release、Boost1.92、32-bit limb、stack64 MiB。時間/peakは並行負荷下の診断値で、性能比較ではない。定義イベントは全call/全limb/CPU命令の上界ではない。

全target Dataの動的複製、Input/Matrix/factor/LCM準備は監査外であり、この複製方式を一般製品へ採用しない。全Data準備後、scratch構築前に48 GiBと照合し、構築後に実capacityとの一致を確認した。固定保持に監査前の動的Input/Matrix/factorの生存分は含まず、全準備入場を認定しない。

この切断入力での全候補アルゴリズム成立を確認した。IM01-7/#140には、残る未実証ケース、全準備の固定化/共有、全call/limb、ResourcePlan/取消/故障/並列、辞書/分割との同時48 GiB、製品開始を残す。親Issueは開いたままとする。#204の教訓は同PRの訂正記録/規則へ反映済みであることを次branch前に振り返り、追加の方式/要求変更を行っていない。
