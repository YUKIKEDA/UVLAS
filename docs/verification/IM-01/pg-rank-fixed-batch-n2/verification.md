# PG固定バッチ：n2四次・全背景の全730列

2026-11-20。子Issue #201、親#140、IM01-7、SH-D03。`solver-rd` と `natural-japanese` を使用した。[n1全列試作](../pg-rank-fixed-batch/verification.md)の核と運用を保ち、幅と規模が大きい1ケースへ伸ばした。方式・数値精度・資源枠を変えず、全730列でのアルゴリズム成立を確認した。製品開始や実務受入は認定しない。

## 同じ核で幅の大きい全列を通した

最初の1従属列は5008活動行・368座標を検査し、監査1 ms、固定保持661,037,288 B、追加new0で成功した。その後の全列実行も、既存の8192bit・128桁・定義イベント20億・48 GiB以内で完了した。各processのtimeoutは従来どおり60秒である。

| 項目 | n2四次・全背景 |
|---|---:|
| 全活動行 / 候補 / rank | 5,008 / 730 / 368 |
| 選択列の単位座標 / 従属列の持上げ | 368 / 362 |
| 全行のexact零残差検査 | 3,655,840 |
| 独立certificateへのexact座標比較 | 268,640 |
| 定義イベント | 204,815,279 |
| 固定保持（計画と実capacity一致）B | 1,536,520,048 |
| 復元座標の最大bit / 全行証明の中間幅上界 | 1,411 / 1,592 |
| batchの追加operator new試行 | 0 |
| 監査時間 / プロセスpeak B | 1,941 ms / 1,523,195,904 |

従属列の成功桁数は、1桁213列、2桁20列、4桁30列、8桁57列、64桁10列、128桁32列だった。既存n1の最大695bitに対し1411bitまで全列で処理した。切断ケースの最大座標708bitよりも、この有理復元幅の負荷を優先してfullを選んだ。観測からcapや固定slotを縮めていない。

## 旧sourceは変更せず、mainの実験selectorだけを派生した

`derive-driver.py.txt` は旧C++原本のSHA-256を固定し、mainのn=1 guardと診断文をn=2へ、ちょうど1箇所だけ置き換える。逆置換して原本とbyte一致することも確認する。有限体factor準備、lift、全行証明、controls、監査、資源計画、full領域のselectorは全て同じbyteである。生成sourceのSHAとサイズは `driver.json` に保存した。原本の末尾へ1 byte足した負対照を拒否し、別出力の生成も同じSHAになった。

旧runnerのcontrols検査、列別lift診断の型/範囲検査、invokeのtimeout、atomic公開をimportして再利用する。n2ケース用の入口が、n2/quartic/full/sparseの要求から凍結study・certificateを先に一意選定し、header/SHA、basisを照合する。全730列の順序と選択bitmap、各列の全5008行・368座標、桁・幅・イベントの型と範囲を確認した後だけ結果を公開する。

記録量を抑えるため、`results.json` は集計、選択ordinal、桁分布、最大幅、全列diagnosticsのcanonical JSONのSHA/byte数を保存する。実行時は省略せず全列・全行・全座標を検査している。列diagnostics全文は同じrunnerで再生成でき、各列の値には時間を含まない。canonical化はkey順を固定し、区切りをcommaとcolonだけにしたUTF-8 JSONである。

## 新caseの拒否と出力保持を軽量に検査した

新caseのgate probeは、実際の捕捉結果を数値transportのstubとして使った。通常と `python -O` で各38系列を通し、正常、行/座標不足、誤選択、別caseの行数、診断欠落・負値・float/bool、総整数積のfloat、時間/peakの誤型・負値/零、heap試行、timeout、置換失敗、平坦入力との交換を、旧出力あり/なしの両方で検査した。拒否時に旧byteまたは未作成状態を保ち、tempを残さなかった。別入力との交換は数値transport前に拒否した。 暗黙のraw書込みは行わず、旧raw名を正式出力に指定した場合も置換失敗で旧byteを保った。

共通controlsは今回のbinaryから取得し、非pivot行のsquare成功後の拒否、座標変異、未来basis、参照不一致、資源境界、失敗後の同scratch再成功という既存の必須観測を、型と値まで旧validatorで検査した。共通probeやn1/48有限体系列/n4最悪2列/無関係suiteを重複実行していない。

## 再現と未認定範囲

1. 凍結exporter/certificate producerの原本を準備する。n2原本とfull certificateのSHAは索引と照合する。
2. `derive-driver.py.txt --base <旧fixed-batch.cpp.txt> --output <build先のcpp> --record <driver.json>` で生成する。MSVC `/TP /std:c++latest /O2 /EHsc /MD /DNOMINMAX`、固定Boost include、旧sourceのdirectoryをincludeへ追加し、`/STACK:67108864` でビルドする。
3. `run-case.py.txt` へ `--inputs --certificates --batch --base-runner --capacity-producer --study --certificate-index --driver-record --output` を渡す。最初の1列、全列の順で検査し、合格後だけ結果を公開する。公開runnerの出力は指定した完成結果1個だけとする。raw diagnosticsが必要なら、同じbinaryへ係数原本/certificate/`all` を直接渡し、stdoutを別途捕捉する。
4. `gate-probe.py.txt --root <repo> --output <別出力>` を通常と `-O` で実行する。これは捕捉済み結果の公開gateだけを検査する。

Windows x64、MSVC19.51.36260 Release、Boost1.92、32-bit limb、stack64 MiB。計測値は並行負荷下の診断で、性能比較ではない。イベントは既存の定義した課金点であり、全call/全limb/CPU命令の上界ではない。

前回と同じく、動的Input/Matrix/factor/LCMと全target Dataの複製準備は監査外である。全Dataを準備した後、scratch構築前に48 GiBを照合する。この固定保持量に監査前の動的Input/Matrix/factorの生存分は含まず、全準備のResourcePlan入場や全プロセス48 GiBを認定しない。複製保持方式を一般製品へ採用しない。

IM01-7/#140には、n4・切断領域の全候補接続、全準備の固定化/共有、全call/limb課金、ResourcePlan/取消/故障/並列、辞書/分割との同時48 GiBと製品開始を残す。親Issueは開いたままとする。旧source・結果・採用仕様・製品コードを変更していない。
