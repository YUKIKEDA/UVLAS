# PG局所rank：全候補の固定持上げと元式のexact検査

2026-11-19。子Issue #193、親#140、IM01-7、SH-D03。`solver-rd`、`solver-design-contracts`、`natural-japanese` を使用した。[全候補の有限体探索](../pg-fixed-discovery/verification.md)と[最悪2列の固定持上げ](../pg-rank-fixed-lift/verification.md)を、実n1四次・全背景の全候補でつなぐ。

## 全383列で生成と元式の検査が完了した

選択183列は単位座標を作り、残る200列は同じ固定scratchを再利用してp-adic持上げと有理復元を行った。全383列について、生成座標を使って元の1252活動行を整数exact再構成し、全座標を独立certificateへ照合した。全batchを連続監査し、追加operator new試行は0回だった。

| 観測 | 結果 |
|---|---:|
| 候補・選択・従属列 | 383・183・200 |
| 全活動行の零残差検査 | 479,516 |
| 全座標の独立参照との比較 | 70,089 |
| 全行証明の整数積 | 479,374 |
| 定義したイベント | 27,754,320 |
| 計画と実capacityが一致した固定保持B | 399,889,676 |
| 元式を整数化した係数の最大bit | 90 |
| 観測した復元座標の最大bit | 695 |
| 全行証明の事前中間幅上界の最大bit | 857 |
| 全batchの追加new試行 | 0 |

従属列の成功桁数は、1桁133列、2桁4列、4桁25列、8桁16列、64桁22列だった。試行スケジュールは1/2/4/…/128桁のままであり、観測から固定8192bit profileや128桁capを縮めていない。選択列にも全行と全座標の検査を行い、従属列だけを全候補とは呼んでいない。

## 監査前の準備と固定区間

要求するn=1・quartic・full・sparseから、凍結studyとcertificate索引の記録を一意に選ぶ。その後に入力headerとSHA、certificate全413,605 BとSHAを検査する。入力やSHAの自己申告から別ケースを探さない。certificateのbasis列とpivot行を、有限体factorの実結果と比較する。

監査前に、元係数の展開、有限体factor、行LCM、全200従属targetのData、独立参照を準備する。この動的処理は無heapとは呼ばない。targetのDataは正方系を複製しており、一般製品向けの保持方式として採用しない。

Data準備の後、scratchを作る前に、全targetの実capacity、元式整数係数、独立参照、各vectorのobjectと予約領域、383件のreport、固定scratchと全行residualを計上する。checked uint64の和と積で48 GiB以内を確認し、構築後の全capacityと一致させる。固定保持量に、監査前に使った動的Input/Matrix/factorの生存分は含めていない。従って、全準備の入場や全プロセス48 GiBを認定した値ではない。

監査はbatchの先頭で有効にし、全候補の生成・再構成・独立比較を終えてから解除する。派生lift核は `owns_audit=false` で呼ぶため、列の途中で監査やnew計数を解除・初期化しない。bodyではtarget・scratch・report配列を拡張しない。scratchは全候補で同じものを使う。

`fixed-batch.cpp.txt` は旧固定liftから明示派生した。新しい全行処理とcaller管理の監査を加え、旧main、使わない旧certificate読取り・旧hand対照を除いた。元のsourceや結果は変更しない。checked/allocator void/et_off、Karatsuba cutoff 1,048,576、512bitのcurrentと8192bitの復元slotを維持する。

## 全活動行で何を証明するか

元候補係数をC、選択元をbとする。監査前に、活動行iの全候補係数の正分母LCMをeとし、整数係数を作る。

$$
E_{ij}=e_i C_{ij}.
$$

全Eを1536bitのchecked固定slotへ保持し、収まらなければ準備を拒否する。これはpivot行だけの正方系とは別の、元の全活動行の係数である。

復元結果を正分母の $$ n_k/q_k $$ とし、正方系検査で得た共通分母Dと整数化した座標yを使う。

$$
y_k=D n_k/q_k,\qquad
\sum_k E_{i,b_k}y_k-E_{ij}D=0\quad\text{（全活動行i）}.
$$

選択列では対応成分だけが1の単位座標、D=1を使う。検査前に、Eの最大bitをh、yの最大bitをv、rankをrとして、左辺と右辺の中間幅を別に囲う。

$$
B_{\mathrm{proof}}=
\max\left(h+v+\operatorname{bit\_width}(r)+2,
 h+\operatorname{bits}(D)+1\right)\le8192.
$$

上界を満たさない場合はProofBits、零でない行があればNonzeroResidualとし、成功へ置き換えない。全行検査の後、暗黙の零も含む全183座標を凍結certificateへexact比較する。従属列の未来basis非零、正方系の不成立、復元不足、p可除性の不成立も既存lift核で拒否する。有限体の候補選択だけで一般のexact rankを認定しない。

## 正方系で見えない不整合を拒否した

手書き系の先頭2行は、以前の独立Cramer参照を維持する。

$$
\begin{bmatrix}154&77\\77&231\\11&22\end{bmatrix}
x=\begin{bmatrix}-31\\72\\37/7\end{bmatrix},
\qquad x=(-3/7,5/11).
$$

runnerがPython FractionのCramer式から先頭2行の解を求め、第3行も37/7になることを別に計算する。3候補・rank2の正常対照は全9行・全6座標、83定義イベントで成功し、new0だった。

- 第3非pivot行の右辺だけを整数化後に変更した。正方系のliftは成功し、第3行の元式検査で拒否した。正方系だけの検査を全行証明へ誤認する変異を検出する。
- 生成座標を変更し、元式の零残差検査が拒否することを確認した。
- targetの選択元を未来列へ変更し、FutureBasisを観測した。
- 比較用座標を変更し、ReferenceMismatchを観測した。
- 固定payload等号と1 B不足、8192bitと8191bit、83イベントと82イベント、桁cap0を確認した。不足時の失敗を完走へ数えない。
- 拒否を繰り返した後、同じ入力と同じscratchで再成功し、イベント・整数積・全座標比較の件数が一致した。

これらは今回の試作に閉じた対照であり、製品ResourcePlan、取消、故障注入、再入・並列の受入を意味しない。

## 完成結果だけを公開する

runnerはcontrolsの版、件数のint型、非pivot拒否などのbool Trueを確認してから実列へ進む。最初の1従属列を検査し、全383列の順序、certificateの選択bitmap、各列の全1252行・183座標、成功status、new0、資源capを確認した場合だけ、隣接tempを閉じて一度の置換で結果を公開する。

軽量publication probeは実際の捕捉結果をstubとして返し、通常と `python -O` で各36系列を検査した。正常、旧controls、float/整数で偽装した型、非pivot拒否の欠落、失敗status、new試行、列/行/座標不足、順序重複、誤選択、桁型違い、未完了列、payload超過、timeout、置換失敗、平坦入力との交換を、出力あり/なしで試した。拒否時は旧byteまたは未作成状態を保ち、tempを残さなかった。このprobeは数値再実行ではない。

## 費用と再現の範囲

`defined_events` は既存liftの課金点に、候補、全行residualのfill要素、選択成分、元式の項、全行零比較、全座標比較を加えた回数である。fill/copyの一部、gcd/sqrt/Euclid内のlimb処理、bit観測、メタ比較、report、ループ制御などを含まない。全call・全limb・CPU命令数の上界ではない。batch上限20億の前に課金し、超過なら未完とする。

Windows x64、MSVC19.51.36260 Release、Boost1.92、32-bit limb、stack予約64 MiB。最初の1従属列は0 ms表示（ms分解能未満）で完了し、全batchは283 ms、プロセスpeakは397,213,696 Bだった。並行作業下の診断値であり、性能比較や全ソルバーの費用予測には使わない。最初の短い実行後に全batchを実行し、各プロセス60秒で停止する。

再現手順は次のとおり。

1. 凍結exporterとcertificate producerで係数6原本とcertificateを再生成するか、索引SHAと一致する既存原本を使う。
2. 新sourceをMSVC `/TP /std:c++latest /O2 /EHsc /MD /DNOMINMAX`、固定Boost include、`/STACK:67108864` でビルドする。監査profileは旧固定liftからincludeする。
3. `run-batch.py.txt` に `--inputs`、`--certificates`、`--batch`、`--capacity-producer`、`--study`、`--certificate-index`、`--output` を渡す。
4. `publication-probe.py.txt` に `--runner`、`--capture results.json`、`--root`、`--output` を渡し、通常と `-O` を別出力で実行する。

独立レビュー担当は、batchの連続監査、正方系から全行検査への接続、全383列のroster・全行・全座標集計、非pivot対照を通読した。重い数値の重複実行は行っていない。

## 採否と親Issueに残すもの

n1四次・全背景について、有限体で選んだ基底を入力に、全候補の固定領域での座標生成・復元と元式検査を接続した。この入力でのアルゴリズム成立を確認した。全targetのsquare Data複製を一般製品方式として採用せず、製品実装開始やR-02実務受入を認定しない。

IM01-7/#140には、n2/n4・切断領域の全候補接続、全準備の固定化と共有、全call/limb課金、ResourcePlan/取消/故障/並列、辞書・分割と同時保持する48 GiBを残す。これらを終えていない親Issueは開いたままにする。今回は数学・精度・製品範囲・48 GiBを変更していない。
