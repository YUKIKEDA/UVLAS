# PG局所rank：全候補の有限体discoveryを固定領域で生成する

2026-11-19。子Issue #183、IM01-7、親#140。`solver-rd` と `solver-design-contracts` を適用する。[固定証明器](../pg-rank-fixed-certificate/verification.md)と[2列の固定持上げ](../pg-rank-fixed-lift/verification.md)に先行する、全候補の生成を実行した。

## 判定と境界

4個の採用31-bit素数と、n1/2/4・平坦/四次・全背景/U形の12ケースを組み合わせた48系列すべてで、固定領域のdiscoveryが完了した。生成した選択列・全剰余座標は凍結exact certificateと一致し、pivot・lower/upperは旧動的eliminateとも一致した。固定整数のlimbを剰余へ変換するところから照合の終了まで、追加operator newの試行は0回だった。

有限体の生成側は、この固定配列構成で引き継げると判断する。元係数の展開、exact入力の保持、旧動的oracleの作成は監査区間の外である。有理持上げ・exact証明までの全候補接続、ResourcePlanによる全準備、取消・故障・並列の製品I/F、辞書とpartitionを同時保持する48 GiBは親#140に残す。有限体だけで一般入力のexact rankを認定しない。

## 入力と比較先

凍結6係数原本から、要求n/shape/domain/row_orderで参照を一意に選び、headerとSHAを照合した。全12certificateのサイズとSHAも照合する。対象素数は2,147,483,647 / 2,147,483,629 / 2,147,483,587 / 2,147,483,579。

元候補はG→原列→Jの順を保つ。監査前に、凍結比較器の入力展開から得た有理係数を、分子・分母各1,536 bitの固定整数へ保持する。参照座標は各2,048 bit。ここでの動的な展開を無heapとは呼ばない。

生成核は期待selectedや座標を選択に使わない。すべての候補を消去した後、暗黙の零も含めた全columns×rankを、凍結exact座標の同じ素数での剰余へ照合する。さらに、別に実行した旧動的eliminateの全pivotとlower/upperの位置・値へ照合する。入力展開は共有するが、核は疎vectorから固定配列へ移し、元有理数からの剰余化もBoostの除算からlimb走査へ分けた。

## 固定配列と手順

行数をm、候補数をN、原列数をF、容量上のrankを $$ R=\min(m,F) $$ とする。実rankをrとし、rを予約量の根拠へ転用しない。

- work[m]、rhs/inverse/selected/pivot[R]、座標[N×R]。
- 正規化したreducedの行・値[m×R]。
- upperは $$ R(R+1)/2 $$、lowerは $$ R(R-1)/2 $$ 個の位置・値。
- 各offset[R+1]とpivotの逆引き[m]。

sizeof、実input capacity、各shapeを検査付きuint64で積み上げ、scratchを作る前に計画する。作成後のcapacityも照合する。bodyでresize/pushをしない。payloadは比較用oracleも含む固定消費区間の保持量で、動的準備物や製品全体の保持量とは区別する。

各候補について、作業列を零にし、元の固定整数の分子/分母を同じ素数で剰余化する。分母がその素数で非可逆なら拒否する。既存pivot順に消去し、残りが非零なら最初の非零行を選ぶ。Jで新しい独立列が必要になった場合は、basisを拡張せず拒否する。従属候補はupperで後退代入して座標を得る。最後にreducedからlowerを作り、独立な参照へ照合する。

剰余化は32-bit subword単位で行う。累積値はp未満なので、次のshift/or/modは64 bitに収まる。有限体積も $$ (p-1)^2<2^{62} $$ である。監査区間では多倍長の加減乗除やgcdを使わず、固定backendからlimbと符号を読む。

## word課金と訪問のモデル

`word_charge` は定義したnative数値演算の課金量である。積と剰余を2、正規化差を高々2、limb抽出とshift/or/modを5として先に課金する。逆元はFermatの反復を使い、1回の反復で最大2積とbit判定/shift、最大31反復なので上界186。これはCPU命令数や時間の単位ではない。

`visits` は、明示した配列fillの要素、係数/limb、指数反復、消去・pivot探索・後退代入、lower抽出、参照座標・因子の走査点で加算する。一部の外側添字制御、最終条件判定、report作成などは含まない。wordは数値算術、visitは走査を制限するモデルとして扱い、全C++命令の数と呼ばない。

元係数の項数をI、参照の非零座標数をE、両方の分子・分母の32-bit subword総数をL、 $$ T=R(R-1)/2 $$ とする。wordの有限上界は

$$
W=5L+2(I+E)+186(I+E+R)+2(I+E)
 +4(NmR+NT)+2(mR+NR).
$$

各項は順に、limbと符号、逆元、分子と逆元の積、消去/後退更新、正規化/座標生成を囲う。visit上界は次の段階の和である。

| 段階 | visit上界 |
|---|---:|
| 座標とoffsetのfill | $$ NR+3(R+1) $$ |
| 候補・work/rhsの初期化 | $$ N+Nm+NR $$ |
| 元項・limb・逆元 | $$ I+2L+31(I+E+R) $$ |
| 前進消去 | $$ NR+NmR $$ |
| pivot探索 | $$ Nm $$ |
| 正規化とupper保持 | $$ mR+T $$ |
| 後退代入 | $$ NR+NT $$ |
| lowerの抽出 | $$ m+2R+mR $$ |
| 全参照座標 | $$ NR+E $$ |
| 選択列/pivot/offset/因子 | $$ R^2+4R+2 $$ |

有限上界とruntime capは役割を分ける。各演算・訪問の前にcapを照合し、超えれば未完として拒否する。構造上界以下のcapでも、そのcap内で完了すれば受理できる。構造上界分の値を演算ごとに新規所有したり、全入力の完走を保証したりしない。

今回のruntime capはword50億、visit30億。n4四次全背景の構造上界はword114,947,656,972 / visit28,806,656,384であり、これは観測費用ではない。観測値を使ってこの上界を縮めていない。

## 結果

全48系列のrankは凍結exact certificateと一致した。

| 面・領域 | n1 | n2 | n4 |
|---|---:|---:|---:|
| 平坦・全背景 | 109 | 224 | 490 |
| 平坦・U形 | 77 | 149 | 316 |
| 四次・全背景 | 183 | 368 | 774 |
| 四次・U形 | 128 | 242 | 498 |

n4四次の範囲は次のとおり。4素数の各実行で全座標と因子を照合した。

| 項目 | 全背景 | U形 |
|---|---:|---:|
| 固定保持payload B | 1,141,280,752 | 347,112,052 |
| 座標比較数 | 1,157,904 | 745,008 |
| lower/upperの値比較数 | 57,431 | 29,410 |
| word課金 | 659,523,210〜676,921,298 | 150,102,448〜154,383,120 |
| visits | 218,598,473 | 54,834,336 |
| 追加new試行 | 0 | 0 |
| 監査区間ms | 555〜613 | 133〜140 |

Windows x64、MSVC19.51.36260 Release、Boost1.92、32-bit limb、stack予約64 MiB。別作業と並行した診断値であり、性能比較や全ソルバーの予測ではない。入力展開を含むプロセスpeakは全背景1,557,385,216 B、U形526,331,904 Bだった。動的oracleと入力展開は監査区間の外なので、この結果を全準備の無heapや48 GiB認定へ読み替えない。

## primitive参照と負対照

4素数すべてで、 $$ -3/7 $$、正の $$ 10^{-400} $$、7の逆元、最大剰余同士の積、0−1を、Pythonの独立整数計算へ照合した。追加newは0回。手書き2行5列はrank2で、微小値の座標も正しく保持した。

payload/word/visitの等号・1不足、rank cap、非可逆分母、Jに新しい独立成分がある場合、参照改竄、寸法積overflowを拒否し、失敗後も同じconst入力で成功した。

独立レビューで、当初のbad-prime対照は先頭列だけをp1倍しており、全体の有限体rankが2のままと指摘された。単なるReferenceMismatchではrank不足の拒否を検証したことにならない。row0の全係数をp1倍する対照へ直した。有理数上は非零minorがp1でrank2、有限体ではrow0が零になりrank1となることをPythonでも確認し、実際のScratch.rank==1と拒否を検査した。

修正はcontrols関数内だけで、48実ケースの生成・比較・計数経路は不変であることをsourceの差分と領域hashで確認した。全48の観測行は保持し、小対照だけ再実行した。由来は `control-correction.json` に記録した。時間とpeakは、その48実行時の診断値である。

## 再現と残条件

`fixed-discovery.cpp.txt` をC++としてMSVC `/TP /std:c++latest /O2 /EHsc /MD /DNOMINMAX`、vcpkgのBoost1.92 include、`/STACK:67108864` でビルドする。`run-discovery.py.txt` に `--inputs`、`--certificates`、`--discovery`、`--capacity-producer`、`--study`、`--certificate-index`、`--output` を渡す。全roster・全参照・正常条件が揃ってからだけ結果を原子的に公開する。各ケースは300秒で停止し、未知のstatusやtimeoutを合格にしない。

独立レビュー担当は、生成が期待selected/座標を使わないこと、全候補順、暗黙零を含む全座標照合、pivot/LU照合、48roster・計数・limb・監査境界を確認した。bad-prime対照の訂正も反映した。重い48実行を重複していない。

この子Issueは有限体の全候補生成の成立判断まで。親IM01-7/#140には、入力/LCMの固定化、全候補の有理持上げ・exact証明との一体化、製品の全call資源契約、ResourcePlan、取消/故障/並列、同時保持48 GiBを残す。全pipeline・製品rank・R-02実用受入は未完である。
