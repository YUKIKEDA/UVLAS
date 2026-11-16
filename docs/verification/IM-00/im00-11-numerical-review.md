# IM00-11の数値部分の終了審査

2026-11-16。対象は `origin/main` の `1a3c98f`、Issue #15の完了条件1・2・6。レビューSkill、`spec-checks`、`numeric-tests`、`evidence` に従い、試験式、独立参照、固定予算、保存記録を通読した。数値試験・R&D生成・証拠checkerの一括実行はしていない。

**判定：数値部分の終了は保留。実4096系の製品認定が未達である。** 初稿は小系・一行の拒否witnessを確認した一方、精度契約・密求解契約がIM00-11へ残した全体系の認定を拾えていなかった。以下の既存基準の確認だけでは、Issue #15を終了できない。両OSの最終CI、型・状態・所有保証、次パッケージの入口も別に判定する。

## 終了を止める重大不足

[精度契約](../../specs/foundation/precision-contract.md)の「4096格子の失敗を実数値で再検証する」は、単一行の拒否witnessと全体系の受入を区別し、後者をIM00-11へ残す。[密求解契約](../../specs/foundation/dense-solve-contract.md)末尾と [im00-5.json](im00-5.json)・[im00-6.json](im00-6.json) にも同じ残作業がある。

既存fixtureは4096分割・k=.1のrow4095/pitch-quarter列で、元の右辺列L2ノルムに対する投影残差が $$ 10^{-11} $$ を超えることを示す。この一行で全系の投影不合格は示せるが、全4096行・全右辺列の保持解、原L2残差、必要観測の合格は示せない。元A/bのhashを維持した全系を製品UpdatedPreparedへ接続し、独立参照と元の全列L2基準で保持／投影を判定する追加試験が必要である。ここを満たすまで数値部分の終了を保留する。物理kernelや実務性能の先送りとは別の、既に約束された数値基盤の未達である。

追加試験の最小条件は、精度契約が指定する4096分割・k=.1の全3右辺列について、次を閉じることである。

- 丸めた有効係数 $$ A_e=A_0+v e_{last}^{T} $$ とBを凍結原本のhashへ照合する。転置を共役転置へ替えず、同じ実数baseと複素最終列更新を製品へ渡す。
- 製品の実返却値を独立に全行評価し、各列の $$ \|B_j-A_eX_j\|_2/\|B_j\|_2\le10^{-11} $$ と成分別後退誤差 $$ 10^{-12} $$ を検査する。成分別相対尺度やInfinity normを元の全L2基準へ読み替えない。参照の旧hi/lo bitや旧反復数を製品へ強制しない。
- 同じ固定高低参照の保持と投影を別に検査する。製品の再求解結果がBinary64でも元の全条件を満たせば適格とし、TwoWord昇格を結果の型だけで要求しない。
- 観測は3運動それぞれのF/Mと11固定圧力区間、計39個を両語から消費する。係数vectorを固定した代数consumerとして接続し、物理荷重kernelを実装したとは主張しない。原R&Dの5%圧力差・2%格子差を、求解器の移植差の許容値として使わない。
- 原本と新規生成値、全系参照と一行witness、保持と投影、実行版と資源条件を別に保存する。選択15行のDecimal監査や旧R&Dの成功flagだけで、新製品値の全行資格を代替しない。

## 完了条件1・2の既存根拠（4096全系は未達）

| 基準 | 読んだ試験と判定式 | 保存記録・成立範囲 |
|---|---|---|
| 密実求解・複素更新 | [dense_solve_test.cpp](../../../tests/numerics/dense_solve_test.cpp)。実3次系は独立100桁消去へ比較。複素rank-two更新は独立に列挙した有効行列と2右辺を使い、共役を取らない更新、全成分の既知解と原問題を固定 $$ 10^{-12} $$ 級で比較する | [im00-6.json](im00-6.json)、[im00-62.json](im00-62.json)。精密保持・不適格投影拒否を含む。物理的な後縁・荷重の組立認定ではない |
| 不定KKT | [sparse_solve_test.cpp](../../../tests/numerics/sparse_solve_test.cpp)。literal右辺2本と既知解を比較し、拘束 $$ q_0+q_1=0 $$ と反力を別consumerで固定 $$ 10^{-12} $$ に認定する | [im00-7.json](im00-7.json)。冗長拘束の拒否も別試験にある。実PG接続graphの認定ではない |
| 100k疎 | [sparse_capacity_test.cpp](../../../tests/numerics/sparse_capacity_test.cpp)。次数100000の三重対角、閉形式 $$ x_i=i+1 $$ と各右辺成分を独立に定め、全100000成分の最大相対差を $$ 10^{-12} $$ 以下とする | [im00-7.json](im00-7.json)。原系の許容差は `atol=1e-9, rtol=1e-12`、解の相対比較は別条件。合成容量であり、実graphのfill・RSS上限・実務時間の合格ではない |
| 負値・零を含む一般化全固有値 | [eigen_generalized_test.cpp](../../../tests/numerics/eigen_generalized_test.cpp)。解析根 $$ -1,4 $$ と非2冪尺度を検査。実返却値を有理数化し、元K/M残差・全Gramの行絶対値和を固定 $$ 10^{-12} $$ および報告上界と照合する。微小負・零・正は符号、零数、相対差を別に検査する | [im00-69.json](im00-69.json)、[im00-69-trace.json](im00-69-trace.json)。原M資格とscaled MのLLTを区別。不定Kや微小負根を消さない |
| 一般完全法・複素Arnoldi | [eigen_full_test.cpp](../../../tests/numerics/eigen_full_test.cpp)、[eigen_arnoldi_test.cpp](../../../tests/numerics/eigen_arnoldi_test.cpp)。解析的な実6根・複素6根、非共役更新、左右・部分空間を検査。完全法の元Schur残差とGramは独立多倍長計算、Arnoldiは元作用で再認定する | [im00-68.json](im00-68.json)、[im00-70.json](im00-70.json)、[尺度変換の追加訂正](im00-70-review-scaled-action.json)。訂正後の有理誤差球oracleと旧helperによる偽資格の負対照を確認。部分抽出を物理mode被覆と扱わない |
| 実PG／自由後流の凍結全履歴行列 | [foundation_reference_test.cpp](../../../tests/numerics/foundation_reference_test.cpp)は保存Pから定義した $$ A=I-P $$ を製品密求解へ渡し、独立Decimal参照解と各原残差を固定 $$ 10^{-12} $$ 級で検査。FullはPの全24根を独立SciPy参照24根へ、Arnoldiは明示 `nev=6,ncv=20` の6根を同参照と製品Fullへ照合する | [fixture導入記録](im00-11-fixture-preparation.json)と上記固有値記録。いずれの集合照合も参照の使用済みmaskで重複消費を禁じる。Fullは件数24、Arnoldiは件数6も必須なので、成功時は対象集合の全件一対一対応となる。製品によるPの再組立・時間積分・物理微分の検証ではない |
| EFT・低語消費・投影 | [precision_test.cpp](../../../tests/numerics/precision_test.cpp)。two-sum/two-productの両語和を独立算術で確認。極端な指数差は成分を厳密比較し、表現不能積とFTZ/DAZを拒否。元の4096格子投影失敗witnessも製品原残差へ接続する | [im00-5.json](im00-5.json)。全物理kernelのTwoWord消費は後続担当であり、ここで完了にしない |
| 有理rankの固定基準 | [foundation_precision_rank_test.cpp](../../../tests/numerics/foundation_precision_rank_test.cpp)。 $$ t=2^{-180} $$ の明示determinantとBinary64でのcollapseを区別し、LUの特異拒否へ接続。TwoWord保持では $$ b-Cx=0 $$ 、投影では残差 $$ t $$ を、尺度 $$ t $$ ・予算 $$ 10^{-12} $$ で拒否。零低語の成功対照もある | [im00-77.json](im00-77.json)。低語消失と拒否偽成功の隔離変異を検出。任意有理rank・TwoWord係数solver・PG拘束認定のAPIを実装したとは扱わない |
| 256bit区間平方根 | [linear_metrics_test.cpp](../../../tests/numerics/linear_metrics_test.cpp)。整数平方根から有理囲いを作り、元 $$ \sqrt{3/2} $$ の幅 $$ 2^{-256} $$ を保持。9複素／TwoWord入力と、6個のhalf成分を5指数・2要求の製品L2認定へ接続。元の平方不等式、囲いの両端包含、finite、固定幅上限、達成／未達の両側を検査する | [im00-77.json](im00-77.json)。製品幅は既存の `64eps/8minsub`、256bitは独立oracleの精度。平方根囲いの上下端を逆方向へ動かす変異と上界を内側へ動かす変異を検出 |
| SHA | [state_test.cpp](../../../tests/core/state_test.cpp)は製品model/payload fingerprintを独立Python hashlibの固定値と比較。[precise_arithmetic_test.cpp](../../../tests/toolchain/precise_arithmetic_test.cpp)は空文字・abc・56byteの標準3既知値を検査する | [im00-9.json](im00-9.json)。依存の既知値と製品のencoding・所有経路を区別する。新しいSHA APIは不要 |

GMRESも [gmres_solve_test.cpp](../../../tests/numerics/gmres_solve_test.cpp) と [im00-8.json](im00-8.json) に製品経路がある。固定前処理の再始動を用い、literal右辺と既知解の全成分を $$ 10^{-12} $$ で照合し、元残差を認定する。縮約Hessenberg解だけを成功にしない。凍結試作 `solver_contracts.cpp` の非線形chartを、この数値基盤gateのために先取りする必要はない。

## 完了条件6：参照と版の扱い

PG参照IDは `j1-pg-free-history-24-v1`、元recordは `numerical-foundation-adoption-audit`、SHA-256は `e1c9ca1f6103e0651c416bc77eea6280d434a1934970f4779e93c26124ff3ad8`。fixtureには3原artifactのID/hash/bundleと原生成sourceのhashが保存され、工具は索引・audit・原本を照合する。製品試験も保存artifactのhashを検査する。Decimal参照は製品求解を使わず120桁で生成し、160桁との照合と保存参照の原残差検査を分けている。

精密基準原本は `rd/foundation/numerics/experiments/cpp-numerics/precision_hash.cpp`、SHA-256は `6c512fa9951f143627b5ac1be847625bf15000e5123b6d5529372d21af9217ce`。rank・sqrtの参照保持と製品消費の範囲は [fixture契約](../../specs/foundation/reference-fixture-contract.md) と [im00-77.json](im00-77.json) に明記されている。

各保存記録のsource commit、候補byte、隔離変異、後段のコメント／タグ／const修飾の差は、その記録の時点のまま保持する。特に一般化の追跡タグ追加、Arnoldi尺度変換の訂正、rank変異後のconst修飾を、旧記録のsource-mapの書換えで現在版へ付け替えない。新終了recordは、これらの記録ID/hashと現在版の参照を別に所有する。

## 終了判断に残す境界

本審査は数値試験・参照・既存記録の内容の審査であり、最新CIの成功宣言ではない。CI採取と両OSの最終実行版の照合は統合担当へ渡す。実4096系の製品認定はIM00-11の残作業として解消する。物理PG/UVLMの組立、全履歴P/P*の製品生成、chart、全consumerのTwoWord対応、物理mode被覆、実験妥当性、製品規模性能、実務受入は該当する後続IMのgateに残す。
