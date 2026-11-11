# 曲面対応ANCFシェルの比較と次の選定対象

版：0.1、2026-10-29。ST-01〜04の候補比較。初期曲面・積層・指定初期ひずみの要求を維持する。主比較対象を3833型へ移すが、製品要素の最終採用とは区別する。

最新の判断：[薄肉・歪み試験](records/3833-bending-distortion-experiment.md) で補正なし3833型の主シェル採用を見送り、[9節点・厚さ二次モデル](q9-thickness-quadratic-design.md) を主比較対象にした。以下はその判断に至る候補・式の記録として保持する。

## 1. 候補を同じ「高次要素」とまとめない

| 候補 | 節点・座標 | 確認できたこと | 採用前に残ること |
|---|---|---|---|
| 3423型・現ANS/EAS形 | 4節点、各位置＋厚さ勾配、24座標 | 補間・参照写像・質量を定式化。層別EASと曲面ANSに解析反例 | 既存のままは採らない。曲面の再現性を保つロッキング対策と接続を別に確立する必要 |
| 8節点・厚さ一次のQ8型 | 8節点、各位置＋厚さ勾配、48座標 | Hyldahlの学位論文に具体的な導出がある | 同論文は二重曲率の厚さロッキングと検証の不足を残す。3833と同一定式化ではない |
| Chrono 3833型 | 8節点、各位置＋厚さ一次・二次微分、72座標 | 固定版の補間を取得。厚さ二次項と整合する参照写像を持つ | 薄肉・二重曲率・歪み・層別初期ひずみ・応力・固有特性・接続を審査 |
| 3443型 | 4節点、各位置＋3方向勾配、48座標 | 完全な勾配表現を持つ | 既に指摘されているロッキングと、表面・接続の自由度対応が未解決。優先順位を下げる |

Hyldahl, *Rectangular shell elements based on the absolute nodal coordinate formulation* (2015)、[公開全文](https://ebooks.au.dk/aul/catalog/download/97/90/303?inline=1) の§16、式(16.44)–(16.45)は48座標の構成を示す。§21は二重曲率・厚さロッキングと追加検証の課題を残している。この論文を3833型の72座標モデルの検証済み証拠として流用しない。

## 2. CS-01：3833型の補間を固定版ソースと照合

[3833実装の固定版](https://github.com/projectchrono/chrono/blob/d26ed26d8344ae9aee4c1d26bf6d481076750b1a/src/chrono/fea/ChElementShellANCF_3833.cpp) の `Calc_Sxi_compact` と座標取得経路を照合した。ファイルを追加取得し、[ソースマニフェスト](../ancf-3423/records/chrono-shell-source-manifest.json) にURL・バイト数・SHA-256を追記した。

節点は四隅A,B,C,Dと辺中央E,F,G,H。面内の8節点serendipity補間を $N_i(\xi,\eta)$ とすると、ソースの係数は次の形へ整理できる。

$$
\mathbf x_h(\xi,\eta,z)=\sum_{i=1}^{8}N_i(\xi,\eta)
\left(\mathbf r_i+z\mathbf d_i+\frac{z^2}{2}\mathbf c_i\right).
$$

$$
\mathbf d_i=\left.\frac{\partial\mathbf x_h}{\partial z}\right|_i,
\qquad
\mathbf c_i=\left.\frac{\partial^2\mathbf x_h}{\partial z^2}\right|_i.
$$

$\mathbf r_i$ はm、 $\mathbf d_i$ は無次元、 $\mathbf c_i$ は $\mathrm{m}^{-1}$ 。2本の勾配を2本の独立回転ベクトルや、2方向の面内勾配へ読み替えない。

$$
\mathbf q_e=
[\mathbf r_A^{\mathsf T},\mathbf d_A^{\mathsf T},\mathbf c_A^{\mathsf T},\ldots,
\mathbf r_H^{\mathsf T},\mathbf d_H^{\mathsf T},\mathbf c_H^{\mathsf T}]^{\mathsf T}
\in\mathbb R^{72}.
$$

ソースは層内の親厚さ座標から、要素全厚さ基準の $z$ へ変換して係数を評価する。層ごとに $z=0$ を置き直す補間ではない。これはEASの層内モード座標を零平均にする決定とは役割が異なる。

## 3. CS-02：8節点serendipityの具体形

隅点の符号を $(\xi_i,\eta_i)\in\{-1,+1\}^2$ として、

$$
N_i=\frac14(1+\xi_i\xi)(1+\eta_i\eta)
(\xi_i\xi+\eta_i\eta-1).
$$

辺中央は、E=(0,−1)、F=(1,0)、G=(0,1)、H=(−1,0)の順に

$$
\begin{aligned}
N_E&=\frac12(1-\xi^2)(1-\eta),&
N_F&=\frac12(1+\xi)(1-\eta^2),\\
N_G&=\frac12(1-\xi^2)(1+\eta),&
N_H&=\frac12(1-\xi)(1-\eta^2).
\end{aligned}
$$

節点補間・分割の一・必要微分を独立に確認する。9節点Lagrange型と同じ多項式空間ではなく、節点数が増えただけで任意の曲面を厳密表現できるとは主張しない。

## 4. CS-03：曲面でのアフィン変形再現性

参照配置も同じ補間で定義する。

$$
\mathbf X_h=\sum_{i=1}^{8}N_i
\left(\mathbf r_{i0}+z\mathbf d_{i0}+\frac{z^2}{2}\mathbf c_{i0}\right).
$$

一定行列 $\mathbf F_*$ と並進 $\mathbf t_*$ によって全節点を

$$
\mathbf r_i=\mathbf F_*\mathbf r_{i0}+\mathbf t_*,\qquad
\mathbf d_i=\mathbf F_*\mathbf d_{i0},\qquad
\mathbf c_i=\mathbf F_*\mathbf c_{i0}
$$

と変換する。 $\sum_iN_i=1$ なら、

$$
\mathbf x_h=\mathbf F_*\mathbf X_h+\mathbf t_*,\qquad
\mathbf F=\mathbf x_{h,p}\mathbf X_{h,p}^{-1}=\mathbf F_*.
$$

これは適格な参照写像上での補間の代数的性質。曲面でも一定物理ひずみを表すための出発条件になる。ひずみを後からANS等で変更する場合、その変更後にも成立するかを別に確認する。

厚さ二次項が存在することだけでロッキングが除去される、という証明ではない。層材料・求積・境界条件・初期ひずみ・接線まで含めて審査する。

## 5. CS-04：次に閉じる採用条件

1. 3833の全必要微分、参照Jacobian、エネルギー、質量、内力・接線、材料軸を具体化する。
2. 同じ曲面のアフィン変形・剛体運動、層分割・初期ひずみ、質量特性を独立解析参照と照合する。
3. 平板と一重／二重曲率で、薄肉化・格子歪み・求積・自由要素の零モード・予荷重固有特性を確認する。
4. 生の積分経路と高速化経路、数値減衰と物理減衰、応力回復のモデル一致を確認する。
5. 72座標を含む勾配・厚さ二次座標の支持、シェル交差部材・梁接続、外部格子からの参照場生成を定義する。

局所接線の成分数は72²で、3423の24²の9倍。しかし、同じ精度に必要な格子や疎構造は異なるため、この比を総計算時間の倍率とはしない。製品の10万自由度・メモリ・フラッタ探索目標に対して別途容量を見積もる。

判断：曲面の主比較候補を3833型へ移す。48座標Q8型は別候補として残すが、既知課題を解消済みとは扱わない。新しい補間補正を即席で採用せず、原式・固定版・反例・受入条件を対応付けて選定する。

進捗：[3833型の数式仕様](ancf-3833-formulation.md) と [解析参照ケース](records/3833-reference-cases.md) を作成し、必要微分・参照写像・材料・初期ひずみ・内力・接線・求積を定義。[代数試験](records/3833-algebraic-experiment.md) で平坦自由要素のランクと剛体条件を確認した。質量と弾性の求積は分ける。曲面・薄肉・歪み・接続・製品適用の採用判定は未完了。
