# PG-Ritz：一般参照面と実務力学の内部契約案

2026-11-05。採用前の候補仕様。柱面R4.1の成立範囲を拡張したとみなさない。一般曲面での精度・薄肉増幅・空間収束の判定を別に要求する。

## 参照幾何と厚さ

参照面を正則な区分C²写像 $$ \boldsymbol r_0(\xi^1,\xi^2) $$ とする。必要入力は位置、一次・二次微分、向き、物理厚さと一次微分である。形状のC²は要素内部で要求し、境界で必要な片側値を保持する。折れ角は同一の滑らかな面へ平均せずパッチ接続へ渡す。

$$
\boldsymbol a_\alpha=\boldsymbol r_{0,\alpha},\quad
\boldsymbol n=\frac{\boldsymbol a_1\times\boldsymbol a_2}{\|\boldsymbol a_1\times\boldsymbol a_2\|},\quad
z=\frac{H(\xi)\zeta}{2},\quad
\boldsymbol X_0=\boldsymbol r_0+z\boldsymbol n,\qquad -1\le\zeta\le1.
$$

$$
\boldsymbol X_{0,\alpha}=\boldsymbol a_\alpha+z\boldsymbol n_{,\alpha}
 +\frac{\zeta H_{,\alpha}}2\boldsymbol n,
\qquad \boldsymbol X_{0,\zeta}=\frac H2\boldsymbol n.
$$

Hは実法線方向の全厚。参照ヤコビアンの向き、最小特異値、正の厚さを求積点と指定の幾何検査点で検査し、反転・退化・評価不能を理由付きで拒否する。点検だけで全領域の非退化を数学的に証明したとは扱わない。

## 参照Ritz持上げと二階微分

荷重・現在状態・材料係数に依存しない曲率種から位置場の種を作る。接線補正は、参照の膜ひずみの物理テンソルノルムを最小化する固定線形問題とする。固定端など、種に定義した値・一階微分の条件を同じ問題に含める。剛体成分を除き、尺度付き平方根作用素のrankを複数閾値で照合する。準備後の係数は現在状態で更新しない。

持上げたベクトル場vの厚さ勾配dは、参照だけの線形方程式で定義する。

$$
L=\begin{bmatrix}\boldsymbol a_1^T\\\boldsymbol a_2^T\\\boldsymbol n^T\end{bmatrix},\qquad
L\boldsymbol d=
\begin{bmatrix}-\boldsymbol n\cdot\boldsymbol v_{,1}\\-\boldsymbol n\cdot\boldsymbol v_{,2}\\0\end{bmatrix}
=\boldsymbol b.
$$

$$
L\boldsymbol d_{,\alpha}=\boldsymbol b_{,\alpha}-L_{,\alpha}\boldsymbol d,
\qquad
b_{\beta,\alpha}=-\boldsymbol n_{,\alpha}\cdot\boldsymbol v_{,\beta}
 -\boldsymbol n\cdot\boldsymbol v_{,\beta\alpha},\quad b_{3,\alpha}=0.
$$

この評価にはvの二階微分と参照面の二階微分が必要。現在配置を差分して代用しない。幾何の点別曲率を柱面式へ代入する方法にも置換しない。スカラー形状関数は各成分 $$ \phi_j=v_j+zd_j $$ を全3Cartesian方向に用い、回転で閉じた位置空間にする。独立の基本r/d/cと組み合わせ、同じSで位置、速度、ひずみ、質量、面仕事を評価する。

現行の実行証拠は二次グラフ面に限る。上記の一般C²入力契約を記したことは、非二次面での品質の証拠ではない。

## 全トレースと拘束

鉛直面ではr、d、cの全厚さ多項式、上下面では同じSからの位置トレースを用いる。接続の仮想仕事は両側の実面積・向き・厚さ対応を含めて定義する。

$$
\delta W_\Gamma=\int_\Gamma\boldsymbol\lambda\cdot
(N_+\delta q_+-N_-\delta q_-)\,dA,
\qquad f_{\Gamma,+}=B_+^T\lambda,\quad f_{\Gamma,-}=-B_-^T\lambda.
$$

固定端の零トレース空間は、解析的に零と分かる種と基本場の同次成分を先に保持する。係数フィットの微小残差に任意精度消去を適用し、その残差を独立の物理拘束へ格上げしてはならない。既知の同次空間の包含、残余の全係数拘束、独立点での全厚さ零トレースを別々に検査する。中立面だけの一致で全トレースの合格を代用しない。

法線の反転には厚さ座標と層順の対応を伴わせる。位置と厚さ座標の対応を先に確定し、力には同じ写像の転置を使う。多パッチの完全な平衡解を検査していない場合、面作用素の仕事照合だけから接続要素の品質を主張しない。

## 材料・質量・状態

$$
F=I+\nabla_{X_0}u,\quad E=\frac12(F^TF-I),\quad
U=\sum_\ell\int_{\Omega_{0,\ell}}W_\ell(E-E_{0,\ell})\,dV_0,
\qquad M=\sum_\ell\int_{\Omega_{0,\ell}}\rho_{0,\ell}N^TN\,dV_0.
$$

材料点には参照材料軸・層ID・片側評価・指定初期ひずみを渡す。層界面をまたいだ応力の平均で片側値を置換しない。内力はUの勾配、接線は材料項と幾何項を含む同じ勾配の微分。参照固定線形Sと固定支持写像ではMは一定である。非線形拘束や状態依存基底ではこの定数質量の契約を流用しない。

確定・試行・棄却・保存・再開は[既存契約](../q9-dynamics/nonlinear-step-restart-contract.md)に従う。基底、参照、支持、材料、初期ひずみ、荷重、方式履歴の同一性を照合する。厳密モード伝播は再開検査の参照であり、製品積分器や連成J1の採用を意味しない。

## 主要処理の境界

| 内部処理 | 入力 → 出力 | 不変条件・拒否 |
|---|---|---|
| `prepare_reference` | 参照の二階jet、Hと勾配、種、支持規約 → 不変な基底係数・rank診断・同一性 | 荷重・材料・現在qから独立。幾何、rank、条件数の不適格を拒否 |
| `evaluate_reference` | 準備済み基底、物質点 → N、物理勾配G、X₀、J、測度 | 位置場と微分の一致、実法線厚さ、片側jet |
| `restrict_trace` | 面対応、向き、層／厚さ対応 → 両側トレース、力の転置写像 | 全厚さ一致、仕事保存、既知の同次場を保持 |
| `material_point` | F、参照軸、層、E₀ → W、第二PK、整合接線 | 客観性、層ごとの回復、有限値 |
| `assemble_mechanics` | 同じN/G、材料点、外力 → U、f、K、M | 物理体積測度、材料＋幾何接線、外力を内力から生成しない |
| `linearize_modes` | 適格な平衡、同じ支持・M/K → 符号を保つ固有対 | 負の固有値を除去しない。6剛体は無応力自由体で別確認 |
| `try/commit/reject/load` | 不変な確定状態、方式履歴 → 新確定状態または理由付き拒否 | 旧世代・二重確定・破損・履歴不足を拒否。全検査前に状態を置換しない |

これは内部の処理境界であり公開APIではない。採用可否と未検証範囲は結果文書で別に判定する。
