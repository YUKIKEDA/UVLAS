# 3833型シェルの微分・材料・内力・接線・求積仕様案

版：0.1、2026-10-29。ST-01〜05、ST-10の設計。比較対象は [固定版3833実装](https://github.com/projectchrono/chrono/blob/d26ed26d8344ae9aee4c1d26bf6d481076750b1a/src/chrono/fea/ChElementShellANCF_3833.cpp)。ソースの読解とUVLAS側の導出を区別し、最終採用・実装済みとはしない。

## 1. 3833-D01：保持する運動学

[候補比較](curved-shell-comparison.md) の8節点補間を使う。節点ごとに位置 $\mathbf r_i$ [m]、厚さ一次微分 $\mathbf d_i$ [無次元]、厚さ二次微分 $\mathbf c_i$ [ $\mathrm{m}^{-1}$ ]、計72座標を持つ。

$$
\mathbf b_i(z)=\mathbf r_i+z\mathbf d_i+\frac{z^2}{2}\mathbf c_i,
\qquad
\mathbf x_h=\sum_{i=1}^8N_i\mathbf b_i=\mathbf S\mathbf q_e.
$$

$$
\mathbf S_i=[N_i\mathbf I_3,\ zN_i\mathbf I_3,\ \tfrac12z^2N_i\mathbf I_3].
$$

層の物理厚さ範囲を $[z_\ell^-,z_\ell^+]$ 、層内親座標を $\zeta_\ell\in[-1,1]$ とする。

$$
h_\ell=z_\ell^+-z_\ell^->0,\qquad
z_{c\ell}=\frac{z_\ell^-+z_\ell^+}{2},\qquad
z=z_{c\ell}+\frac{h_\ell}{2}\zeta_\ell.
$$

全層で同じ基準面 $z=0$ を使う。層内座標の原点を、節点の位置・微分座標の基準面へ読み替えない。厚さ方向オフセットは下端・中央・基準面のどれかを明記し、曖昧な一つのoffset引数へまとめない。

## 2. 3833-D02：面内形状関数の必要微分

隅点の符号を $a_i=\xi_i\in\{-1,1\}$ 、 $b_i=\eta_i\in\{-1,1\}$ と置く。

$$
N_i=\frac14(1+a_i\xi)(1+b_i\eta)(a_i\xi+b_i\eta-1).
$$

$$
\begin{aligned}
N_{i,\xi}&=\frac{a_i}{4}(1+b_i\eta)(2a_i\xi+b_i\eta),\\
N_{i,\eta}&=\frac{b_i}{4}(1+a_i\xi)(a_i\xi+2b_i\eta),\\
N_{i,\xi\xi}&=\frac12(1+b_i\eta),\qquad
N_{i,\eta\eta}=\frac12(1+a_i\xi),\\
N_{i,\xi\eta}&=\frac{a_i b_i}{4}(2a_i\xi+2b_i\eta+1).
\end{aligned}
$$

辺中央の順はE=(0,−1)、F=(1,0)、G=(0,1)、H=(−1,0)。形状関数自体は候補比較のCS-02と同じ。

$$
\begin{array}{c|ccccc}
i&N_{i,\xi}&N_{i,\eta}&N_{i,\xi\xi}&N_{i,\xi\eta}&N_{i,\eta\eta}\\\hline
E&-\xi(1-\eta)&-\tfrac12(1-\xi^2)&-(1-\eta)&\xi&0\\
F&\tfrac12(1-\eta^2)&-\eta(1+\xi)&0&-\eta&-(1+\xi)\\
G&-\xi(1+\eta)&\tfrac12(1-\xi^2)&-(1+\eta)&-\xi&0\\
H&-\tfrac12(1-\eta^2)&-\eta(1-\xi)&0&\eta&-(1-\xi)
\end{array}
$$

## 3. 3833-D03：厚さ微分と曲面参照写像

$$
\begin{aligned}
\mathbf x_{,\xi}&=\sum_iN_{i,\xi}\mathbf b_i,&
\mathbf x_{,\eta}&=\sum_iN_{i,\eta}\mathbf b_i,\\
\mathbf x_{,z}&=\sum_iN_i(\mathbf d_i+z\mathbf c_i),&
\mathbf x_{,zz}&=\sum_iN_i\mathbf c_i,\\
\mathbf x_{,\xi z}&=\sum_iN_{i,\xi}(\mathbf d_i+z\mathbf c_i),&
\mathbf x_{,\xi\eta}&=\sum_iN_{i,\xi\eta}\mathbf b_i.
\end{aligned}
$$

親座標 $\mathbf s=(\xi,\eta,\zeta_\ell)$ を用いる場合、厚さ微分には $h_\ell/2$ を掛ける。参照配置にも同じ式・層写像を使い、

$$
\mathbf J_0=\frac{\partial\mathbf X_h}{\partial\mathbf s},\qquad
\mathbf J=\frac{\partial\mathbf x_h}{\partial\mathbf s},\qquad
\mathbf F=\mathbf J\mathbf J_0^{-1},\qquad
\mathrm dV_0=\det\mathbf J_0\,\mathrm d\xi\,\mathrm d\eta\,\mathrm d\zeta_\ell.
$$

層ごとの厚さ係数は $\mathbf J_0$ に含む。積分重みへ再度 $h_\ell/2$ を掛けない。逆行列の数学記号を、明示逆行列を必ず生成する実装指示と扱わない。向き・条件数・非退化の診断を行う。

ここではANS/EASによるひずみ修正を加えない適合ひずみ経路を基準候補とする。厚さ二次座標があるだけで薄肉ロッキングが解消したとは判定しない。

## 4. 3833-D04：材料軸と指定初期ひずみ

材料の参照直交基底を $\mathbf A_m=[\mathbf a_1,\mathbf a_2,\mathbf a_3]$ とする。部材基底と繊維角から構成し、曲面では評価点に対応する基底を用いる。現在の試行変形に合わせて参照材料軸を再設定しない。

$$
\mathbf E=\frac12(\mathbf F^{\mathsf T}\mathbf F-\mathbf I),\qquad
\mathbf E^m=\mathbf A_m^{\mathsf T}\mathbf E\mathbf A_m,
\qquad
\mathbf E_{\mathrm{el}}^m=\mathbf E^m-\mathbf E_0^m.
$$

小弾性ひずみ・大回転の材料候補は

$$
W_\ell=\frac12\mathbf E_{\mathrm{el}}^m:\mathbb C_\ell^m:\mathbf E_{\mathrm{el}}^m,
\qquad
\mathbf S^m=\mathbb C_\ell^m:\mathbf E_{\mathrm{el}}^m.
$$

これはQ27 Bの入力範囲を具体化する材料候補で、任意の有限初期ひずみや任意の残留応力場との等価性を保証しない。自然配置型との比較はST-03に残る。

体積直交異方弾性を使うなら、3方向のYoung率・3方向のせん断弾性率・独立なPoisson比と相反関係・正定値性が必要。面内の薄板材料定数だけから厚さ方向定数を黙って補完しない。平面応力へ縮約するなら別の材料・要素仕様として審査する。

## 5. 3833-D05：内力・接線をエネルギーから導く

節点座標成分 $q_j$ に対して

$$
\mathbf F_j=\frac{\partial\mathbf F}{\partial q_j}
=\frac{\partial\mathbf J}{\partial q_j}\mathbf J_0^{-1}
$$

はqに依存しない。 $\operatorname{sym}\mathbf H=(\mathbf H+\mathbf H^{\mathsf T})/2$ と置くと、

$$
\mathbf E_j=\operatorname{sym}(\mathbf F^{\mathsf T}\mathbf F_j),\qquad
\mathbf E_{jk}=\operatorname{sym}(\mathbf F_j^{\mathsf T}\mathbf F_k).
$$

材料系へ同じ $\mathbf A_m$ で変換し、

$$
\Pi_e=\sum_\ell\int_{V_{0\ell}}W_\ell\,\mathrm dV_0,
\qquad
f_j=\sum_\ell\int_{V_{0\ell}}\mathbf E_j^m:\mathbf S^m\,\mathrm dV_0,
$$

$$
K_{jk}=\sum_\ell\int_{V_{0\ell}}
\left(\mathbf E_j^m:\mathbb C_\ell^m:\mathbf E_k^m
+\mathbf S^m:\mathbf E_{jk}^m\right)\,\mathrm dV_0.
$$

初期ひずみに由来する応力も第二項の幾何剛性へ含む。外力・減衰・拘束・写像の接線はこの弾性内力接線と区別する。

## 6. 3833-D06：質量・単位・尺度

$$
\mathbf M_e=\sum_\ell\int_{V_{0\ell}}\rho_{0\ell}\mathbf S^{\mathsf T}\mathbf S\,\mathrm dV_0,
\qquad
\mathbf Q_g=\sum_\ell\int_{V_{0\ell}}\rho_{0\ell}\mathbf S^{\mathsf T}\mathbf g\,\mathrm dV_0.
$$

参照補間・密度・層厚さが固定ならMは定数で、慣性項は $\mathbf M_e\mathbf a_e$ 。72座標を回転角へ置き換えた式ではない。

| ブロック | 単位 |
|---|---|
| $\mathbf M_{rr}$ | kg |
| $\mathbf M_{rd}$ | kg·m |
| $\mathbf M_{rc}$ 、 $\mathbf M_{dd}$ | kg·m² |
| $\mathbf M_{dc}$ | kg·m³ |
| $\mathbf M_{cc}$ | kg·m⁴ |

一般化荷重は位置・厚さ一次・厚さ二次の順でN、N·m、N·m²。混合ベクトルへ単位のない一つの絶対ノルムを適用しない。

長さ尺度 $L_*$ に対して

$$
\mathbf q_e=\mathbf D_q\widehat{\mathbf q}_e,\qquad
\mathbf D_{q,i}=\operatorname{diag}(L_*\mathbf I_3,\mathbf I_3,L_*^{-1}\mathbf I_3)
$$

で無次元座標を作り、接線は $\widehat{\mathbf K}=\mathbf D_q^{\mathsf T}\mathbf K\mathbf D_q$ 、質量も同じ合同変換を用いる。大きさを揃えずにraw固有値の絶対閾値でランクを判定しない。

## 7. 3833-D07：質量と弾性の求積を分離する

固定版の `ComputeMassMatrixAndGravityForce` は各層で5×5×5点、内力の既定は3×3×2点。ヘッダーのNP/NTを質量へ流用しない。

**質量基準：** 各層一定密度、上記の多項式参照補間、適格な向きの場合、 $\mathbf S^{\mathsf T}\mathbf S$ は各親座標で最大4次、 $\det\mathbf J_0$ は最大5次という上界を持つ。積は最大9次なので、各方向5点Gaussで多項式質量を厳密積分できる。任意の密度場・非多項式幾何へこの結論を拡張しない。

**弾性基準：** 参照Jacobianと材料軸が一定、Cと初期ひずみが層ごとに一定なら、適合StVKエネルギーの各親座標の次数上界は8次で、5点Gaussを比較の基準経路にできる。曲面では参照逆Jacobianや材料基底が位置依存するため、5点を厳密解と呼ばない。

3×3×2点の弾性経路は低減積分候補として別の離散エネルギーを定義する。同じ点・重みで内力・接線を導き、完全積分経路との一致・差・薄肉感度・零モードを評価してから採否を決める。質量の低減積分は別の判断であり、厚さ2点をそのまま採らない。

曲面の求積比較は5、7、9点等で行い、出力量・接線・固有特性への影響を記録する。最終2条件の差だけで誤差上限と断定せず、非漸近・別の支配誤差も点検する。

## 8. 3833-D08：応力回復・初期化・処理境界

応力出力は内力と同じ材料式・初期ひずみ・層・評価点から得る。第二PK応力を共通参照系へ戻し、

$$
\mathbf S_I=\mathbf A_m\mathbf S^m\mathbf A_m^{\mathsf T},\qquad
\boldsymbol\sigma=\frac1{\det\mathbf F}\mathbf F\mathbf S_I\mathbf F^{\mathsf T}.
$$

現在配置に剛体回転だけを重ねた場合、参照の第二PK成分は変わらず、Cauchy応力と空間力が回転する。参照・現在の両座標系を変える共変性試験とは区別する。

要素評価は、不変Model、試行q・v、Workspaceを分ける。出力はエネルギー、72成分内力、72×72接線、同じ状態の回復値と診断。初期化時のM・幾何・変換・求積のキャッシュはモデル世代を検査する。内力を直前に呼んだことへ暗黙依存せず、確定Stateを書き換えない。

参照形状が既に曲面でも、厚さ二次座標の初期値を必ずゼロと仮定しない。所望の参照体積をどの初期微分で表すかをモデル作成・外部格子変換で定義する。

## 9. 採用の現在地

補間・必要微分・参照写像・材料候補・内力・接線・質量・求積の役割を式で定義した。解析的な採否ケースは [3833参照ケース](verification/3833-reference-cases.md) に置く。

未完了は、曲面材料基底・入力変換の具体契約、低減積分の安定性、薄肉・歪み・積層初期ひずみ・独立FEM比較、シェル交差／梁接続と性能。従って、3833を検証済みの主要素とする判定と本体実装開始はまだ行わない。

[自由単要素の代数試験](verification/3833-algebraic-experiment.md) では、5×5×5と3×3×2の弾性接線はrank66、厚さ2点の質量はrank48と確認した。今回の平坦・無応力状態に限った証拠として、質量への厚さ2点流用と安定化なしの2×2×2弾性経路を除外。3×3×2弾性は次の比較対象に残す。

続く [薄肉・歪み試験](verification/3833-bending-distortion-experiment.md) では、台形で純曲げ場の節点補間に人工せん断が生じ、求積削減だけでは除けないことを確認した。今回の適合ひずみ・補正なし3833型を主シェルとしてそのまま採用せず、[9節点比較モデル](q9-thickness-quadratic-design.md) を優先する。上記のランク成立はこの歪み感度を否定する証拠ではない。
