# 局所圧力・後縁極限・剛体ピッチの独立参照

版：0.1、2026-10-30。AE-02/04とTR-02の比較用。ゼロ循環平板まわりの2D・線形・非圧縮・調和問題。

## 1. 参照の独立性と座標

[Floryan–Rowley著者公開稿](https://dfloryan.github.io/publications/floryan2018clarifying.pdf) 付録A、式(A8)〜(A13)の加速度ポテンシャル展開を剛体運動へ特殊化する。同稿の半弦座標を物理座標へ戻し、UVLASの上向き圧力を正とする。揚力・モーメントは [NASA/TM–2017–219667](https://ntrs.nasa.gov/api/citations/20170009599/downloads/20170009599.pdf) の上下変位・回転軸・正方向を変換した式とも照合する。以下の特殊化・観測規約はUVLAS側の導出。

半弦b=c/2、 $X=(x-b)/b$ 、回転軸 $a=(x_p-b)/b$ 、換算振動数 $k=\omega b/U$ 。上向きH、正の鼻上げピッチαに対し、変位の複素振幅は

$$
\widehat z(x)=H-\alpha(x-x_p),\qquad
\eta=H/b,\qquad C(k)=\frac{K_1(ik)}{K_0(ik)+K_1(ik)}.
$$

参照は修正Bessel K、既存の荷重参照はHankelを使い、算術も別表現で照合する。実空間の虚数と時間調和の虚数を混同しない。k=0ではC=1の極限を明示する。

## 2. 剛体運動に対する圧力

無次元変位 $Y_0=\eta-\alpha(X-a)$ に対して、 $D\Psi=-(ik+D)^2Y_0$ と定数・一次の面速度から

$$
B_0=-2ikC\eta+\left[2C\{1+ik(1/2-a)\}-ik\right]\alpha,
$$
$$
B_1=k^2\eta+(2ik+k^2a)\alpha,\qquad B_2=-k^2\alpha/4.
$$

従って、 $R=\sqrt{(1-X)/(1+X)}$ 、 $S=\sqrt{1-X^2}$ として

$$
\widehat p_z(x)=\rho U^2\left[B_0R+2B_1S+4B_2XS\right],\qquad 0<x<c.
$$

このB係数は渦輪Γでもポテンシャル差μでもない。UVLMの循環解を入力にせず、運動・流れ・解析循環関数だけから計算する。前縁の逆平方根特異性を有限点圧力として保存しない。後縁の値は解析極限0である。

## 3. 圧力・荷重・後縁平均の照合

独立な全揚力と軸周りの鼻上げモーメントは

$$
\widehat F_z'=\rho\pi U^2b(B_0+B_1),\qquad
\widehat M_y'(x_p)=\rho\pi U^2b^2\left[(a+1/2)B_0+aB_1-B_2/2\right].
$$

式の再代入だけで検査せず、圧力をθ座標でGauss積分した荷重と、展開前のTheodorsen荷重式とも比較する。圧力の区間平均は同じ区間・同じ単位で [弱圧力観測](uvlm-graded-wake-pressure-contract.md) と比較する。

後縁幅w→0の平均には

$$
\lim_{w\to0}\frac{1}{\sqrt{w/c}}\frac1w\int_{c-w}^{c}\widehat p_z\,dx
=\frac23\rho U^2(B_0+4B_1+8B_2)
$$

という独立な漸近係数がある。平均圧力が0へ減るだけでなく、この係数とも比較できる。有限格子のP0点圧力が連続になることを要求する式ではない。

## 4. ピッチの非貫通と荷重点

線形化した非貫通の右辺は

$$
\widehat u_i=i\omega H-\left[U+i\omega(x_i^c-x_p)\right]\alpha.
$$

法線の回転による−Uαを落とさない。基準循環0に限るため、係数の幾何摂動と循環の積、荷重点摂動と力の積は二次。非ゼロ迎角・キャンバ・予変形平衡へこの省略を一般化しない。

循環辺力は既存の辺中点、非循環面力は既存の三角重心を用い、 $M_y=-\sum(x-x_p)F_z$ を計算する。同じP0場の保存則からの別計算は

$$
\widehat M_y'=-\rho\left[d_\Delta\int_0^c(x-x_p)\widehat\mu_h\,dx
+U\left\{(c-x_p)\widehat\mu_{TE}-\int_0^c\widehat\mu_h\,dx\right\}\right].
$$

観測基準点だけをΔx動かす場合は $M_y(x_p+\Delta x)=M_y(x_p)+\Delta xF_z$ 。運動の回転軸を変える操作とは区別する。

## 5. 決定と残る範囲

- LPP-D01：解析圧力は独立な加速度ポテンシャル展開から作り、UVLMの循環解や力から分布を逆算しない。
- LPP-D02：P0圧力は同一区間の弱観測で比較し、後縁平均の平方根係数と局所分布を別に検査する。
- LPP-D03：ピッチの法線項、作用点・基準点、辺／面のモーメントと保存則を一組で検査する。

[事前仕様・結果](records/uvlm-local-pressure-pitch-experiment.md) に限定照合を記録する。実務有限幅・曲面・複数面・有限振幅・自由後流・ANCF強連成・J1は未完了。

## 求解精度とトレースの追加診断

[精度契約](uvlm-residual-and-trace-accuracy.md) と [補償残差・圧力トレース試験](records/uvlm-residual-trace-experiment.md) で、独立圧力からμを積分し、局所平均の差を区間積分と左右トレースへ分解した。反復改良で残差は成立したが圧力の格子差は残り、単一の残差値や総荷重だけでは局所圧力の適格性を判断しない。

