# EAS21残差の構成則・厚さ・平均による分離

2026-11-02。ST-04の原因診断。保存共有解を変更しない。[事前仕様](records/q9-eas21-residual-experiment.md)、[結果](records/q9-eas21-residual-results.md)。以下は本プロジェクトの診断用分解であり、製品の膜／曲げ出力定義を確定しない。

## 構成則の結合を残す分離

誤差ひずみをeとし、中立面の物理直交基底でMandel成分を並べ替える。

$$
p=(e_{11},e_{22},\sqrt2e_{12})^T,\quad
n=e_{33},\quad s=(\sqrt2e_{23},\sqrt2e_{13})^T.
$$

等方材のCについて、面内ブロックA、面内－厚さ結合b、厚さ成分cを用いると

$$
e^TCe=p^TC_{\rm ps}p+c\left(n+c^{-1}b^Tp\right)^2
+2\mu s^Ts,\qquad C_{\rm ps}=A-bc^{-1}b^T.
$$

$$
b=\lambda(1,1,0)^T,\quad c=\lambda+2\mu,\quad
C_{\rm ps}=2\mu I_3+
\frac{2\mu\lambda}{\lambda+2\mu}(1,1,0)^T(1,1,0).
$$

3項は非負で総和が全3D誤差エネルギーとなる。ここで平面応力Schur補行列を使うのは平方完成のためであり、元の全3D構成則を平面応力へ変更しない。厚さ項には面内とのPoisson結合を残す。

## 実体積重みによる厚さ直交分離

各面内点で実体積密度j(z)>0を用い、次を定義する。

$$
J_0=\int j\,dz,\quad
\bar z=J_0^{-1}\int zj\,dz,\quad
J_2=\int(z-\bar z)^2j\,dz,
$$

$$
p_0=J_0^{-1}\int pj\,dz,\quad
p_1=J_2^{-1}\int(z-\bar z)pj\,dz,\quad
p_\perp=p-p_0-(z-\bar z)p_1.
$$

$$
\int p^TC_{\rm ps}pj\,dz
=J_0p_0^TC_{\rm ps}p_0+J_2p_1^TC_{\rm ps}p_1
+\int p_\perp^TC_{\rm ps}p_\perp j\,dz.
$$

これらを膜型・曲げ型・高次厚さの誤差と呼ぶ。中立面値とその微分で定義する漸近膜・曲げとの同一視はしない。横せん断・結合を含む厚さ誤差と合わせて5項となる。離散積分でも同じ正の重みを使うため、分離の等式が成立する。

膜型p₀については、各要素の面内実体積重みで平均と変動へさらに分ける。この平均は位置依存する局所直交基底の成分平均であり、以下の物理固定軸6成分の平均とは区別する。

## 要素平均による別の分離

均質C、物理固定軸Mandel成分で

$$
\bar e_e=V_e^{-1}\int e\,dV,\qquad
\int e^TCe\,dV
=V_e\bar e_e^TC\bar e_e
+\int(e-\bar e_e)^TC(e-\bar e_e)\,dV.
$$

これは5項の分離とは別の軸であり、平均項を5項へ重複加算しない。EAS21の補正は平均ゼロなので、同じ物理変位に対してこの平均は生Bと消去後Bで一致する。現在の平均エネルギーは、その変位を固定して平均ゼロ補正だけを変える場合の下限であり、共有変位空間全体の下限ではない。

正規化はすべて同じ解析参照の全3D二次形式で行う。

$$
\rho_i=\frac{\int q_i(e)\,dV}{\int (\varepsilon^*)^TC\varepsilon^*\,dV},
\qquad \eta=\sqrt{\sum_{i=1}^5\rho_i}.
$$

エネルギー比の百分率と、エネルギーノルム誤差の百分率を混同しない。
