# 材料軸・初期ひずみ・内力・一貫接線の契約案

版0.1、2026-10-29。ST-03/10、CV-03、MAT-01/02、OP-02、INV-03/17。Q9の比較基準となる三次元・適合ひずみ・参照体積の材料経路を具体化する。材料の実験妥当性、定量的な適用限界、平面応力縮約・等価特性入力、製品APIの最終採用は未完了。

## 1. 測度・材料軸・初期ひずみ

参照直交物理基底に対する材料基底の列を $\mathbf A_m=[\mathbf a_1,\mathbf a_2,\mathbf a_3]$ とし、直交・右手系を要求する。参照材料点と層IDに付随し、試行変形で更新しない。曲面上の繊維軸の指定・生成法は別途定義する。

$$
\mathbf F=\mathbf J\mathbf J_0^{-1},\quad
\mathbf E=\frac12(\mathbf F^{\mathsf T}\mathbf F-\mathbf I),\quad
\mathbf E^m=\mathbf A_m^{\mathsf T}\mathbf E\mathbf A_m,\quad
\mathbf E_{\mathrm{el}}^m=\mathbf E^m-\mathbf E_0^m.
$$

初期ひずみE0は参照材料基底のGreen型対称テンソル。正の軸成分は自然長を伸ばす向きで、参照形状を保持すると圧縮応力を生み得る。有限値・対称性・ $\mathbf I+2\mathbf E_0^m$ の正定値を検査する。幾何・初期変位・初期応力の独立入力と混同しない。再開までModelに保持し、重複適用しない。

## 2. 工学成分と三次元材料入力

$$
\mathbf e=\begin{bmatrix}E_{11}&E_{22}&E_{33}&2E_{23}&2E_{13}&2E_{12}\end{bmatrix}^{\mathsf T},\quad
\mathbf s=\begin{bmatrix}S_{11}&S_{22}&S_{33}&S_{23}&S_{13}&S_{12}\end{bmatrix}^{\mathsf T}.
$$

この順序とせん断係数により $\mathbf s^{\mathsf T}\delta\mathbf e=\mathbf S:\delta\mathbf E$ 。応力側のせん断を2倍しない。外部ソフトの成分順は境界で変換する。

直交異方性の基準入力はE1/E2/E3、G23/G13/G12、ν12/ν13/ν23。νijはi方向単軸荷重によるj方向横収縮の比。 $\nu_{ji}=\nu_{ij}E_j/E_i$ を満たす対称コンプライアンスを構成する。

$$
\mathbf H=\begin{bmatrix}
1/E_1&-\nu_{12}/E_1&-\nu_{13}/E_1&0&0&0\\
-\nu_{12}/E_1&1/E_2&-\nu_{23}/E_2&0&0&0\\
-\nu_{13}/E_1&-\nu_{23}/E_2&1/E_3&0&0&0\\
0&0&0&1/G_{23}&0&0\\
0&0&0&0&1/G_{13}&0\\
0&0&0&0&0&1/G_{12}
\end{bmatrix},\qquad\mathbf C=\mathbf H^{-1}.
$$

各弾性率は正、尺度化したHを正定値とする。異方材のPoisson比へ等方材の区間条件をそのまま適用しない。三次元入力の根拠は [NTNUの構成則解説](https://nilspv.folk.ntnu.no/TMM4175/hookes-law.html)、工学定数・相反性は [COMSOLの線形弾性理論](https://doc.comsol.com/6.3/doc/com.comsol.help.sme/sme_ug_theory.06.026.html) と照合する。COMSOLの既定成分順を本書の順序へ流用しない。

面内4定数だけから厚さ方向の定数を補完しない。3D経路では不足入力を拒否し、平面応力縮約や等価シェル特性の経路は、その式・回復量・初期ひずみの扱いを別に審査する。製品で必要な複合材・等価構造の要求を削除した判断ではない。

## 3. エネルギーと応力

$$
W_0=\frac12(\mathbf e^m-\mathbf e_0^m)^{\mathsf T}\mathbf C^m(\mathbf e^m-\mathbf e_0^m),\quad
\mathbf s^m=\mathbf C^m(\mathbf e^m-\mathbf e_0^m).
$$

W0は元の参照体積当たり[J/m³]。第二PK応力SmはEと仕事共役。共通参照基底へ戻したS、第一PK応力P、現在配置のCauchy応力は、

$$
\mathbf S=\mathbf A_m\mathbf S^m\mathbf A_m^{\mathsf T},\quad
\mathbf P=\mathbf F\mathbf S,\quad
\boldsymbol\sigma=\frac1{\det\mathbf F}\mathbf F\mathbf S\mathbf F^{\mathsf T}.
$$

試行配置detF>0を要求する。応力出力に測度・配置・基底・層・評価点を付ける。初期ひずみを差し引く経路は、Greenひずみと第二PKを関係付けるStVK型の比較モデル。[MOOSEのStVK材料](https://mooseframework.inl.gov/source/materials/lagrangian/ComputeStVenantKirchhoffStress.html) は等方限定であり、本書では応力・ひずみ測度の照合にだけ使う。参照材料軸で一定の異方性Cを置く本候補のエネルギーと微分は、本書の式として導出・検証する。同ソフトの異方性対応を根拠にしない。

本モデルの回転客観性は有限な弾性ひずみへの物理的な適用性を保証しない。小弾性ひずみを想定するが、許容ひずみの製品基準は材料データと参照比較で別に定める。

## 4. 要素内力と接線

固定した参照補間・材料軸・初期ひずみに対し、スカラー基底の物理勾配をg、 $\mathbf F_j=\partial\mathbf F/\partial q_j$ とする。

$$
\mathbf E_j=\operatorname{sym}(\mathbf F^{\mathsf T}\mathbf F_j),\qquad
\mathbf E_{jk}=\operatorname{sym}(\mathbf F_j^{\mathsf T}\mathbf F_k).
$$

Ejを材料基底・工学成分へ変換した列がB。局所内部自由度・平面応力の消去がない本経路では、

$$
\Pi=\sum_\ell\int W_{0\ell}\,dV_0,\quad
\mathbf f=\sum_\ell\int\mathbf B^{\mathsf T}\mathbf s^m\,dV_0,\quad
\mathbf K=\mathbf K_{\mathrm{mat}}+\mathbf K_{\mathrm{geo}}.
$$

$$
\mathbf K_{\mathrm{mat}}=\sum_\ell\int\mathbf B^{\mathsf T}\mathbf C^m\mathbf B\,dV_0,\qquad
(K_{\mathrm{geo}})_{ia,jb}=\delta_{ab}\sum_\ell\int\mathbf g_i^{\mathsf T}\mathbf S\mathbf g_j\,dV_0.
$$

初期ひずみはSを通じて幾何剛性にも入る。Cが同じでもKは同じとは限らない。材料Cの正定値と、圧縮初期応力を持つ全接線の正定値を混同しない。外力・追従荷重・拘束・連成の接線は別に追加する。

## 5. 回転・座標変換の不変条件

現在配置だけへ剛体回転Qを重ねるとF→QF。参照材料軸・E0を固定し、E・Sm・W0は不変、P→QP、σ→QσQᵀ、要素内力は空間成分ごとにQで変換する。非ゼロ初期応力の内力を0とは要求しない。

非ゼロ内力を持つ状態で、回転経路の方向 $\mathbf v=\boldsymbol\Omega\mathbf q$ と二階方向 $\mathbf w=\boldsymbol\Omega^2\mathbf q$ について、客観性から

$$
\mathbf f^{\mathsf T}\mathbf v=0,\qquad
\mathbf v^{\mathsf T}\mathbf K\mathbf v+\mathbf f^{\mathsf T}\mathbf w=0.
$$

従って任意の予応力状態でKに6個の剛体零モードを要求しない。無応力平衡のランク試験と区別する。参照・現在の座標系を共に変える場合は、参照材料軸もQで変換する。

厚さ座標の再表現では、物理材料点・層IDと材料軸・E0を保持し、係数Tに対して内力はTᵀ、接線は合同変換で引戻す。局所配列の順序だけで材料を割り当てない。

## 6. 自然配置型との区別

 $\mathbf F_{\mathrm{el}}=\mathbf F\mathbf F_0^{-1}$ を用いるモデルとは、体積測度・弾性定数の基準・材料軸も合わせて比較する。例えば等方・一様な自然膨張 $\mathbf F_0=\lambda_0\mathbf I$ 、 $\lambda_0=\sqrt{1+2e_0}$ で、自然体積当たりに同じCを使うなら、元参照体積当たりのエネルギーは

$$
W_{\mathrm{mult},0}=\lambda_0^3\frac12\left(\frac{\mathbf E-e_0\mathbf I}{\lambda_0^2}\right):\mathbb C:\left(\frac{\mathbf E-e_0\mathbf I}{\lambda_0^2}\right)=\frac{W_{\mathrm{add},0}}{\lambda_0}.
$$

この比較は明示した材料定数・体積の定義での代数例で、材料精度の誤差評価ではない。任意の有限初期ひずみで両モデルが等価とはしない。初期Greenひずみ加算型を今回の比較基準とし、製品の適用限界や自然配置型の採否は別の審査に残す。

## 7. 処理I/F・失敗保証

| 項目 | 契約案 |
|---|---|
| 入力 | 参照幾何世代、層・材料ID、参照材料軸、E0、3D定数または定義済みC、試行q、評価位置と求積 |
| 材料点出力 | E・Eelの測度と基底、W0、Sm/S/P/σ、材料接線、detF・適格性診断 |
| 要素出力 | 同じ状態・同じ求積のΠ・f・Kmat・Kgeo・必要な回復値。別経路の応力で出力を上書きしない |
| 微分 | 参照幾何・層・軸・E0・求積を固定したq微分。内部変数等を追加する場合は消去も含む別接線 |
| 可変性 | Modelと確定状態は変更しない。材料点評価は直前の内力呼出しへ依存しない |
| キャッシュ | 参照勾配・軸・C・E0・層・求積・モデル世代をキーとする。試行依存のF・S・Bを参照キャッシュへ混入しない |
| 失敗 | 定数不足、非有限、材料非正定値、非直交／左手系材料軸、自然計量不適格、detF非正、未対応縮約を区別 |
| 寿命・性能 | 配列所有権・疎組立・スレッド／Python境界・製品依存・性能は未確定。実験スクリプトを製品APIにしない |

対応する [材料・接線試験](../../../rd/structure/main-shell/q9-dynamics/records/q9-material-tangent-experiment.md) で、解析解・テンソル形式・複素ステップ・差分刻み列・回転・厚さ座標変換を照合する。これらは材料の物理的妥当性、平衡・固有値・過渡・再開の受入とは区別する。

## 8. 技術判断と採用の現在地

MAT-D01：参照材料軸・Green初期ひずみ・元参照体積・工学成分順を比較基準として固定。MAT-D02：一貫接線に非ゼロSの幾何剛性を含め、材料の正定値と構造接線の性質を分離。MAT-D03：同じ材料点評価をエネルギー・内力・接線・応力回復へ使い、参照場は試行変形で更新しない。

上記の式・規約は、材料点、一軸例、層断面、RECT/NORMAL-Cの要素微分・変換試験で限定照合した。3D適合ひずみ経路の基準であり、任意の有限初期ひずみ・物理的な材料範囲・平面応力縮約・材料入力の全経路の採用完了ではない。製品の型・依存・性能、初期化・安定性・再開を閉じるまで本体実装の開始条件は未達。
