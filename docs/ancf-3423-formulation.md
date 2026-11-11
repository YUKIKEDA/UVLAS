# 3423型シェル：補間・ANS/EAS・接線の数式仕様案

版：0.1、2026-10-29。ST-01〜05、ST-10の候補定式化。採用・実装開始判定は未完了。式を定義した範囲と、参照実装の挙動、独立検証の不足を区別する。

追記：[層別EAS整合審査](eas-layer-consistency.md) により、層別αへ要素全体ζを使う案と、SH-04を曲面へそのまま適用する案は不採用とした。以下の参照形は比較・反例の対象として残す。現在のEAS設計条件はEAS-D02〜04、曲面ANSの採用は保留。

## 1. 出典と照合範囲

[Chrono技術報告（2016）](https://www.projectchrono.org/assets/validations/FEA/ancfShellValidation.pdf) §2–3、式(1)–(11)を補間・材料・質量の出発点とした。ANS/EASの具体配置・補間・局所消去は、公式ソースの固定revision `d26ed26d8344ae9aee4c1d26bf6d481076750b1a` で照合した。

- [3423実装](https://github.com/projectchrono/chrono/blob/d26ed26d8344ae9aee4c1d26bf6d481076750b1a/src/chrono/fea/ChElementShellANCF_3423.cpp)：ShapeFunctions、CalcStrainANSbilinearShell、Basis_M、ComputeInternalForces、ComputeInternalJacobians。
- [3423ヘッダー](https://github.com/projectchrono/chrono/blob/d26ed26d8344ae9aee4c1d26bf6d481076750b1a/src/chrono/fea/ChElementShellANCF_3423.h)：節点型・順序・内部変数。
- 取得・ハッシュ・照合上の懸案は [固定版ソース審査](verification/chrono-shell-source-audit.md) に記録。

以下の仮想仕事・縮約・初期ひずみの式はUVLAS側で整理した設計案。参照実装と完全一致する要素を採用した、あるいは初期曲面・初期ひずみの全要求を文献が保証した、という意味ではない。

## 2. SH-01：局所座標、節点、自由度

親座標は $\xi,\eta,\zeta\in[-1,1]$ 。節点順はA(−1,−1)、B(+1,−1)、C(+1,+1)、D(−1,+1)。要素厚さ $h>0$ 、層境界 $\zeta_{\ell}$ と材料・密度・繊維角を持つ。

節点iの未知量は絶対位置 $\mathbf r_{i}$ [m]と横方向勾配 $\mathbf d_{i}$ [無次元]。要素座標は

$$
\mathbf q_e=[\mathbf r_A^{\mathsf T},\mathbf d_A^{\mathsf T},\mathbf r_B^{\mathsf T},\mathbf d_B^{\mathsf T},\mathbf r_C^{\mathsf T},\mathbf d_C^{\mathsf T},\mathbf r_D^{\mathsf T},\mathbf d_D^{\mathsf T}]^{\mathsf T}\in\mathbb R^{24}
$$

の順に各3成分を並べる。 $\mathbf d_{i}$ は角度や単位ベクトルではなく、厚さ方向の伸縮・せん断も表現する。現在の $\mathbf d_{i}$ を毎回単位化すると別の運動学になる。

参照値 $\mathbf q_{e0}$ は初期の位置・勾配を同じ順で保持する。初期 $\mathbf d_{i}$ は参照面・厚さの定義に整合させる。交差部材で位置が一致しても、異なる方向の $\mathbf d_{i}$ をそのまま共有しない。

## 3. SH-02：補間と参照写像

節点の親座標を( $\xi_{i}$ , $\eta_{i}$ )とし、

$$
N_i(\xi,\eta)=\frac14(1+\xi_i\xi)(1+\eta_i\eta)
$$

$$
z=\frac h2\zeta
$$

$$
\mathbf x_h=\sum_iN_i(\mathbf r_i+z\mathbf d_i)=\mathbf S\mathbf q_e
$$

$$
\mathbf X_h=\sum_iN_i(\mathbf r_{i0}+z\mathbf d_{i0})=\mathbf S\mathbf q_{e0}
$$

とする。Sは3×24、節点ごとのブロックは $[N_{i} I_{3}, z N_{i} I_{3}]$ 。従って位置と勾配への微分はそれぞれ $N_{i}$ とz $N_{i}$ になる。

物理単位を持つ親座標を $p=(a\xi /2,b\eta /2,h\zeta /2)$ と置く。a,b>0は親座標の尺度で、歪んだ実要素を長方形と仮定するものではない。

$$
\frac{\partial N_i}{\partial p_1}=\frac2a\frac{\partial N_i}{\partial\xi}\qquad \frac{\partial N_i}{\partial p_2}=\frac2b\frac{\partial N_i}{\partial\eta}\qquad \frac{\partial(zN_i)}{\partial p_3}=N_i
$$

を用い、 $A=\partial x_{h}/\partial p$ 、 $A_{0}=\partial X_{h}/\partial p$ 、 $F=A A_{0}^{-1}$ とする。 $\mathbf A_{0}$ を単位行列に置き換えない。初期曲面もこの参照写像で表す。

参照体積要素は

$$
\mathrm dV_0=\det(\mathbf A_0)\frac{abh}{8}\,\mathrm d\xi\,\mathrm d\eta\,\mathrm d\zeta
$$

となる。適格な向き・非退化を確認し、負のdetを絶対値で隠さない。求積点だけで要素内部全域の適格性が証明されたとはしない。

## 4. SH-03：適合ひずみと材料基底

親座標基底の共変Greenひずみを

$$
\widetilde{\mathbf E}=\frac12(\mathbf A^{\mathsf T}\mathbf A-\mathbf A_0^{\mathsf T}\mathbf A_0)
$$

とする。本書の中間成分順は参照実装との照合用に

$$
\widetilde{\mathbf e}=[\widetilde E_{11},\widetilde E_{22},2\widetilde E_{12},\widetilde E_{33},2\widetilde E_{13},2\widetilde E_{23}]^{\mathsf T}
$$

を使う。既存の外部出力順`11,22,33,23,13,12`とは境界で変換し、同じ配列と扱わない。

参照材料基底の列を $\mathbf a_{m}$ とし、 $\beta _{m}=A_{0}^{-1} a_{m}$ なら

$$
E^m_{mn}=\boldsymbol\beta_m^{\mathsf T}\widetilde{\mathbf E}\boldsymbol\beta_n
$$

で材料軸へ変換する。この関係から工学せん断を含む6×6変換Tを構成する。参照基底・繊維角・材料はqの試行に依存しない。厚さを含む位置依存の $\mathbf A_{0}$ と材料軸を評価する。

## 5. SH-04：ANSの評価点と補間

適合ひずみのうち、 $\widetilde{e}_{11}$ , $\widetilde{e}_{22}$ ,2e~12は評価点で直接計算する。厚さ垂直・横せん断は次のサンプリング案に置き換える。全点でζ=0。

| 成分 | サンプリング点 | 補間 |
|---|---|---|
| $\widetilde{E}_{33}$ | A,B,C,Dの4隅 | $\sum _{i} N_{i}(\xi ,\eta ) \widetilde E_{33}(\xi _{i},\eta _{i},0)$ |
| 2E~13 | (0,−1,0)、(0,+1,0) | ηに対する線形補間 |
| 2E~23 | (−1,0,0)、(+1,0,0) | ξに対する線形補間 |

これを $\widetilde{e}_{\mathrm{ANS}}(q)$ とする。サンプル値と補間重みは座標で対応付け、配列の添字順が一致していると仮定しない。サンプル点のB・二階微分にも同じ補間を適用する。

曲面では共変成分の補間と材料軸変換の順序が重要。 $T(\xi ,\eta ,\zeta )\widetilde{e}_{\mathrm{ANS}}$ と、各点の材料ひずみを先に補間する操作を同一とみなさない。薄肉・歪み・曲面・剛体運動のケースでこの案を審査する。

参照実装の隅点列と重みには対応上の懸案があるため、本書は上表の座標対応を仕様候補とする。[審査記録](verification/chrono-shell-source-audit.md) のSRC-V01に反例を記載。

## 6. SH-05：EASと指定初期ひずみ

参照実装では層ごとに5個の内部パラメータ $\alpha_{\ell}$ を持ち、親座標基底でのモード行列は

$$
\mathbf P(s)=\begin{bmatrix}
\xi&0&0&0&0\\
0&\eta&0&0&0\\
0&0&\xi&\eta&0\\
0&0&0&0&\zeta\\
0&0&0&0&0\\
0&0&0&0&0
\end{bmatrix}
$$

という形で、中心の変換 $\mathbf T_{c}$ 、参照det比を用いて $G(s)=T_{c} P(s) \det(A_{0c})/\det(A_{0}(s))$ としている。この参照形のままの層別モードはEAS-V01の反例により不採用。以降のGは、層内厚さ座標 $t_{\ell}$ 、共通参照系での零平均と材料系への点別変換を使う [EAS-C02/03](eas-layer-consistency.md) の候補へ置き換えて比較する。

UVLASの初期ひずみを含むエネルギー候補は、材料軸で

$$
\mathbf e_{\mathrm{el}}=\mathbf T(s)\widetilde{\mathbf e}_{\mathrm{ANS}}(\mathbf q)+\mathbf G(s)\boldsymbol\alpha_\ell-\mathbf e_0^m(s)
$$

$$
W_\ell=\frac12\mathbf e_{\mathrm{el}}^{\mathsf T}\mathbf C_\ell^m\mathbf e_{\mathrm{el}}
$$

$$
\Pi_e(\mathbf q,\boldsymbol\alpha)=\sum_\ell\int_{V_{0\ell}}W_\ell\,\mathrm dV_0
$$

とする。Cのせん断係数は工学せん断に対応させる。初期ひずみは指定した材料場として積分点で評価し、ANSやEASで無断に平均・除去しない。特に層差を失わないことを検証する。

ここでCは体積弾性の候補。平面応力モデルと同じ式とは扱わず、厚さ応力を緩和する設計と入力可能な初期ひずみの関係を比較する。

## 7. SH-06：局所消去、内力、一貫接線

$\mathbf B=\frac{\partial\mathbf e_{\mathrm{el}}}{\partial\mathbf q}$ （αを固定）、 $s=C e_{\mathrm{el}}$ として

$$
\mathbf f_q=\frac{\partial\Pi_e}{\partial\mathbf q}=\sum_\ell\int_{V_{0\ell}}\mathbf B^{\mathsf T}\mathbf s\,\mathrm dV_0
$$

$$
\mathbf r_\alpha=\frac{\partial\Pi_e}{\partial\boldsymbol\alpha}=\int\mathbf G^{\mathsf T}\mathbf s\,\mathrm dV_0=\mathbf0
$$

を満たすαを求める。固定したqと線形材料、q非依存Gでは、αに関する式は線形である。一般の非線形材料へ拡張する場合に局所非線形反復が必要になる。参考実装がNewtonループを使うことを、UVLASにも同じループが必要な根拠にはしない。

$$
\mathbf K_{qq}=\sum_\ell\int_{V_{0\ell}}\left[\mathbf B^{\mathsf T}\mathbf C\mathbf B+\sum_js_j\frac{\partial^2e_{\mathrm{el},j}}{\partial\mathbf q^2}\right]\,\mathrm dV_0
$$

$$
\mathbf K_{q\alpha}=\int\mathbf B^{\mathsf T}\mathbf C\mathbf G\,\mathrm dV_0\qquad \mathbf K_{\alpha\alpha}=\int\mathbf G^{\mathsf T}\mathbf C\mathbf G\,\mathrm dV_0
$$

局所平衡を満たす縮約接線は

$$
\mathbf K_{\mathrm{red}}=\mathbf K_{qq}-\mathbf K_{q\alpha}\operatorname{solve}(\mathbf K_{\alpha\alpha},\mathbf K_{\alpha q})
$$

である。逆行列を明示生成する仕様にしない。 $\mathbf K_{\alpha \alpha}$ の適格性、尺度付き局所残差、非収束・特異・非有限値の失敗保証を定義する。局所平衡が未成立のまま上式だけを返さない。

残留局所残差があると大域Newtonの右辺にも対応する凝縮補正が必要となる。UVLASの最初の契約候補では、要素成功の条件として局所平衡を要求し、大域へ未収束の局所状態を正常値として返さない。

Q27 Bの初期ひずみはsを通じて幾何剛性へ影響する。安定性の接線でも同じ $e_{0}$ と局所消去を使う。

## 8. SH-07：質量と重力

固定した参照補間S、密度 $\rho_{0}$ 、層厚さに対し

$$
\mathbf M_e=\sum_\ell\int_{V_{0\ell}}\rho_0\mathbf S^{\mathsf T}\mathbf S\,\mathrm dV_0
$$

$$
\mathbf Q_g=\sum_\ell\int_{V_{0\ell}}\rho_0\mathbf S^{\mathsf T}\mathbf g_I\,\mathrm dV_0
$$

である。本候補の絶対位置・勾配座標ではSと参照体積がqに依存しないため $\mathbf M_{e}$ は定数。T=1/2 $\mathbf v_{e}$ ^T $\mathbf M_{e}$ $\mathbf v_{e}$ から、この要素の慣性項は $\mathbf M_{e}$ $\mathbf a_{e}$ となる。偏心質量や別の梁・拘束座標へ同じ結論を無条件に拡張しない。

非対称な層・密度では位置と勾配の質量連成が残り得る。零と仮定せず積分する。総質量・重心・慣性と剛体運動の運動エネルギーを独立に照合する。

参照実装は層ごとに三方向2点の積分を使う。UVLASではこれを初期候補とし、曲面・歪み・材料場・初期ひずみについて求積感度を確認する。公式一覧の「4 integration points」を体積積分全体の点数へ読み替えない。

## 9. SH-08：応力・ひずみ回復

回復で使う $e_{el}$ 、α、材料・層・評価点・初期ひずみは、内力計算と同じ定義で再評価する。出力に以下を区別して持たせる。

- 適合ひずみ、ANS/EAS後のモデルひずみ、指定初期ひずみ、弾性ひずみ。
- エネルギー共役の第二PK応力と、必要に応じて $\sigma =(1/detF) F S F^{\mathsf T}$ で得るCauchy応力。後者を出す場合、実際のFと仮定ひずみモデルとの関係を明記。
- 層ID・厚さ位置・材料基底・成分順・平均化の有無。

内力用のEASを出力だけで落とす、または初期ひずみを出力で再び差し引く、という不一致を許さない。仮定ひずみの応力を用いた断面力回復と一般化内力の整合も別ケースで確認する。

## 10. 要素処理の契約案

| 項目 | 仕様化する内容 |
|---|---|
| 不変入力 | $\mathbf q_{e0}$ 、層・材料・密度、初期ひずみ場、求積・参照幾何・変換・EASの定義 |
| 試行入力 | $\mathbf q_{e}$ [24]、必要なら $\mathbf v_{e}$ [24]。時刻依存材料入力は今回の必須範囲として追加しない |
| 出力 | $\Pi_{e}$ 、 $\mathbf f_{q}$ [24]、 $\mathbf K_{red}$ [24×24]、局所状態α、尺度付き局所残差・診断。質量はモデル構築時に評価可能 |
| 副作用 | Model・確定Stateを書き換えない。αや中間行列は呼出側Workspace／試行結果へ格納 |
| 呼出順 | 内力の直前呼出しを前提にせず、接線・応力の各要求が対応するq・世代を検証して再評価 |
| 成功条件 | 幾何・材料・局所解が適格、全要求出力が同じq・α・ $e_{0}$ に対応 |
| 失敗 | 要素・層・評価点・失敗分類を返す。未収束の値を有効な接線・応力として公開しない |

これは要素評価の責任・形状の具体化。C++の型・所有権・配列APIの確定は共通IF設計と合わせて行う。

## 11. 採用前の未完了項目と現在の判断

| 項目 | 判定・必要な作業 |
|---|---|
| 補間・参照写像・質量 | 式と固定版実装の対応を整理。曲面・剛体質量特性の具体ケースへ展開する |
| ANS | 座標順序の懸案に加え、ANS-V02で曲面の一定物理ひずみ再現性に反例。SH-04の単純な曲面適用は不採用 |
| EAS | 要素全体ζの層別案は不採用。層内座標・共通参照系での直交・線形局所消去を設計条件とした。要素全体の安定性・薄肉等は未確認 |
| 材料・初期ひずみ | 候補エネルギーを明示。測度・厚さ応力・自然配置型との比較はST-03で継続 |
| 接続・輸入 | 交差部材の勾配関係、梁接続、三角形、Nastran回転DOFとの変換は未定式化 |
| 証拠 | 数式整理とソース読解まで。シェル数値計算・独立FEM比較は未実施 |

判断：3423型は補間・質量を具体化できる候補として残す。ANS/EAS・曲面・接続・層別初期ひずみの課題が閉じていないため、主シェルとしての採用と計算本体の実装開始は保留する。3443／3833との比較も継続する。
