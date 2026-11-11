# Q9横せん断の局所混合投影候補 P-S36

2026-10-31。ST-04の比較候補を数式と処理契約へ落とす。**製品への採用決定ではない。** [成分診断](records/q9-distorted-energy-experiment.md) と [適合性反例](q9-distorted-shear-compatibility.md) を動機とする。まず固定参照・平坦・一定厚さ・等方材の範囲で構成する。

## 出典と独自導出の区分

[Bucalem & Bathe, Higher-order MITC general shell elements (1993)](https://doi.org/10.1002/nme.1620362109) の公開情報は、共変ひずみ成分の混合補間を9／16節点要素へ用いることを示す。今回のP-S36は下記のエネルギー直交投影を独自に定義したもので、公刊MITC9の忠実な実装とは称さない。論文全文の取得は今回のブラウザでは失敗したため、原式との一致や同論文の収束・安定性保証は主張しない。

## 対象と共変成分

Q9×厚さ二次の81座標、物理配置、質量、面荷重、支持を維持する。面内写像の参照ヤコビアンJ₂は物理XY行・親座標s,t列。横せん断の工学成分順は(XZ,YZ)。

$$
\widetilde{\boldsymbol\gamma}=\mathbf J_2^{\mathsf T}
\boldsymbol\gamma,\qquad
\boldsymbol\gamma=\begin{bmatrix}\gamma_{XZ}&\gamma_{YZ}\end{bmatrix}^{\mathsf T}.
$$

親座標は無次元なので、共変成分は長さを持つ。厚さζ=2z/hでLegendre多項式を使い、次の36次元空間Sを固定する。

$$
\widetilde\gamma_s\in Q_{1,2}(s,t)\otimes P_2(\zeta),\qquad
\widetilde\gamma_t\in Q_{2,1}(s,t)\otimes P_2(\zeta).
$$

Qijは各変数について次数i,j以下のtensor積。P₂は厚さ次数2以下。列順は、s成分／t成分、s次数、t次数、厚さ次数の順で厚さ次数を最速にする。Sは2×36。各列は単位の共変ベクトルとLegendre積。物理表現は次式。

$$
\mathbf V=\mathbf J_2^{-\mathsf T}\mathbf S.
$$

Vの単位は1/m、局所係数βの単位はm。単純な物理成分の点サンプリングや、積分点を減らす方式とは異なる。

## 局所投影と凝縮

μは参照等方材のせん断弾性係数。体積・材料・S・J₂は参照状態に固定する。

$$
\mathbf M=\int_{V_0}\mu\mathbf V^{\mathsf T}\mathbf V\,dV_0,
\qquad
\mathbf h(\mathbf q)=\int_{V_0}\mu\mathbf V^{\mathsf T}
\boldsymbol\gamma(\mathbf q)\,dV_0.
$$

$$
\mathbf M\boldsymbol\beta=\mathbf h(\mathbf q),\qquad
\overline{\boldsymbol\gamma}=\mathbf V\boldsymbol\beta.
$$

Mは36×36のGram行列。局所対角合同尺度化とCholeskyで解き、明示逆行列は作らない。階数・正定値が不成立なら候補の失敗とし、対角正則化を加えて通さない。

全エネルギーから元の横せん断項のみを外した部分をWnとする。投影後の縮約エネルギーは次式。

$$
\Pi_r(\mathbf q)=W_n(\mathbf q)
+\frac12\mathbf h(\mathbf q)^{\mathsf T}
\operatorname{solve}(\mathbf M,\mathbf h(\mathbf q))-W_{\rm ext}(\mathbf q).
$$

同値な局所混合汎関数は、βについて最大、qについて停留する次の鞍点形。

$$
\Pi_m=W_n+\boldsymbol\beta^{\mathsf T}\mathbf h
-\frac12\boldsymbol\beta^{\mathsf T}\mathbf M\boldsymbol\beta-W_{\rm ext}.
$$

投影は横せん断のμ重み付きL2直交投影で、連続の同じ積分なら投影後せん断エネルギーは元以下になる。柔らかくなること自体を精度向上の証拠としない。定数の物理せん断は、双一次J₂の共変変換が選んだSに入るため再現対象になるが、要素全体の零モードは別に検査する。

## 内力と一貫接線

hのJacobianをHとすると、固定Mについて次式が必要。

$$
\mathbf r=\nabla W_n+\mathbf H^{\mathsf T}\boldsymbol\beta-\mathbf f,
$$

$$
\mathbf K=\nabla^2W_n+\mathbf H^{\mathsf T}\operatorname{solve}(\mathbf M,\mathbf H)
+\sum_{a=1}^{36}\beta_a\nabla^2h_a-\frac{\partial\mathbf f}{\partial\mathbf q}.
$$

線形比較ではh=Hqなので最後のhの二次微分は0。非線形のGreenひずみへ拡張するなら、その項を省略しない。上式は固定参照・固定Mの場合の微分で、今回非線形の数値検証を済ませたことにはしない。参照状態やモードをqに合わせて更新する別方式には自動流用しない。

## 主要処理・不変条件

- prepare：参照幾何・厚さ・μ・求積・モード定義を受け、Mの尺度化因子・分解・版付きキーを生成する。幾何反転、非有限値、非正定値、尺度不成立を失敗として返す。部分的なキャッシュを確定しない。
- evaluate：同じキーの試行qからh,Hと必要な二次微分を評価し、β、縮約エネルギー・内力・接線を返す。入力qや確定状態を変更しない。βを前ステップの値で代用しない。
- βは各要素の代数的な局所量で、独立の慣性・履歴は持たない。構造の81座標と整合質量を維持する。動解析で節点自由度を静的縮約する判断とは区別する。
- Mの再利用キーは参照配置、厚さ座標、材料μ、求積、Sの次数・列順・版。いずれかが変われば無効。荷重だけの変更ではMを再利用可能だがhは再評価する。
- 定数ひずみの内力、剛体6モード、要素rank75、M正定値、Kの対称性、元と投影のエネルギー順序、局所停留と凝縮の一致を必須条件とする。

## 採用へ進めない境界

この候補はまだ、一般曲面・異方材・層別材料の法線／せん断結合、厚さ軸の再パラメータ化、一般形状での一定物理ひずみ、非線形接線、全体の動的安定性を閉じていない。μだけを変えて異方材へ適用しない。製品候補へ昇格するには全構成則との仕事共役性、既存EASの直交条件、材料軸と初期ひずみ、曲面パッチと層分割の整合を別途導く。

最初に [局所適格性ゲート](records/q9-shear-projection-gate-experiment.md) を実施する。通過してから既存の歪み板曲げと同じ8ケースを比較する。ゲート未達なら大域計算は実施せず、理由を保存してこの候補を見直す。
