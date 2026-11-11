# 初期ひずみを含む平衡・初期化・微小応答の契約案

版0.2、2026-10-30。TM-01/03、ST-08/10、OP-08/09/10、IS-V04の具体化。構造保存系の限定経路を定義する。空力・後流・減衰・時間依存支持を含むJ1共同決定や、製品の一般非線形解法・時間積分の採用完了ではない。

## 1. 同じモデル・状態を引き継ぐ条件

参照幾何、層・材料軸・初期ひずみ、質量、支持・拘束、外力の種類と微分をモデル世代で識別する。平衡の変形形状を無応力参照へ置き換えず、初期ひずみを消去しない。固有値・過渡へ渡す状態は、平衡残差と拘束を評価した同じ状態とする。

初期化は次を区別する。

- 平衡から開始：指定した分岐・荷重・拘束の平衡を求め、残差・反力・接線を検査して受理する。
- 意図した非平衡から開始：入力q/vを勝手に静的緩和しない。運動方程式・拘束から初期加速度と必要な反力を求める。

平衡でない状態の接線を固有値計算へ渡し、平衡周りの安定性と表示しない。残差の小ささと物理的な安定性を区別する。

## 2. 制約された座標と接線

ANCF全体座標qの整合質量Mが定数でも、固定した非線形の運動学的制約を $\mathbf q=\mathbf g(\mathbf x)$ で表すと、低次座標xの質量は状態依存になり得る。参照・初期ひずみはgの変化で更新しない。

$$
\mathbf G=\frac{\partial\mathbf g}{\partial\mathbf x},\quad
\mathbf H_{ab}=\frac{\partial^2\mathbf g}{\partial x_a\partial x_b},\quad
\mathbf r_q=\mathbf f_{\mathrm{int}}-\mathbf f_{\mathrm{ext}},\quad
\mathbf r_x=\mathbf G^{\mathsf T}\mathbf r_q.
$$

$$
\mathbf K_x=\mathbf G^{\mathsf T}\mathbf K_q\mathbf G+
\left[\mathbf r_q^{\mathsf T}\mathbf H_{ab}\right]_{ab},\qquad
\mathbf K_q=\frac{\partial\mathbf r_q}{\partial\mathbf q}.
$$

非平衡反復中に第二項を捨てない。線形拘束写像ではH=0。乗数表示 $\mathbf r_q+\mathbf C_q^{\mathsf T}\boldsymbol\lambda=0$ の平衡では、対応する拘束曲率を含むLagrangian接線と整合させる。一般に構造内力の接線だけを全接線と呼ばない。

## 3. 非線形平衡の限定解法

今回の安定分岐の比較には、全ポテンシャルΦと一貫Kxを用いたNewton法・Armijo線探索を使う。残差と候補補正量の両方を検査し、未収束を受理しない。Newton方向が降下方向でない、接線解が不適格、幾何が不適格、線探索・反復上限に達した場合は未達で返す。

限界点・不安定分岐・弧長法・多重解探索の解法をこれで採用した扱いにしない。平衡分岐、初期推定、荷重経路、反復・線探索履歴を記録する。空力の非保存荷重がある系に同じポテンシャル線探索を無条件で流用しない。

全ポテンシャルの減少量が丸め誤差へ埋もれる場合、総量同士の差ではなく、同じ材料エネルギーと外力の増分を直接評価する。実際の候補qの差を用い、非線形写像のGΔxだけでエネルギー差を近似しない。判定方法と精度を記録し、許容差を緩めて未達を合格へ変更しない。

## 4. 運動方程式と初期加速度

$$
\mathbf M_x(\mathbf x)=\mathbf G^{\mathsf T}\mathbf M\mathbf G,\qquad
\mathbf h_x(\mathbf x,\dot{\mathbf x})=\mathbf G^{\mathsf T}\mathbf M\sum_{ab}\mathbf H_{ab}\dot x_a\dot x_b,
$$

$$
\mathbf M_x\ddot{\mathbf x}+\mathbf h_x+\mathbf r_x=0.
$$

非平衡の初期加速度は同じMx・h・rから求める。qとvの拘束適合を検査し、時間依存拘束の追加項は別途必要。MがANCFで定数という理由だけで縮約後のhを消さない。

## 5. 同じ平衡からの固有値・微小応答

静止平衡x*・速度0では、

$$
\mathbf M_*\delta\ddot{\mathbf x}+\mathbf K_*\delta\mathbf x=0,\qquad
\mathbf K_*\boldsymbol\phi_j=\omega_j^2\mathbf M_*\boldsymbol\phi_j.
$$

今回の保存・安定ケースでは対称一般化固有値問題を解く。質量正定値、固有対残差、質量正規化・直交性を確認する。負のω²を実数周波数へ絶対値で変換せず、零・近接・不確かさと区別する。減衰・追従荷重・空力のある問題へ対称性・実固有値を一般化しない。

微小応答の比較は同じM*・K*・支持を用いる。有限振幅の非線形応答との差は振幅を変えて調べ、時間離散化誤差とは別に評価する。線形化だけで有限振幅の妥当性を主張しない。

## 6. 時間離散化との接続

今回の検証用線形時間積分は陰的中点則。これは製品の積分方式の採用ではない。 $\mathbf z=(\delta\mathbf x,\delta\dot{\mathbf x})$ に対し、

$$
\underbrace{\begin{bmatrix}\mathbf I&-\Delta t\mathbf I/2\\\Delta t\mathbf K_*/2&\mathbf M_*\end{bmatrix}}_{\mathbf E_\Delta}\mathbf z_{n+1}
=\underbrace{\begin{bmatrix}\mathbf I&\Delta t\mathbf I/2\\-\Delta t\mathbf K_*/2&\mathbf M_*\end{bmatrix}}_{\mathbf A_\Delta}\mathbf z_n.
$$

離散固有対は $\mathbf A_\Delta\mathbf v=\mu\mathbf E_\Delta\mathbf v$ 。各振動モードの理想的な増幅率は

$$
\mu_\pm=\frac{1\pm i\omega\Delta t/2}{1\mp i\omega\Delta t/2},\quad
\omega_\Delta=\frac2{\Delta t}\arctan(\omega\Delta t/2).
$$

単位円上でも位相誤差がある。μの偏角/Δtを連続ωと同一視しない。符号・対数の枝・別名化・刻み収束を明記する。連成系のJ1と接続する際も同じ離散残差から線形化する。

## 7. 処理I/Fと失敗保証

| 処理 | 必須の入出力・保証 |
|---|---|
| 平衡試行 | モデル・荷重・支持世代、分岐・初期推定、残差尺度・補正尺度、Newton・線探索上限。候補x/q、同状態のr/K、反復履歴、受理／未達を返す |
| 平衡受理 | 幾何・有限性・残差・補正・拘束・反力・モデル一致を点検。受理前に確定状態を更新しない |
| 初期化 | 平衡／非平衡を明示。q/vの適合、初期加速度・反力・空力履歴の責任を識別。非平衡を暗黙緩和しない |
| 固有値 | 平衡状態・M/K・支持と入力世代、対象モード、残差・正規化・追跡法。連続／離散時間を区別 |
| 微小応答 | 同じ平衡とM/K、初期摂動・時刻・刻み・積分則。試行／受理を分け、平衡を無応力参照へ変更しない |

所有権・疎演算子・並列化・再開形式・一般の拘束／連成解法は未確定。残差評価・返す解の表現は [線形求解精度契約](linear-solve-accuracy-contract.md) と整合させる。

## 8. 証拠と残る範囲

[平衡・固有値・微小応答試験](verification/q9-prestress-equilibrium-experiment.md) は、Q9に埋め込んだ4座標のアフィン運動拘束を用いる。全81自由度の一般平衡や、実翼の固有振動の検証ではない。IS-V04の構成要素の証拠として用い、一般構造・連成・後流・フラッタ・再開の採用条件は残す。

INIT-D01：平衡・意図した非平衡を明示し、モデル世代と同じ状態のM/Kを受け渡す。INIT-D02：非線形運動拘束の接線曲率・慣性の速度項を含める。INIT-D03：連続固有値と同じ時間残差の離散固有値・位相を区別する。これらを今回の限定証拠に基づく技術判断とする。Newton線探索・中点則・DOP853を製品の標準解法へ採用した判断ではない。
