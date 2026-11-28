# Case で指定するシェルの全厚側面への面力

[FD-05](../foundation/case-schema-fields.md#fd-05-接続支持質量荷重) のシェル `Load.traction` で、`Region` を `kind: edge`、`side: interior` とした場合、指定された物理外周辺に沿う側面の全層・全厚を選ぶ。荷重の次元は Pa である。面の選択は `area`、`frame`、空間場や時刻倍率の指定に共通とし、`boundary` は Case の参照幾何にある名付き境界の番号を使う。 [採用した分布荷重と仮想仕事](structural-adopted-specification.md#3-荷重の測度追従接線)・[ケースの参照と拘束](pg-case-chart-contract.md)

## 名付き境界から全層の側面を選ぶ

型付き SI Case の Shell ID と Geometry ID を、その Shell が実際に持つ認定済み PG パッチと照合する。`boundary` の元線分ごとに、同じ [共通分割](pg-integration-partition-contract.md#3-体積辺上下層は同じ片を参照する) に由来を持つ原子辺を取得する。外側に物理片が存在せず、元領域の境界由来を持つ原子辺だけが `LateralFace` を持つ。元の各線分を有理数の端点で区切り、原子辺の合併が端から端まで厳密に覆い、正の長さの重複がないことを確認する。対象の原子辺ごとに、材料層の全区間に対応する `LateralFace` を一度ずつ展開する。

選択順は Shell の Case 順、Load の Case 順、名付き境界の元線分順とその始点からの位置、層の下から上、辺節点、深さ節点とする。元線分の向きとパッチ側面の向きが逆なら、辺の局所座標を $u\mapsto 1-u$ とする。同じ点を共有する隣接面片では、面片ごとに積分点と重みの所有を記録する。異なる Load ID が同じ面へ明示した力は、由来を保って加算する。

Nastran と Gmsh では、認証済み参照から準備した同じ Shell の物理境界へ、型付き Case の同じ荷重を適用する。荷重の出所は Case の Load ID と名付き境界である。 [外部面と Case の対応](../foundation/external-mesh-sidecar-contract.md)

## 参照の面積当たりの力を元の列へ返す

最初の実装入口は `area: reference`、`frame: world`、`spatial.kind: constant`、`multiplier.kind: constant` の組を受理する。この入口で適格な別のタグを指定すると `Unsupported` を返す。一定の世界 XYZ 面力 $t_{0,a}$ [N/m²] と無次元倍率 $m$ は同じ SI Case から得る。状態に依存する測度や方向を使う場合の作用と接線には、構造採用仕様の現在面積と追従荷重の規約を適用する。

各面片の規則は、辺と層内深さに対する二つの正の単位区間則から作る。各則は空でなく、節点は $[0,1]$ で厳密昇順、各重みは厳密に正、総和は厳密に1とする。すべての面片に同じ則を展開する。`PgPreparedPatch::evaluate_reference` に `PgReferenceLateralPoint` として面片、辺座標、層内深さを渡し、同じ所有点から無次元の $N_i$ と `PgReferenceLateralMeasure.area.squared_norm` を得る。[$A_0$ [m²] とその二乗ノルム](pg-reference-geometry-contract.md) には辺と層区間の写像の因子を含む。有限則による外力の $(i,a)$ 成分は次で与える。

$$
f_{i,a} = \sum_p N_i(p)\,t_{0,a}\,m\,w_{p,u}w_{p,v}\sqrt{A_0(p)\cdot A_0(p)}.
$$

成功結果は、拘束消去前の条件付け済み参照基底の列順、世界 XYZ 順の $r\times3$ [N] の所有包絡である。元残差には $-f$ を加える。固定参照面積と固定世界方向を使うこの経路では、$\partial f/\partial q$ は厳密に零となる。一定倍率でも評価時刻は有限でなければならない。拘束と反力は、この荷重を同じ元の自由度へ加えた平衡から求める。 [元の構造残差](structural-adopted-specification.md#2-構造の共通残差と符号)

## 所有と失敗を記録する

準備結果は Case、Load、物理境界、パッチと基底の識別、元線分と側面片、層、座標、厳密重み、点の順を保持する。各作用は同じ所有結果を読み、失敗しても元の Case、基底、試行を保持する。欠けた面片や途中までの荷重を成功として公開しない。

同じ Case 内の対象や幾何の不一致、未知の辺番号、節点と重みの不備は `InvalidInput`。非有限の時刻は `Nonfinite`、別 Case・パッチ・基底の準備結果または世代の取り違えは `StaleGeneration` とする。物理外周の被覆欠け・重複や、側面の面積二乗が非正と確定した場合は `GeometryInvalid` とする。面積二乗 $a+b\sqrt g$ の符号は参照幾何の有理比較で決める。平方根を含む計算は入力精度と算術の丸めを外向きに囲う。正値や指定された出力幅を証明できなければ `PrecisionUnresolved`、要求した表現を扱えなければ `PrecisionUnsupported` を返す。出力幅と連続面に対する有限則の誤差は別々に照合する。

Load 数、元線分・原子辺・層・節点数を検査し、面片数と二つの則の積を確保前に検査する。全点と $r\times3$ の元列に対する算術、平方根精密化、計画領域と所有配列、累積確保を予算に含める。checked 上限を超えると `ResourceExceeded`、取消は `Cancelled` とする。失敗には可能な範囲で元の Load ID、境界番号、面片または層を付ける。

## 全厚と載荷位置を独立に検査する

独立参照は、全厚2.54 mm、載荷端 $y=4\,\mathrm m$、幅方向 $z\in[3,5]\,\mathrm m$ の矩形とする。面積と面心から求めた合力、移動原点回りのモーメント、剛体仕事を Cartesian 座標の有理積分で照合する。外周・層を分割しても物理観測が一致すること、載荷端と厚さ方向を変える対照が失敗することを確かめた。 [独立式と製品試験の受入条件](../../verification/IM-02/pg-lateral-traction/verification.md)

製品試験は型付き Case を準備して元列の荷重を評価し、全成分を独立な元基底の面上積分と照合する。準備した基底への変換は参照側で最後に適用する。剛体並進と回転、軸方向と深さ方向の非剛体仕事を、同じ元列の双対積から回復する。元線分の向きを逆にした入力、非等分の側面と層、異なる Load ID の重ね合わせ、所有・資源・時刻・取消の失敗も含める。Case の面力を平衡・反力・固有対へ渡す検証と、三つの入力経路を通した物理比較は IM02-0 と IM06-0 の R-05 に続く。
