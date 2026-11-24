# 凍結R4の全入力でPGの参照位置とJacobianを照合する

R4の品質系列24配置と追加係数系列4配置を使った。柱面、楕円面、ねじれ面、鞍面、混合4次面、歪みを加えた混合4次面、11項の任意係数面について、原本が指定する二つの格子と二つの厚さをそろえた。必要な配置では $H=h(1+s/10+t/20)$ とした。原本の7形状、格子と厚さの組合せを保持した。 [品質系列の入力](../../../rd/structure/main-shell/paired-geometry/experiments/check_mass_paired_quality.py)・[追加係数の入力](../../../rd/structure/main-shell/paired-geometry/experiments/check_mass_paired_coefficient_input.py)

製品には、物理座標の多項式を有理係数で $x=s+a t$ に展開したチャートを渡す。独立な参照は、物理座標のまま高さ $f$ と勾配を評価し、参照面に $\zeta H/2$ と単位法線の積を足す。Jacobianはこの写像の三方向に複素ステップを与え、単位法線を含む写像全体から求める。先頭と末尾の積分領域で、位置3成分とJacobian9成分の差を各 $2\times10^{-14}$ 以下と確認した。 [製品試験](../../../tests/structure/pg_reference_geometry_test.cpp)

この系列で照合した量は参照位置とJacobian。切断した後退・テーパの体積と全外面荷重はR5の系列で別に照合した。 [R5の独立荷重](pg-r5-swept-taper-reference.md)
