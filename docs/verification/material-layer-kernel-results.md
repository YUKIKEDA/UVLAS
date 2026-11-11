# 1.2 材料点・層契約の追加代数照合

2026-11-03。[材料・層採用契約](../material-layer-adoption-contract.md)の新しい境界を、[事前仕様](material-layer-kernel-protocol.md)と[準備ハッシュ](material-layer-kernel-preparation.json)を固定して照合した。[結果](material-layer-kernel-result.json)の6群は全て成立。既存のQ9材料接線試験は再実行していない。

| 検査群 | 最大規定誤差 | 結果 |
|---|---:|---|
| Mandel_engineering_bridge | 1.040e-16 | 合格 |
| curved_reference_fiber_axis | 7.536e-16 | 合格 |
| parallel_fiber_negative_control | 0.000e+00 | 合格 |
| normal_stress_Schur_and_recovery | 2.157e-16 | 合格 |
| layer_moments_vs_tensor_quadrature | 2.711e-19 | 合格 |
| homogeneous_layer_subdivision | 1.770e-16 | 合格 |

Mandelと工学成分の仕事・剛性・compliance、参照層法線と繊維方向からの右手軸と回転変換、法線平行繊維の拒否、SzzだけのSchur消去と別compliance回復、層別のA/B/D・初期項とテンソル直接求積、均質層の3分割不変性を確認した。Schur消去を省略した剛性との差も負の対照として検出した。

検査材料は合成値。数学的なエネルギー・回復・積分の整合を示すもので、実材料の有効範囲、縮約カーネルと要素厚さ運動学の整合、一般構造の実務精度を検証した結果ではない。

**閉じた範囲**はMAT-D04〜07の材料点の構成式・成分・参照軸・層帰属・必要微分。1.2全体には、採用要素の層別積分／回復、厚さ法線消去の要素アダプタ、等価特性入力、実材料データの適用範囲が残る。主シェルの成否と無関係に再利用する設計判断として保存する。
