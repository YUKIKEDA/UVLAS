# IM-06 Nastran の厳密な位置座標と平面オフセット

対象ソースは `b136310bf036557cf9444cb3f0de4335b62f7601`。[IO4-03 と IO4-04](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran系の限定読込み) に従い、型付きデッキから基本CartesianのGRID位置と要素参照面を作り、[既存PG平面の照合](external-planar-pg.md)へ渡す。

## 採用式と検査

CORD2R の三点 $A,B,C$ を長方形の参照系RIDで読む。$e_3=(B-A)/\lVert B-A\rVert$、$e_2=(e_3\times(C-A))/\lVert e_3\times(C-A)\rVert$、$e_1=e_2\times e_3$ とする。各GRIDに対して $A+e_1x+e_2y+e_3z$ を適用し、RIDをたどって基本系へ移す。元の有向ベクトルで軸の符号を決め、分子と分母の整数平方根を二乗して有理の正規化を証明する。SI倍率は基本系に移した位置へ一度適用する。退化する三点、正規化に非有理数を要する三点、曲線座標系への参照は別の失敗値とし、元のCP/RIDまたは三点の座標欄を示す。[MSC Nastran 2022.1 Quick Reference Guide](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

数値ZOFFSは元のG1、G2、G3の順で作る右手法線に沿う符号付き距離である。非零なら $P_i^\mathrm{ref}=P_i^\mathrm{GRID}+ZOFFS\,n$ とし、四角形では元の第4点も平面上にあることを検査する。PGのチャートに宣言した向きから元のZOFFS符号を変更しない。認定済み結果は元GRIDのSI位置、移動した位置、元CP、元ZOFFSと欄の位置を保持する。PCOMPの層底位置Z0は別の量であり、既存の材料、厚さ、質量、接続をこの幾何比較から置換しない。[MSC Nastran 2022.1 Getting Started Guide, Chapter 6](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2022.1_Getting_Started_Guide/raw/resource/enus/MSC_Nastran_2022.1_Getting_Started_Guide.pdf)

## 独立対照

入れ子のCORD2Rへ与えた局所座標 $(x,y,z)$ は、元単位mmで $(1220+4x/5-3y/5,-540-3x/5-4y/5,700-z)$ になる。負の軸と $3/5,4/5$ を使う回転を別の有理計算で確認した。正順の四角形と三角形には −300 mm、逆順の四角形には +300 mmを与え、元平面 $Z=700$ mmから同じ $Z=1000$ mmの参照平面を得た。正順の符号を +300 mmに変えると $Z=400$ mmになり、元ZOFFS欄を報告した。第4点が平面を外れる四角形、非有理の正規化、共線の三点、円筒座標の負対照も元デッキから読んだ。

斜め法線 $(0,3/5,4/5)$ に +0.5 mを適用する対照では、元点 $(1,4/5,-3/5)$ mが $(1,11/10,-1/5)$ mへ移る。

## 適用範囲と残る設計

入力の有理座標を既存PGへ許容差0で照合する範囲を確認した。IO4-03の非有理CORD2R、および位置を円筒・球から変換するCORD2C/Sを扱うには、入れ子の回転と三角関数から生じる量の代数表現または有向丸め付き包絡を決める必要がある。包絡の重なりだけで厳密な一致とは判定できない。[FD-03の参照幾何](../../specs/foundation/case-schema-fields.md#fd-03-幾何と外部source) とIO4-04に従う近似の受入には、距離・法線の予算と独立観測を使う。形状の調整、積層の側やZ0、接続点や質量の変換を同じ物理モデルとして認定するのはR-05の範囲。

## 実行結果

MSVC 19.51.36260 Debug と Ubuntu 24.04.5 の GCC 14.2.0 Debug（ASan・UBSan）で、`io.external_planar_pg`、`io.nastran_coordinate_cards`、`io.nastran_shell_cards`、`io.gmsh_surface_topology` の各4件が合格した。主試験は両環境で4ケース・105アサーション。固定ソースの刻印と未変更状態を照合し、Linuxのclang-tidy 21.1.1 と clang-format 21.1.2 を通過した。

Linuxの分離コピーで元法線の世界X/Y成分を零にすると、斜め法線の照合がGeometryInvalidで失敗した。元に戻して再実行した主試験は合格した。
