# IM-06 Nastran の基本系座標の包絡

対象ソースは **a24d72ecbb8ba4b04ff068dffb23cecbecbfd016**。[IO4-03 の限定読込み](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran系の限定読込み) で得た型付きデッキから、指定したGRIDの基本Cartesian系SI座標を閉区間で返す。

## 入力と計算

GRIDのCPとCORD2のRIDをたどる。CORD2の定義三点 $A,B,C$ は参照先RIDの座標種別に従って変換する。円筒の $(r,\theta,z)$ は $(r\cos\theta,r\sin\theta,z)$、球の $(r,\theta,\phi)$ は $(r\sin\theta\cos\phi,r\sin\theta\sin\phi,r\cos\theta)$ になる。GRIDのCDは変位と自由度の座標系であり、位置にはCPを使う。[MSC Nastran 2022.1 Linear Static Analysis User's Guide, Chapter 3](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Linear_Static_Analysis_User_Guide/raw/resource/enus/MSC_Nastran_2022.1_Linear_Static_Analysis_User_Guide.pdf)

ネイティブの度は元の十進数から有理数として読み、厳密に360度で剰余を取り、$\pi$ の倍数にする。軸は $e_3=(B-A)/\lVert B-A\rVert$、$e_2=(e_3\times(C-A))/\lVert e_3\times(C-A)\rVert$、$e_1=e_2\times e_3$ で作る。各段で変換した位置はCartesianのまま親へ渡す。親が円筒系や球系であっても、中間結果を角度として読み替えない。長さ単位の倍率は基本系に移した位置へ適用する。三角関数と平方根には既存の256ビット有向丸めを使い、下端と上端を分母が2の冪の厳密有理数として保持する。[MSC Nastran 2022.1 Getting Started Guide, Chapter 5](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2022.1_Getting_Started_Guide/raw/resource/enus/MSC_Nastran_2022.1_Getting_Started_Guide.pdf)

親を再取得せず、結果は元の参照位置、ハッシュ、GRID番号、元CPと座標欄の位置を所有する。軸または面を作るノルムの二乗の上端が零ならGeometryInvalid、零を含み正の上端が残るならPrecisionUnresolvedを返す。選択点数、作業計数、有理数の寸法、確保容量を検査し、失敗時は出力用の確保を解放する。

## 対照と適用範囲

認証したデッキから異種のR・C・Sの親子関係、非有理の軸正規化、非常に大きな負角、90度をわずかに越えた十進角を読み込んだ。軸と位置の正負を別計算で決め、平方根を含む値は有理の両端を二乗して包絡を照合した。90度近傍では符号と幅も照合した。CPとCDの違い、厳密に縮退する軸と包絡で判定できない軸、親の破棄後の所有、資源不足後の再試行を実入力で確認した。

この入口が返すのはGRID点の座標包絡である。既存PGの厳密照合には [有理入力の経路](nastran-planar-placement.md) を使う。座標包絡から面やZOFFSを含む距離と法線の誤差を評価するには、[FD-03の参照幾何とIO4-04](../../specs/foundation/io-execution-adoption-specification.md#io4-04-gmsh系の限定読込み) の変換予算と独立観測による判定が必要になる。

## 実行結果

MSVC 19.51.36260 Debug と Ubuntu 24.04.5のGCC 14.2.0 Debug（ASan・UBSan）で固定ソースの刻印を照合した。両環境で io.nastran_coordinate_cards、io.nastran_shell_cards、io.external_planar_pg の3件が合格した。座標試験はそれぞれ7ケース・334アサーション。Linuxのclang-tidyと Windowsのclang-format 22.1.3を通過した。

Linuxの分離コピーで球座標の方位成分のsinとcosを交換すると、対応する試験の4アサーションが失敗した。元の式に戻すと、座標試験の7ケース・334アサーションが再び合格した。
