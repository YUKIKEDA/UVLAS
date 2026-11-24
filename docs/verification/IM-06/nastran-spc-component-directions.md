# IM-06 NastranのSPC成分方向

対象ソースは **0382e59862106d597dc6dacd22aaf4176763df15**。[IO4-03](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran系の限定読込み) の認証済みデッキから、既に読んだSPCまたはSPC1の一つの元項と対象GRIDを選ぶ。並進1–3、回転4–6のうち、元項で指定した成分について基本Cartesian系の単位方向を返す。同じ軸の並進と回転は同じ方向を使う。[MSC Nastran 2022.1 Linear Static Analysis User’s Guide, Chapter 3, p. 63](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2022.1_Linear_Static_Analysis_User_Guide/raw/resource/enus/MSC_Nastran_2022.1_Linear_Static_Analysis_User_Guide.pdf)

## 基底と原デッキ

CPから基本系に移した点を、CDの定義直交枠へ戻す。点の変換では原点を使い、ベクトルには回転だけを適用する。CPとCDが同じ場合は、その枠内のCartesian位置を直接使う。直交CDではGRID位置にかかわらず定義軸を使う。

円筒・球のCD内の点を $(x,y,z)$ とし、$\rho=\sqrt{x^2+y^2}$、$r=\sqrt{\rho^2+z^2}$ とする。円筒の三軸は $(x,y,0)/\rho$、$(-y,x,0)/\rho$、$(0,0,1)$ の順。球の三軸は $(x,y,z)/r$、$(xz,yz,-\rho^2)/(r\rho)$、$(-y,x,0)/\rho$ の順。正の半径が有向丸めで成立する領域で、この局所基底を基本系へ回転する。結果の区間端点は、無次元の厳密有理数である。

円筒CDの3・6のみを選ぶ場合は、極軸でも一意な定義直交枠の第3軸を返す。位置に依存する基底を求める場合、極軸が確定していればUnsupported、半径の包絡が零を含む場合はPrecisionUnresolvedとして、元GRIDのCD欄を示す。商用版の極軸近傍での自動切替を再現するには、切替の条件を別に確定する必要がある。[MSC Nastran 2022.1 Linear Static Analysis User’s Guide, Chapter 8, p. 320](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2022.1_Linear_Static_Analysis_User_Guide/raw/resource/enus/MSC_Nastran_2022.1_Linear_Static_Analysis_User_Guide.pdf)・[MSC Nastran 2021 Reference Guide, Grid Point and Coordinate System Definition](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2021_Reference_Guide/raw/resource/enus/MSC_Nastran_2021_Reference_Guide.pdf)

SPCの元項とGRIDの所有者は親、参照宣言、元バイトの指紋で照合する。SPC1のTHRUは範囲内の既存GRIDに適用する。返値は元セットID、GRID、CD、カードと元欄、指定値の十進字句を所有する。Case Controlの組合せ選択、指定値の単位換算、ANCFへの拘束組立てはその先の工程となる。

## 独立対照

認証した実デッキからSPCとSPC1を読み、CDの長方形・円筒・球を選んだ。長方形CDは平行移動を含む45度の回転で、成分1と4の方向が同じになる。円筒CDでは別のCPと移動したCD原点を使い、半径、接線、軸を独立に照合した。球の対照点はCD内で $(1,2,2)$ であり、三方向は $(1,2,2)/3$、$(2,4,-5)/(3\sqrt{5})$、$(-2,1,0)/\sqrt{5}$ になる。有理の両端から平方と符号を比較し、向きを照合した。

確定した球の極、丸め包絡で極を判定できない円筒点、確定した極と包絡で未確定な極における円筒軸成分、THRU内の存在しないGRID、異なる元バイトの所有者、作業上限ちょうどと1不足、親の破棄後の所有を確認した。

## 実行結果

固定ソースの刻印と変更なしの状態を照合した。MSVC 19.51.36260 DebugとUbuntu 24.04.5のGCC 14.2.0 Debug（ASan・UBSan）で io.nastran_coordinate_cards、io.nastran_shell_cards、io.external_planar_pg の各3件が合格した。座標の主試験は両環境で9ケース・486アサーション。Linuxのclang-tidyとWindowsのclang-format 22.1.3を通過した。

Linuxの分離コピーで球の子午線方向を求める外積の順序を逆にすると、独立対照の9アサーションが失敗した。元の式に戻して、座標試験の9ケース・486アサーションが再び合格した。
