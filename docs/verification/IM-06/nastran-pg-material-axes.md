# IM-06 Nastran と Case の層材料軸を照合する

Case と認証済みの Nastran 平面の対応を使い、一定全厚が確認された各元面で、Case の各層の材料軸と元の材料軸を照合する。層面の法線には厚さ勾配の影響があるため、材料軸の準備では変厚の面を適用外とする。[MAT-D05 材料軸の契約](../../specs/structure/material-layer-adoption-contract.md)

## 同じ向きの法線と軸で角度を照合する

THETA は元の G1→G2 を基準とした度の回転である。基本座標系の MCID は基本 x 軸を、基本 RID に直接置いた CORD2R の MCID はその x 軸を元面に射影する。PCOMP では元 ply の角度を加える。元法線と Case 法線が反対なら、Case 法線を軸に採り、角度の符号を反転する。これにより第一軸の世界方向を保ち、第二・第三軸が法線の向きに従って変わる。[MSC Nastran 2022.1 Quick Reference Guide、CQUAD4 pp.1609–1611、CTRIA3 pp.1667–1669、PCOMP pp.2764–2768](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

一定の projected 軸では、Case の方向を同じ層面に射影する。元と Case の角度差が90°の整数倍なら、元の基準方向をその分だけ回し、外積が零かつ内積が正かを厳密に調べる。両側の角度が同じ30°なら、差が厳密に零になる。射影と四分回転は正の倍率を残した有理ベクトルで表せるため、正規化の丸めは生じない。一定の explicit 軸では右手直交基底を厳密に検査し、Case 法線と第一軸を照合する。幾何変換の距離・法線予算はこの軸照合の許容差に使わない。

結果は元の EID・PID・Case LayerId、要素欄と ply 角の欄、世界座標の基準方向・Case 法線の有理ベクトル、厳密な度の角度を所有する。ラジアン表記は角度零を扱う。この方法で差を認定できない一般角、空間変化する Case 軸、円筒・球 MCID、入れ子の CORD2R は適用外を返す。

## 実デッキと独立な軸方向を照合する

四辺形の元 G1→G2 を +x とし、THETA 30° に二層の ply 角 +60° と −30° を加えた。Case の第一軸を下層 +y、上層 +x とすると合格し、上層だけを +y にすると第2層と元 ply 欄を報告した。もう一つの対照では元20°＋10°とCaseの30°が相殺することを確かめた。

元リングを逆にした対照では、元の基準方向が +y、元法線が −z になり、90°で第一軸は +x を向く。Case 法線 +z に対しては合成角 −90°を保持した。平行移動を持つ CORD2R の x 軸が +y の対照では、明示基底 (+y, −x, +z) と照合した。射影すると零の方向、差が四分回転にならない角度、変厚、ちょうどの件数上限と一件不足、取消も検査した。

## 検証に使った構成

対象ソースは `8c8f082768bc156fbf7a7eeabb7afa323a028154`。Windows 11 の MSVC 19.51.36260 Debug と Release、Ubuntu 24.04 の GCC 14.2.0 Debug（ASan・UBSan）で、`io.external_pg_regional_thickness`、`io.external_planar_pg`、`io.nastran_laminate_cards`、`io.reference_input` が各4件合格した。各構成で、実ファイルの対応表試験は1404アサーションだった。

隔離コピーで反向きの角度符号をそのまま使うと、90°の反向き対照が失敗した。元へ戻すと、実ファイルの対応表試験の1404アサーションが再び合格した。変更した翻訳単位の clang-tidy と clang-format 21.1.2 も通過した。
