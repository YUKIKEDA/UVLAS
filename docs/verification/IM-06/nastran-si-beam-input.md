# IM-06 Nastran CBEAM・PBEAMのSI入力

対象ソースは **b2aba145d884899f3c1b004b2c1dc5c229e21cb3**。[IO4-03](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran系の限定読込み) の認証済みデッキのCBEAMから、同じPBEAMとMAT1を参照する。せん断中心軸の端、向き、長さ、材料と宣言断面のSI物量を所有する。デッキの宣言位置、指紋、元欄を保持する。

## 軸、端部、偏心

GRIDのCPから基本系位置を包絡で求める。向きベクトルはGAからG0への差、基本系の実数成分、またはGAのCDの実数成分を使う。CBEAMのOFFTでGと指定したWA・WBは各GRIDのCDに従う。Oと指定した成分の基底は、偏心を加える前のGAからGBへの向きをx軸とする。いずれも成分ベクトルを基本系へ回転してから、元GRIDに足す。

偏心後の二端間を最終要素のx軸とし、その単位方向と向きベクトルの外積を正規化してz軸を得る。y軸はzとxの外積とする。長さは偏心後の両端間で測る。元のGA/GBとせん断中心軸の端は別々に保持する。使用する軸の確定した縮退と、包絡だけでは非零を証明できない場合を区別する。[MSC Nastran 2022.1 Linear Static Analysis User's Guide, Chapter 4, pp. 119–120](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2022.1_Linear_Static_Analysis_User_Guide/raw/resource/enus/MSC_Nastran_2022.1_Linear_Static_Analysis_User_Guide.pdf)

## 材料と断面

MAT1の二つの弾性定数から $E=2G(1+\nu)$ で残る一つを厳密に求める。三つを供給する場合は、同じ等方材としてこの式に厳密に一致することを確認する。RHOが明示された場合の線質量は $\rho A+\mathrm{NSM}$ とする。RHOの空欄は保持し、その場合の線質量は未確定として返す。

断面積の単位はm²、断面二次モーメント・積・ねじり定数の単位はm⁴、NSMはkg/mとする。端Bの空欄は端Aの同じ欄を使う。途中ステーションの空欄は、端Aと端BをX/XBで線形補間する。明示した途中の欄はその値を使い、結果は元の宣言順で保持する。二つの二次モーメントと積がつくる行列の正定性も検査する。せん断係数K1・K2は明示値か規定値1を使う。個々の $EA$、$K_iGA$、$GJ$、$EI_i$、$EI_{12}$ を宣言面番号のスカラーとして保持する。[MSC Nastran 2022.1 Quick Reference Guide, pp. 2694–2698](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

NSMの重心Mと中立軸Nは、せん断中心からの局所y・zの距離として別々に保持する。NSIはkg mで保持する。有限のせん断入力として表せないK1/K2の明示0、非零テーパ補正、端解放、反り指定、明示した応力回復は元欄を示して拒否する。断面の包絡寸法、完全な重心慣性、原点を移して組み立てる6×6の構成則、Caseへの対応と実組立ては後続で扱う。[採用梁仕様](../../specs/structure/beam-adopted-specification.md#4-断面入力原点反り)

## 独立対照

元のGAとGBを $(1,2,3)$ mと $(4,2,3)$ mに置いた。同じ向きをG0、基本系の成分、円筒CDの成分で指定した。GとOの偏心を使い分け、三経路とも偏心後の端 $(1,3,3)$ mと $(5,6,3)$ mを得た。独立な3–4–5三角形を使い、長さ5 m、x軸 $(4/5,3/5,0)$、y軸 $(-3/5,4/5,0)$、z軸 $(0,0,1)$ を包絡と照合した。

元デッキのmmとgから、端Aの $EA=2.1\times10^8$ N、せん断の2値 $6.3\times10^7$ N・$4.2\times10^7$ N、$GJ=2.52\times10^5$ N m²を別計算した。端Bを途中より先に宣言し、元順のまま既定値を照合した。途中の面積は0.0015 m²、明示のI1は $2.0\times10^{-6}$ m⁴、線質量は14.2 kg/mとなる。断面二次モーメントと整合しない積、等方定数の不一致、非零反り、回復、共線性、元デッキ不一致、資源上限も確認した。

## 実行結果

固定ソースの刻印と変更なしの状態を照合した。MSVC 19.51.36260 DebugとUbuntu 24.04.5のGCC 14.2.0 Debug（ASan・UBSan）で、座標、シェル、CBEAM、PBEAM、PG平面照合の各5件が合格した。SI梁の絞り込み試験は両環境で3ケース・312アサーション。Linuxのclang-tidyとWindowsのclang-format 22.1.3を通過した。

Linuxの分離コピーで最終要素のx軸を偏心前のGRID間方向に置き換えると、三経路の長さと姿勢の計15アサーションが失敗した。元の式に戻して、SI梁の3ケース・312アサーションが再び合格した。
