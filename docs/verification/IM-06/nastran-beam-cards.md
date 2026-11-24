# CBEAM の節点、向き、端部の指定を読む

SHA-256 照合済みの Nastran デッキから CBEAM を型付きで取り出し、参照先 PBEAM の宣言を照合する。結果は元の宣言位置と正確な十進値を所有する。元デッキは読取り専用で受け取り、梁以外のカード読取りにも続けて使える。[IO4-03](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran)

## 梁の両端と向きの指定を区別する

GA と GB は別々の GRID を参照する。整数で指定された G0 は GA・GB と異なる GRID に照合し、実数で与えた成分とは区別する。成分を指定する場合は X1〜X3 の元欄と空欄を保つ。BEAMOR は採用カードの範囲外なので、向きの指定が一つもなければ拒否する。PID の空欄は EID を用いて PBEAM を照合する。

OFFT の空欄は GGG の規則を表す。明示指定は QRG 表の8通りに限る。X1〜X3 の成分には GA の変位座標系か基本座標系を、端部ごとのオフセットには各 GRID の変位座標系かオフセット座標系を区別して保持する。採用表にない指定は元の欄を示して拒否する。[MSC Nastran 2022.1 Quick Reference Guide、pp.1470・1472–1473](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

## 端部の解放、ずれ、反り変数を記録する

PA と PB は両端の要素成分の解放を表し、各端で1〜6のうち重複しない最大5成分を保持する。WA と WB は GRID からせん断中心軸の端までの成分で、明示した零と空欄を区別する。

SA と SB の空欄、明示0、正の GRID 参照を区別する。GRID の場合は外部成分 T1 の反り変数を指定する。SPOINT は採用範囲外である。ここでは欄と参照を確定し、座標変換、断面とせん断中心の物性、端部解放と反りを含む物理モデルの適合性、K6 と単位換算は後段で扱う。[MSC Nastran 2022.1 Quick Reference Guide、pp.1471–1473](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

## 検証

small・large・free の CBEAM を同じ入口で読み、両端、G0 と成分指定、8通りの OFFT、空の第1継続を経た SA と SB を照合した。重複した PBEAM、参照切れ、値域外の解放、未知の指定、予約欄、確保量と作業量、取消を試験した。成功と失敗の後に元デッキを読めることと、元デッキを破棄した後にも梁の型付き値を参照できることを確かめた。

Windows の MSVC 19.51.36260.0 Debug と Linux の GCC 14.2.0 Debug（ASan・UBSan）で、Nastran の対象6件がそれぞれ通過した。梁カードは各5ケース・304アサーション。Linux の clang-tidy 21.1.1 で本体の通常構成・確保計測構成と試験を検査した。
