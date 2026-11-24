# MAT1 と PSHELL の物性欄と材料参照を検査する

SHA-256 照合済みの Nastran デッキから MAT1 と PSHELL の欄を型付きで取り出す。結果は元のデッキを所有し、シェルの PID から PSHELL、各 MID から材料の宣言と元の欄をたどれる。[IO4-03](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran)

## MAT1 の値と空欄を保持する

E、G、NU、RHO、熱膨張係数、参照温度、減衰係数、引張・圧縮・せん断の応力限界を正確な十進値として保持する。E と G には少なくとも一方の指定が必要であり、明示値は正とする。NU は指数を含む元の値から開区間 (−1, 0.5) を照合する。MCSID の指定は座標系の宣言に照合する。補完前の E、G、NU と密度を保持し、3D 材料定数の確定と単位換算は後段で行う。

MAT1 と MAT8 の MID は同じ名前空間で重複を検査する。MAT8 の識別子と元のカードも保持し、異方性の参照を識別できるようにする。[MSC Nastran 2022.1 Quick Reference Guide、pp.2169–2172・2185](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

## PSHELL の材料用途を区別する

MID1 は膜、MID2 は曲げ、MID3 は横せん断、MID4 は膜と曲げの連成を表す。明示した正の番号は材料の宣言に照合する。MID2 の空欄、明示 0、平面ひずみを表す −1 は別々に保持する。平面ひずみには MID1 の MAT1 が必要である。MID2 が空欄または −1 のときは MID3 の指定を拒否する。明示 0 と正の MID3 の組合せは保持し、3D モデルでの扱いは後段で確認する。MID4 には正の MID1 と MID2 が必要で、いずれとも同じ番号にはできない。

基準厚さ、曲げ倍率、横せん断の厚さ比、面積当たりの非構造質量、応力計算位置 Z1・Z2 は入力のとおり保持する。明示した曲げ倍率と横せん断の厚さ比は正とする。材料番号の異なる組合せや MAT8 を含む組合せを識別して渡す。厚さと Z1・Z2 の既定値、せん断・曲げ倍率や異なる材料の 3D モデルとの等価性、解析で使う減衰と強度評価は、後段の変換と認定で扱う。[MSC Nastran 2022.1 Quick Reference Guide、pp.2901–2905](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

## 検証

取得済みのスナップショットから、節点を持つシェルと PSHELL、MAT1 の定義まで順にたどった。small・large・free の継続、空欄と明示値、指数を含む元の十進値と欄位置を照合した。材料番号の重複・未解決参照、用途の組合せ、NU の境界、予約欄、確保量と作業量の境界、取消を試験した。

Windows の MSVC 19.51.36260.0 Debug と Linux の GCC 14.2.0 Debug（ASan・UBSan）で、Nastran の対象5件がそれぞれ通過した。材料カードは各5ケース・295アサーション。Linux の clang-tidy 21.1.1 で本体の通常構成・確保計測構成と試験を検査した。
