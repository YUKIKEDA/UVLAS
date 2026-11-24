# PBEAM の断面と軸方向ステーションを読む

SHA-256 照合済みの Nastran デッキから構造用 PBEAM を型付きで読み取る。CBEAM の検査結果と材料の読取り結果は、元の宣言位置、親マニフェストの識別情報と版、デッキの SHA-256 が一致する場合に組み合わせる。CBEAM の参照番号から同じ PBEAM をたどり、MID から MAT1 の宣言を参照できる。両方の入力は読取り専用で使用する。[IO4-03](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran)

## 断面と応力回復の指定を区別する

端 A の断面積と二つの断面二次モーメントは正の実数を必須とする。断面相乗モーメント、ねじりのパラメータ、単位長さ当たりの非構造質量と、応力回復点 C〜F の局所 y・z 座標も元の欄で保持する。端 A のねじりパラメータ J は空欄と零を保持し、負値は元の欄を示して拒否する。

ステーションを入力した場合は最大10件とし、宣言順で保持する。位置 X/XB の範囲は 0 より大きく1以下で、端 B を表す1を含む。指数の綴りが違っても同じ十進値は重複として拒否する。SO は YES、YESA、NO を区別する。YES に続く座標行は空欄を含めて保持し、YESA と NO は座標行を消費しない。継続行がある PBEAM では、直後に SO があるときだけ端 A の座標行を省ける。SO がなく末尾属性を指定するなら、先に端 A の座標行を置く。回復点を指定しない場合は、継続行を明示してデータ欄を空ける。[MSC Nastran 2022.1 Quick Reference Guide、pp.2694–2696](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

## せん断中心を基準とする位置を別々に保つ

末尾のせん断係数 K1・K2、テーパによるせん断補正 S1・S2、端 A・B の非構造質量の慣性と反り係数をそれぞれ保持する。非構造質量の重心 M と中立軸 N は、せん断中心からの局所 y・z 座標として区別する。断面テンソルの成立、欠けた断面値の補完・補間、反りと応力回復の適格性、物理単位と K6・質量への変換は後段で扱う。[MSC Nastran 2022.1 Quick Reference Guide、pp.2695–2698](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

## 検証

同じデッキから CBEAM の参照先と PBEAM の断面をたどり、元バイトと宣言位置を照合した。small・large・free の継続、端 B を先に置いた場合の宣言順、10件の境界、指数表記の違う同値、応力回復行の空欄と省略、末尾の属性を確認した。材料の未宣言と種類違い、元の識別情報やデッキが異なる組合せ、確保量・作業量の不足、取消も試験した。端 A の J には負値、絶対値の小さい負値、空欄、正負の零を入力した。後続の非構造質量欄も不正な場合に J の欄が診断されることを確認した。

Windows の MSVC 19.51.36260.0 Debug と Linux の GCC 14.2.0 Debug（ASan・UBSan）で、Nastran の対象8件がそれぞれ通過した。PBEAM は各7ケース・276アサーション。Linux の clang-tidy 21.1.1 で本体の通常構成・確保計測構成と試験を検査した。
