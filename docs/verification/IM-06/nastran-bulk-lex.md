# MSC Nastran の欄と継続行を元の位置付きで読む

取得済みの `nastran-deck` スナップショットから、MSC Nastran 2022.1 の Bulk Data のカード名とデータ欄を読む。成功したときにだけスナップショットの所有を移し、元バイトの SHA-256、宣言元、行・桁・バイト位置を保持する。欄の内容は元の文字列を参照する。読み直しやファイル全体の複製は発生しない。[IO4-03](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran)

## 固定幅と自由形式で空欄を保つ

small format は1行から8桁幅のデータ欄を8個取り出す。large format は1行ごとに16桁幅の欄を4個取り出し、識別子で結ばれた二つの半行を一組とする。名前付き継続は物理行の順序に依存せず、親側の第10欄の先頭文字を除いた識別子と、継続側の第1欄の識別子を照合する。small format では、親の末尾欄と次のデータ行の先頭欄が空なら自動継続としてつなぐ。

comma 区切りの free format では、前の行が8欄の区切りまで埋まり、次の行が comma で始まれば自動継続としてつなぐ。80桁未満なら、1行内に8個を超えるデータ欄を保持する。連続した comma の間は空欄として記録する。固定幅の空欄も省かない。値がある欄では綴りを保存し、桁数、指数、単位、既定値の解釈は後段へ渡す。[MSC Nastran 2022.1 Linear Static Analysis User’s Guide、pp.51–54](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Linear_Static_Analysis_User_Guide/raw/resource/enus/MSC_Nastran_2022.1_Linear_Static_Analysis_User_Guide.pdf) [Getting Started Guide、pp.83–85](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2022.1_Getting_Started_Guide/raw/resource/enus/MSC_Nastran_2022.1_Getting_Started_Guide.pdf)

## 継続の不整合と資源不足を区別する

継続に使う識別子の重複、親の欠落、循環は失敗にする。物理行が順不同でも、カードの先頭行の順と各欄の位置を保持する。試験では三つの形式、CRLF、途中のコメント、順不同の継続、空欄、元の指数表記を照合した。バイト、行、カード、欄、作業量、確保サイズの上限と取消も検査した。

現在の範囲は Bulk Data 部分の字句構造である。固定幅の混在継続、名前付き free 継続、8欄の途中から次の物理行へ続く free 形式、等号を区切りに使う free 形式、継続か注記か確定できない末尾欄は `Unsupported` となる。INCLUDE、Executive／Case Control、カードの物理変換、sidecar、単位は、各処理の契約に従って後段で扱う。

Windows の MSVC 19.51.36260.0 Debug と Linux の GCC 14.2.0 Debug（ASan・UBSan）で、Nastran と既存の Gmsh の対象2件がそれぞれ通過した。Nastran は各4ケース・300アサーション。Linux の clang-tidy 21.1.1 は本体の通常構成・確保計測構成と回帰試験を対象とした。
