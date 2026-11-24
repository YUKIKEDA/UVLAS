# MSC Nastran の large 欄と固定幅の混在継続

取得した deck の元バイトから、small と large の固定幅を、8個の論理欄が揃ったところで切り替えられる。large は4欄ずつの半行を二つ使う。名前付きの継続は、物理行が前後していても識別子で照合する。 [Linear Static Analysis User’s Guide、pp.52–54](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Linear_Static_Analysis_User_Guide/raw/resource/enus/MSC_Nastran_2022.1_Linear_Static_Analysis_User_Guide.pdf)

large free の `GRID*,` と `*,` は4欄ずつ読み、4欄目の後の区切りを data field と数えない。固定幅で隣接する `*` および末尾の `+` と次行の `*+` も照合する。省略した座標や SEID の欄は空欄で保持する。余分なデータ欄や後半の行が欠けた場合は拒否する。 [MSC Nastran 2022.1 Quick Reference Guide、pp.1118–1119・2079](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

独立試験で論理欄の順序、実数の綴り、元の桁位置を比較した。large free で宣言した GRID の正確な位置を、CQUAD4 の節点参照からたどれることも確認した。名前付き free と free 内の幅切替えは別の入力形式として残る。

Windows の MSVC 19.51 Release と Linux の GCC 14.2 Release は Nastran の対象4件が各通過した。Linux の clang-tidy 21.1.1 は字句読取りと2本の試験翻訳単位を対象とした。Linux の分離コピーで以前の読取りに置き換えると、字句とシェルの対象2件が失敗した。現行版に戻すと4件とも通過した。
