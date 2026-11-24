# Nastran デッキを開き直さずに子参照を検査する

`uvlas validate` は、SHA-256 が宣言と一致した Nastran デッキを末尾まで読む。Bulk Data の継続行を検査し、採用仕様が列挙する17種類のカードだけを含むときに、子参照がないと判定する。`INCLUDE` や範囲外のカードを検出すると、デッキ内のバイト位置とケース上の宣言位置を報告する。型付きカードの解釈と sidecar の検査は IM-06 の別工程で行う。 [参照の開始契約](../../specs/foundation/reference-snapshot-start-contract.md) [採用カード](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran)

CLI の回帰は、子参照のないデッキを通過して sidecar の検査位置へ進むことを確かめる。`INCLUDE` を加え、ダイジェストの不一致が先に返ることを検査する。ダイジェストを合わせると、次のファイルを開く前に `INCLUDE` の位置を報告する。C++ の回帰は採用カード名、コメント、未採用カード、終端後の行、容量と作業枠の境界を確かめる。

Windows の MSVC 19.51 Release と Linux の GCC 14.2 Release で、CLI、Nastran 字句、参照入力、Gmsh メッシュの CTest が各4件通過した。clang-tidy 21.1.1 で変更した5翻訳単位を検査した。GCC の分離コピーで未採用の `CBUSH` を受理するように変えると Nastran の回帰が失敗し、復元後は4件とも通過した。
