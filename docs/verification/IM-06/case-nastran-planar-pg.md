# IM-06 Case と Nastran の平面・厚さ・材料軸をまとめて照合する

`prepare_case_nastran_planar_pg` は SI 変換済み Case、対応表で認証された外部面、同じデッキの型付き物性欄を受け取る。既存の平面照合で Case の所有 PG を作り、元面、厚さ、積層と材料の照合を通す。続けて同じ参照と厚さ証拠から材料軸を照合し、すべてが成立したときに所有 PG・平面の照合結果・材料軸を返す。既存の共通平面入口には従来の呼出し契約を適用する。[材料軸の厳密照合](nastran-pg-material-angle.md)

## 失敗した段階と出典を保持する

実デッキの正常例は、同じ PG の識別と対応表ダイジェストを平面・軸の両結果に保持し、二層の元角を Case 層順で返した。第2層の Case 軸だけを変えると、合流入口は元要素・PID・第2層・元の要素角と ply 角の欄を示して失敗した。その後、同じ Case とデッキで共通平面入口を使うと、厚さと二層の照合に成功した。合流入口を再試行すると同じ第2層の不一致を返した。

材料軸の件数上限ちょうど、一件不足、取消も検査した。入力は読取り専用で使う。成功した割当は、後段で不一致を返した場合も共有 ResourcePlan に累積する。材料軸の適用範囲は、一定全厚が証明された平面と一定の Case 軸である。

## 検証に使った構成

対象ソースは `98cfb03817359f16fd21d538607647794d9dd793`。Windows 11 の MSVC 19.51.36260 Debug と Release、Ubuntu 24.04 の GCC 14.2.0 Debug（ASan・UBSan）で、`io.external_pg_regional_thickness`、`io.external_planar_pg`、`io.nastran_laminate_cards`、`io.reference_input` が各4件合格した。各構成で、実ファイルの対応表試験は1559アサーションだった。変更した翻訳単位の clang-tidy と clang-format 21.1.2 も通過した。
