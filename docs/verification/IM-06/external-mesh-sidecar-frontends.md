# CLI と Python で版付き対応表を検査

<!-- verification: FD-03 IO4-01 IO4-03 IO4-04 R-05 -->

共通の Case 読取りが、認証済み Gmsh・Nastran 対応表を版付き JSON として検査する。CLI では両形式の不正な対応表を拒否し、正しい形式を受理した。Python の `load_case` と `case_from_dict` では同じ Nastran 対応表を読み、成功と InvalidInput が CLI の結果に対応することを確認した。認証後の元面割当と物性照合は[実ファイルによる検証](external-mesh-sidecar-role.md)に記録した。

元面の一時領域は、借用参照と面の並びだけを計上済み領域に構築する。MSVC Debug でデバッグ用イテレータが計上済み領域を超えていた構成でも、元面と Case の両方を保持し、同じ厚さ・材料照合が成立した。

統合コミット `cbe2739b299ce4c63f8aac44515818569feb4b3d` で Windows 11 の MSVC Debug は `io.reference_input`、`io.external_planar_pg`、`cli.validate` の3件、MSVC Release と Python 3.12 は `cli.validate`、`python.case_file_operations` の2件が通った。Ubuntu 24.04 の GCC 14.2 Release と Python 3.12 ではこれらと親変更のPG回帰を合わせて6件が通った。
