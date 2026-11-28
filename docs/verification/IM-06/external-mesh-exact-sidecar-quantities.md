# Gmsh と Nastran の明示全厚を同じ Case で厳密照合する

全厚 2.54 mm と 0.00254 m は厳密に一致した。Gmsh の対応表と Nastran の PSHELL に同じ全厚を指定し、一つの平面 Shell を準備して確かめた。MSH、BDF、両対応表は実ファイルとして Case から取得し、元バイトを個別に SHA-256 で照合した。独立した単位換算は $(254/100)(1/1000)=127/50000$ m となる。両経路が返す厳密値はこの分数と一致し、同じ PG 参照識別子を得た。[対応表の採用契約と最小例](../../specs/foundation/external-mesh-sidecar-contract.md)

Case の最後の厚さ頂点だけを `0.00254000000000000000000000000001` m に変えると差は $10^{-32}$ m になる。元メッシュと対応表は同じファイルを使い、両経路が `PrepareModel / InvalidInput` を返した。診断は元要素 701 と Case の厚さ三角形 1 を示す。Gmsh は対応表の `full_thickness`、Nastran は PSHELL 90 の T 欄も示した。

Gmsh の材料対照では、同じ物性値を持つ別の Case 材料 ID を対応表に指定した。結び付けの失敗は物理グループ、対応表の行、Case Shell を示した。Nastran では MAT1 の E を `200000.00000000001` N/mm² に変え、E 欄のバイト位置を照合した。厚さの差と Nastran の E の差は binary64 に丸めると消える。負対照の後に正常な Case と認証済み参照を再利用すると、両形式で元と同じ PG 参照識別子を得た。

Nastran の対応表だけに存在しない EID 999 を指定すると、認証には成功し、元面の結び付けで拒否された。診断は未解決の EID と対応表の元行を示した。

## 設計の完了条件と証拠

| 完了条件 | 証拠 |
| --- | --- |
| 対応表の形式、元ファイルと Case の識別、タグと元診断 | 実ファイルの認証、型と元タグの照合に加え、異なるダイジェストを拒否する共通取得試験が通過した。[対応表の所有と認証](external-mesh-sidecar-role.md) |
| 元面と Case Shell／PG の対応 | 元節点順、有理チャート、被覆、向き、内部辺と外周を検査する。曲面では独立式による位置差 1/40 m と法線成分を使い、全域の上界と変換予算を検査した。[平面の対応と測度](external-shell-cross-format.md)・[元境界の照合](external-pg-edge-incidence.md)・[曲面の変換予算](external-mesh-conversion-budget.md) |
| 材料と厚さの厳密 SI 照合 | 上記の両形式による一致・拒否と再準備に加え、厚さの違う隣接領域、共通細分、元法線を反転した場合の層順を照合済み。[領域の厚さ](external-pg-regional-thickness.md)・[層順と材料](nastran-pcomp-pg-quantities.md)・[厚さの異なる接続の両側](../IM-01/pg-prepared-patch-trace.md) |
| 近似一致を受け入れる予算と失敗時の保証 | 独立式で得た距離と法線を使い、各予算が不足する対照を拒否した。元領域が一致しないときは適格な変換予算と観測 ID を要求する。資源上限ちょうど、一件不足、取消と元 Case の保持は対応表試験で確認した。[変換予算](external-mesh-conversion-budget.md)・[資源と取消](external-mesh-sidecar-role.md) |

認定の範囲は外部面の対応と明示物性の照合である。変換前後の質量、仕事、変位、固有モードを含む実務規模の物理検証は R-05 の条件に従う。

## Ubuntu Release の再検証

ソースは未コミットの変更がない `50d7c911ebef83a26c4a7d71558d964f900c222a`。Ubuntu 24.04.5、GCC 14.2.0 Release の `linux-x64-gnu14.2-sse2` で、対象 CTest 4 件が通過した。新しい試験は 67 アサーションで、登録済みの形式間試験全体は 3 ケース、1862 アサーションだった。

```bash
cmake --build --preset linux-release --target uvlas_test_io_external_shell_cross_format uvlas_test_io_external_pg_regional_thickness uvlas_test_io_reference_input
ctest --preset linux-release -R '^io[.](external_shell_cross_format|external_case_chart|external_pg_regional_thickness|reference_input)$'
```

隔離した複製で、`external_pg_regional_thickness.cpp` の一定厚さ比較を `!= 0` から `< 0` に変えた。Case の厚さが外部指定を上回るときの拒否が消え、Gmsh と Nastran の両方の拒否アサーションが失敗した。同じ厚さ検査を含む最終版は 67 アサーションが通過した。

## Debug の CTest 登録

既存の `io.external_shell_cross_format` は、変更前の PR #382 の Linux Debug と、この変更を追加した PR #387 の Windows Debug で、60 秒の上限に達した。
厳密比較は `io.external_shell_cross_format.exact_sidecar` として独立に実行する。二つの CTest の上限を Debug で300秒、Release で60秒とする。
分離後の Ubuntu 24.04.5、GCC 14.2.0 Release で、既存2ケースの1,795アサーション、厳密比較の67アサーションを含む対象5件が通過した。
同じ環境の Debug では、既存試験が39.95秒、厳密比較が0.81秒で通過した。アサーション数は Release と一致した。

```bash
ctest --preset linux-release -R '^io[.](external_shell_cross_format([.]exact_sidecar)?|external_case_chart|external_pg_regional_thickness|reference_input)$'
ctest --preset linux-debug -R '^io[.]external_shell_cross_format([.]exact_sidecar)?$'
```
