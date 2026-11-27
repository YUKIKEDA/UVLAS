# 凍結V-MATの曲面材料点

凍結V-MATの変厚曲面と3層を使う。原本の係数状態を、元実験の基底で物理点の変形勾配へ写した。再現時の参照係数の差は0、内力の相対差は $4.79\times10^{-15}$、質量行列は $2.25\times10^{-16}$ だった。比較の許容差は元実験の $10^{-8}$ を使った。 [凍結結果](../../../../rd/structure/material/records/material-shell-integration-v2-result.json)

3層の内点と2つの層界面の両側を合わせて7点とした。原本の材料入力から工学成分の3D剛性を作り、参照の層面接線から材料軸を作った。Greenひずみと初期ひずみから、世界座標の第二・第一Piola–Kirchhoff応力、Cauchy応力、参照体積当たりのエネルギーを独立に計算した。 [抽出コード](oracle.py) · [製品試験fixture](../../../../tests/structure/fixtures/pg-frozen-material-probes-v1.json)

製品試験は、同じ4次曲面と変厚を有理数で構築し、7点の参照位置を $10^{-12}\,\mathrm{m}$ 以内で照合する。各層の材料則で、材料軸・応力・エネルギーの包絡全体が独立値の $10^{-10}\max(1,|y|)$ 以内に収まることを調べる。製品側の変形勾配には凍結係数から写した値を与える。製品基底への係数変換、組立済み接線・質量、平衡・固有振動の凍結照合は別の検証を要する。 [製品試験](../../../../tests/structure/pg_frozen_shell_material_test.cpp)

## 再現

原本NPZを `tools.rd.evidence` のキャッシュに置く。抽出コードは原本、凍結結果、結果が引用する生成ソースのSHA-256を確認し、新しい出力ファイルを作る。

```text
uv run --python 3.12 --no-project --with-requirements rd/requirements-reference.txt --with python-flint==0.8.0 python -X utf8 docs/verification/IM-02/pg-frozen-material-probes/oracle.py out/pg-frozen-material-probes-replay.json
ctest --test-dir <build> --output-on-failure -R "^structure[.]pg_frozen_shell_material$"
```

同じWindows環境で独立に2回抽出したfixtureのSHA-256は `83ea890fbc83d651d63c51477c197b5e2ed8c167b08e6714758fa3715d51df97` で一致した。
GCC 14 Releaseでは製品試験の577アサーションがすべて合格した。
