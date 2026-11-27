# 凍結V-MATの曲面材料点

凍結V-MATの変厚曲面と3層を使う。原本の係数状態を、元実験の基底で物理点の変形勾配へ写した。再現時の参照係数の差は0、内力の相対差は $4.79\times10^{-15}$、質量行列は $2.25\times10^{-16}$ だった。比較の許容差は元実験の $10^{-8}$ を使った。 [凍結結果](../../../../rd/structure/material/records/material-shell-integration-v2-result.json)

3層の内点と2つの層界面の両側を合わせて7点とした。原本の材料入力から工学成分の3D剛性を作り、参照の層面接線から材料軸を作った。Greenひずみと初期ひずみから、世界座標の第二・第一Piola–Kirchhoff応力、Cauchy応力、参照体積当たりのエネルギーを独立に計算した。 [抽出コード](oracle.py) · [製品試験fixture](../../../../tests/structure/fixtures/pg-frozen-material-probes-v1.json)

製品試験は、同じ4次曲面と変厚を有理数で構築し、7点の参照位置を $10^{-12}\,\mathrm{m}$ 以内で照合する。各層の材料則で、材料軸・応力・エネルギーの包絡全体が独立値の $10^{-10}\max(1,|y|)$ 以内に収まることを調べる。材料点試験の変形勾配には凍結係数から写した値を与える。 [製品試験](../../../../tests/structure/pg_frozen_shell_material_test.cpp)

## 変位係数の写像

凍結時の幾何質量による正規化を元の精度で再計算し、アーカイブの変位を元の有理多項式基底に戻す。変換後に有限精度で固定した係数を使い、7点の現位置と変形勾配の差がそれぞれ $10^{-12}\,\mathrm{m}$、$10^{-10}$ 以内であることを調べる。

凍結時と製品で共通のセル局所Legendre係数を使う。171本の製品の選択基底について有理数で求解し、4セルの $R$、$P$、$C$ の全1,096行を3つのCartesian成分で再構成する。各残差は有理数で0だった。出力には選択基底での変位と、凍結内力・接線・質量の同じ変位に沿った縮約値を含める。製品が組み立てる作用、平衡、固有振動の凍結照合はIM-02の残りの検証に含まれる。 [係数の写像と照合](transfer.py)

復元後の7点における最大絶対差は、現位置が $1.11\times10^{-16}\,\mathrm{m}$、変形勾配が $4.44\times10^{-16}$ だった。GCC 14 Releaseで出力した係数とWindowsの再現環境を用いた変位写像の出力SHA-256は `565c4d2dd0b29cc7fe25c00bd867fd4f2b3dc5a73f5eba13769d4a1dc109c095` で、2回の実行で一致した。

## 再現

原本NPZを `tools.rd.evidence` のキャッシュに置く。抽出コードは原本、凍結結果、結果が引用する生成ソースのSHA-256を確認し、新しい出力ファイルを作る。

```text
uv run --python 3.12 --no-project --with-requirements rd/requirements-reference.txt --with python-flint==0.8.0 python -X utf8 docs/verification/IM-02/pg-frozen-material-probes/oracle.py out/pg-frozen-material-probes-replay.json
ctest --test-dir <build> --output-on-failure -R "^structure[.]pg_frozen_shell_material$"
<build>/tests/uvlas_test_structure_pg_frozen_shell_material "[.fixture]" > out/pg-frozen-material-basis.txt
uv run --python 3.12 --no-project --with-requirements rd/requirements-reference.txt --with python-flint==0.8.0 python -B -X utf8 docs/verification/IM-02/pg-frozen-material-probes/transfer.py out/pg-frozen-material-basis.txt out/pg-frozen-material-state.json
```

同じWindows環境で独立に2回抽出したfixtureのSHA-256は `83ea890fbc83d651d63c51477c197b5e2ed8c167b08e6714758fa3715d51df97` で一致した。
GCC 14 Releaseでは製品試験の577アサーションがすべて合格した。
