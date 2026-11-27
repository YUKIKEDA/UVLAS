# 凍結V-MATの曲面材料点

凍結V-MATの変厚曲面と3層を使う。原本の係数状態を、元実験の基底で物理点の変形勾配へ写した。再現時の参照係数の差は0、内力の相対差は $4.79\times10^{-15}$、質量行列は $2.25\times10^{-16}$ だった。比較の許容差は元実験の $10^{-8}$ を使った。 [凍結結果](../../../../rd/structure/material/records/material-shell-integration-v2-result.json)

3層の内点と2つの層界面の両側を合わせて7点とした。原本の材料入力から工学成分の3D剛性を作り、参照の層面接線から材料軸を作った。Greenひずみと初期ひずみから、世界座標の第二・第一Piola–Kirchhoff応力、Cauchy応力、参照体積当たりのエネルギーを独立に計算した。 [抽出コード](oracle.py) · [製品試験fixture](../../../../tests/structure/fixtures/pg-frozen-material-probes-v1.json)

製品試験は、同じ4次曲面と変厚を有理数で構築し、7点の参照位置を $10^{-12}\,\mathrm{m}$ 以内で照合する。各層の材料則で、材料軸・応力・エネルギーの包絡全体が独立値の $10^{-10}\max(1,|y|)$ 以内に収まることを調べる。材料点試験の変形勾配には凍結係数から写した値を与える。 [製品試験](../../../../tests/structure/pg_frozen_shell_material_test.cpp)

## 変位係数の写像

凍結時の幾何質量による正規化を元の精度で再計算し、アーカイブの変位を元の有理多項式基底に戻す。変換後に有限精度で固定した係数を使い、7点の現位置と変形勾配の差がそれぞれ $10^{-12}\,\mathrm{m}$、$10^{-10}$ 以内であることを調べる。

凍結時と製品で共通のセル局所Legendre係数を使う。171本の製品の選択基底について有理数で求解し、4セルの $R$、$P$、$C$ の全1,096行を3つのCartesian成分で再構成する。各残差は有理数で0だった。出力には選択基底での変位と、凍結内力・接線・質量の同じ変位に沿った縮約値を含める。 [係数の写像と照合](transfer.py)

変位写像を単独で検証した時点では、7点での最大絶対差は現位置が $1.11\times10^{-16}\,\mathrm{m}$、変形勾配が $4.44\times10^{-16}$ だった。GCC 14 Releaseで出力した係数をWindowsの参照計算で2回処理し、同じ出力SHA-256 `565c4d2dd0b29cc7fe25c00bd867fd4f2b3dc5a73f5eba13769d4a1dc109c095` を得た。

## 条件付け済み状態と有限則の組立

製品の幾何Gramを正の積則で作り、認定した変換 $T$ に対して $T\delta q=a$ を解く。ここで $a$ は選択基底に写した凍結変位。解の各成分を小数部256ビットの二進有理数で保持し、各式の残差を $10^{-40}\,\mathrm{m}$ 未満で照合する。$q_0+\delta q$ を製品の状態評価に渡し、7点の現位置、変形勾配、材料軸、エネルギー密度と3種類の応力を照合する。 [製品試験](../../../../tests/structure/pg_frozen_shell_material_test.cpp)

試行状態を、三角形の重心と層中央を組み合わせた24点の有限則（各層8点）で積分する。参照側ではこの有限則の座標と重みから、凍結基底の変位・勾配と独立テンソル式で各積分を計算する。製品側では同じ試行からエネルギー、質量、内力仕事 $\delta q^T f$、接線 $\delta q^T K\delta q$、整合質量 $\delta q^T M\delta q$ を求め、包絡全体がそれぞれの参照値の相対 $10^{-8}$ 以内に収まることを調べる。凍結時の高次則への求積収束、平衡、固有振動はIM-02の残りの検証に含まれる。 [組立期待値の生成](transfer.py) · [変位と組立期待値](../../../../tests/structure/fixtures/pg-frozen-material-transfer-v1.json)

### 組立の作業量

当初の組立試験では、弾性作用で作業量が499,992,912に達した時点で次の9,702を確保できず、上限500,000,000の `ResourceExceeded` が返った。この時点の累積所有確保は894,160,224 Bだった。組立の外向き格子丸めは整数除算と積をlimb数で課金する。全24点・171列・Cartesian 3成分の内力と接線について、両端点の演算を含む上限として組立試験の枠に $2\times10^9$ を指定した。組立に用いる整数値と作業の幅は16,384 bitと65,536 bit、有限和の格子は小数部112 bit、包絡幅の上限は $10^{-21}$ とした。

GCC 14 Releaseで凍結試験2件の2,551アサーションと既存の材料作用試験が通った。成功した24点の組立で、独立のResourcePlanに計上された所有領域と作業量は以下のとおり。

| 作用 | 累積所有確保 [B] | 作業量 | 分割訪問 |
| --- | ---: | ---: | ---: |
| 弾性と接線 | 895,667,872 | 511,941,279 | 678,186 |
| 整合質量 | 812,779,024 | 27,227,928 | 669,981 |

### 質量作用で使う条件付け済みの形状

整合質量の体積行は、認定した基底と変換を用いて $N=B T$ と体積Jacobianを求める。質量作用の入力はこの2つと層密度、方向係数である。24点のうち3層を代表する3点で、従来の参照行が返す $g$、Jacobian、171列の $N$ を有理係数と根係数の両方で厳密に照合した。既存の材料作用と凍結試験2件（3,601アサーション）はGCC 14 Releaseで合格した。

同じ24点の質量作用で、累積所有確保は812,779,024 Bから255,333,520 Bに減った。作業量は27,227,928から17,761,728、分割訪問は669,981から665,829になった。力・接線の組立には物理勾配が必要なので、全参照行を使う。

## 再現

原本NPZを `tools.rd.evidence` のキャッシュに置く。抽出コードは原本、凍結結果、結果が引用する生成ソースのSHA-256を確認し、新しい出力ファイルを作る。

```text
uv run --python 3.12 --no-project --with-requirements rd/requirements-reference.txt --with python-flint==0.8.0 python -X utf8 docs/verification/IM-02/pg-frozen-material-probes/oracle.py out/pg-frozen-material-probes-replay.json
ctest --test-dir <build> --output-on-failure -R "^structure[.]pg_frozen_shell_material$"
<build>/tests/uvlas_test_structure_pg_frozen_shell_material "[.fixture]" > out/pg-frozen-material-basis.txt
uv run --python 3.12 --no-project --with-requirements rd/requirements-reference.txt --with python-flint==0.8.0 python -B -X utf8 docs/verification/IM-02/pg-frozen-material-probes/transfer.py out/pg-frozen-material-basis.txt out/pg-frozen-material-state.json
```

同じWindows環境で独立に2回抽出したfixtureのSHA-256は `83ea890fbc83d651d63c51477c197b5e2ed8c167b08e6714758fa3715d51df97` で一致した。
最初の3D材料点の証拠では、GCC 14 Releaseで577アサーションが合格した。
変位写像と中点則を含む新しいfixtureは再生成で一致し、SHA-256は `37214f9ea0c70772712cd82050dbf96901d412524ac2c7e401e67b53fcc1814b` だった。
