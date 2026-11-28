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

### 弾性組立で共有する合成基底

弾性組立で呼出し側が合成を選ぶと、$n=1$ の認定済み基底 $B$ とその変換 $T$ をセル多項式の段階で合成する。既定の評価では、従来と同じ作業枠で作用を計算できる。$R$、$P$、$C$ の群ごとに分母をそろえ、合成の正負を保って整数で積和する。キャッシュは組立中だけ元の認定結果を借り、別の認定結果を渡したときは世代不一致を返す。

準備では保有する $B$ と $T$ の整数を読み、係数幅、更新回数、保存配列の確保枠を先に照合する。係数は各群512 bit以内で保持する。点のセル局所座標の共通分母が32 bit以内なら、値と面内一次微分を1,536 bitの整数で集約する。積和の幅は係数512 bit、総次数21に対する $21\cdot32$ bit、加算項と微分倍率を合わせて1,202 bit以下となる。点が範囲外なら従来の参照評価へ戻り、照合で消費した作業も記録する。

凍結171列の実係数を独立に合成すると、非零項数は34,050から160,528へ増えた。合成後の共通分母と分子の幅は最大285 bit、313 bitだった。4セルの各16点で、各列の $R$、$P$、$C$ とその面内微分を独立な有理積和と照合した98,496成分は、零を含めてすべて一致した。

凍結24点の同じ試行について、直接経路と合成基底経路のエネルギー、内力全列、接線作用全列の有理な区間端点が一致した。
Ubuntu 24.04.5、GCC 14.2.0 Release で37,318アサーションを確認した。この比較には参照行と合成基底の所有検査も含まれる。

| 弾性作用の経路 | 累積所有確保 [B] | 作業量 | 分割訪問 |
| --- | ---: | ---: | ---: |
| 直接 | 895,667,872 | 511,941,279 | 678,186 |
| 合成基底 | 189,631,072 | 678,514,358 | 125,769 |

## 凍結した高次則との比較

凍結時の体積積分は、4背景セルの各軸45点と、3層の深さ各7点を合わせた170,100点を使う。凍結変位と同じ材料則でその点列を独立に再計算し、アーカイブに保存した内力仕事、接線と質量の縮約値を相対 $10^{-8}$ で照合した。表はその再計算を基準とする。面内点数は各軸、深さ点数は各層に対する数である。 [凍結時の求積と組立](../../../../rd/structure/material/experiments/check_material_shell_integration_v2.py) · [収束の再計算](convergence.py)

| 面内 × 深さ | 全点数 | 弾性3観測の最大相対差 | 総質量の相対差 | 質量仕事の相対差 |
| --- | ---: | ---: | ---: | ---: |
| $6^2\times3$ | 1,296 | $1.08\times10^{-4}$ | $3.95\times10^{-16}$ | $1.62\times10^{-6}$ |
| $9^2\times1$ | 972 | $6.04\times10^{-6}$ | $1.71\times10^{-6}$ | $6.20\times10^{-5}$ |
| $9^2\times2$ | 1,944 | $4.97\times10^{-6}$ | $1.97\times10^{-16}$ | $3.15\times10^{-10}$ |
| $16^2\times2$ | 6,144 | $6.72\times10^{-8}$ | $0$ | $1.56\times10^{-10}$ |
| $9\times17\times2$ | 3,672 | $2.34\times10^{-9}$ | $3.95\times10^{-16}$ | $8.78\times10^{-11}$ |
| $17^2\times2$ | 6,936 | $2.46\times10^{-10}$ | $3.95\times10^{-16}$ | $5.37\times10^{-11}$ |
| $24^2\times3$ | 20,736 | $2.50\times10^{-12}$ | $3.95\times10^{-16}$ | $5.44\times10^{-16}$ |

24点の弾性組立から単純に点数比例で見積もると、170,100点で累積所有確保約6.35 TB、同じ許容差を満たした面内9×17、深さ2の3,672点でも約137 GBとなる。これは累積確保の目安であり、プロセスのピークや未計測の実行時間は含めない。48 GiBのプランでの弾性3観測を確認するには、同じ凍結変位、求積精度、作用を保つ実装の追加検証が必要である。

1,944点の質量候補は座標を小数部28 bit、重みを44 bitに丸めた。丸め後の総質量の独立差は $1.04\times10^{-11}$、変位に沿う質量仕事は $1.64\times10^{-10}$ だった。試験には背景セル内の各点を所有三角形へ写し、層幅と三角形の座標変換を重みから除いた局所則を渡す。 [高次質量試験](run-high-mass.cpp.txt)

GCC 14 Releaseで、この1,944点を製品の整合質量入口へ渡した。総質量の包絡と凍結変位に沿う質量仕事の包絡は、どちらもそれぞれの独立参照値の相対 $10^{-8}$ 以内に収まり、8,995アサーションが合格した。作用に要した時間は731.424秒、累積所有確保は19,402,002,064 B、作業量は4,279,433,365、分割訪問は53,931,909だった。計測は作用だけを対象とし、48 GiBの独立のResourcePlanを使用した。

## 再現

原本NPZを `tools.rd.evidence` のキャッシュに置く。抽出コードは原本、凍結結果、結果が引用する生成ソースのSHA-256を確認し、新しい出力ファイルを作る。

```text
uv run --python 3.12 --no-project --with-requirements rd/requirements-reference.txt --with python-flint==0.8.0 python -X utf8 docs/verification/IM-02/pg-frozen-material-probes/oracle.py out/pg-frozen-material-probes-replay.json
ctest --test-dir <build> --output-on-failure -R "^structure[.]pg_frozen_shell_material$"
ctest --test-dir <build> --output-on-failure -R "^structure[.]pg_frozen_shell_material[.]composed$"
<build>/tests/uvlas_test_structure_pg_frozen_shell_material "[.fixture]" > out/pg-frozen-material-basis.txt
uv run --python 3.12 --no-project --with-requirements rd/requirements-reference.txt --with python-flint==0.8.0 python -B -X utf8 docs/verification/IM-02/pg-frozen-material-probes/transfer.py out/pg-frozen-material-basis.txt out/pg-frozen-material-state.json
uv run --python 3.12 --no-project --with-requirements rd/requirements-reference.txt --with python-flint==0.8.0 python -B -X utf8 docs/verification/IM-02/pg-frozen-material-probes/convergence.py out/pg-frozen-material-basis.txt --mass-rule-output out/vmat-high-rule.json
```

3,672点の弾性試験では、面内と厚さ方向の座標を小数部28 bit、重みを44 bitに丸める。独立再計算で、丸め後のエネルギー、内力仕事、接線仕事の凍結則からの相対差は、それぞれ $2.33\times10^{-9}$、$2.27\times10^{-9}$、$2.15\times10^{-9}$ だった。製品試験は3層各1,224点を読み、3観測の包絡全体をそれぞれの凍結則の値の相対 $10^{-8}$ で検査する。組立には48 GiBの累積所有確保枠を使う。
製品試験に渡した点列の SHA-256 は `8e4054f426c3a0fa6810f64fb7a3462d6b1a34efad602a4aadbf7889242336e7` である。

Ubuntu 24.04.5 の GCC 14.2.0 Release で、固定幅の有理数を外向きに丸める試験を実行した。分母が257 bitと9,001 bitの入力について、正負それぞれの上下端点と、作業幅を1 bit下げた場合の拒否を調べた。14アサーションが通過した。
同じ構成の `structure.pg_shell_material`、`structure.pg_frozen_shell_material`、`structure.pg_discrete_gram`、`coupling.case_pg_shell_material` も通過した。凍結24点の比較と Case の材料作用への接続を含む。
`structure.pg_prepared_patch`、`io.external_shell_cross_format`、`io.external_case_chart` では、異なる層・パッチの引渡しと資源不足の拒否、認証済み外部参照から作った Case の作用を照合した。対象3件が通過した。

正式試験の初回は、スカラー演算の上限を5,000億に設定し、8,048.21秒で `ResourceExceeded` を返した。
失敗時の累積所有確保は25,122,895,240 B、プロセスの最大常駐量は435,176 KiBだった。
同じ点列の先頭・層境界・末尾から選んだ14点では、1点当たりのスカラー演算課金が145,742,116〜156,778,567だった。
全点の所属を検査した上で同じ14点を使い、演算枠だけを変えた対照では、21億の枠で `ResourceExceeded` を再現した。
25億の枠では、同じ14点の組立と内部の区間幅検査を通過し、2,198,115,510演算を課金した。
入力と区間幅の条件をそろえた対照から、初回の全点実行は演算枠の不足で止まったと判断した。
観測した最大課金で全3,672点と準備の分を見積もると約5,757億となるため、再試験のスカラー演算上限を1兆に設定した。
点列、丸め精度、判定値、48 GiBの累積所有確保枠は維持する。

```bash
uv run --python 3.12 --no-project --with-requirements rd/requirements-reference.txt --with python-flint==0.8.0 python -B -X utf8 docs/verification/IM-02/pg-frozen-material-probes/convergence.py out/pg-frozen-material-basis.txt --elastic-rule-output out/vmat-elastic-rule.json
cmake --preset linux-release -DUVLAS_FORMAL_VMAT_ELASTIC_FIXTURE="$PWD/out/vmat-elastic-rule.json"
cmake --build --preset linux-release --target uvlas_test_structure_pg_frozen_shell_material
ctest --preset linux-release -R '^structure[.]pg_frozen_shell_material[.]high_elastic$'
```

GCCで1,944点を再実行する場合は、Releaseで構築した3つの静的ライブラリとCatch2をリンクする。以下では `BUILD` にビルドディレクトリ、`EXACT_DEPS_INCLUDE` にBoostとnlohmann/jsonの依存ヘッダー、`CATCH_PREFIX` にCatch2のインストール先を指定する。

```bash
g++-14 -std=c++23 -O2 -pthread -Iinclude -Isrc/structure -Itests/structure \
  -isystem "$EXACT_DEPS_INCLUDE" -isystem "$CATCH_PREFIX/include" \
  "-DUVLAS_PG_FROZEN_MATERIAL_PROBES=\"$PWD/tests/structure/fixtures/pg-frozen-material-probes-v1.json\"" \
  "-DUVLAS_PG_FROZEN_MATERIAL_TRANSFER=\"$PWD/tests/structure/fixtures/pg-frozen-material-transfer-v1.json\"" \
  -x c++ docs/verification/IM-02/pg-frozen-material-probes/run-high-mass.cpp.txt -x none \
  -Wl,--start-group "$BUILD/libuvlas_structure.a" "$BUILD/libuvlas_numerics.a" \
  "$BUILD/libuvlas_core.a" -Wl,--end-group \
  "$CATCH_PREFIX/lib/manual-link/libCatch2Main.a" "$CATCH_PREFIX/lib/libCatch2.a" \
  -o out/vmat-high
(cd out && ./vmat-high '[.high-mass]' --rng-seed 34324)
```

有理点列は再生成で一致した。SHA-256は `bb5ea724af885cbba4bd01212697f9f4cd7101f411c903b4995104043e05d39c` だった。

同じWindows環境で独立に2回抽出したfixtureのSHA-256は `83ea890fbc83d651d63c51477c197b5e2ed8c167b08e6714758fa3715d51df97` で一致した。
最初の3D材料点の証拠では、GCC 14 Releaseで577アサーションが合格した。
変位写像と中点則を含む新しいfixtureは再生成で一致し、SHA-256は `37214f9ea0c70772712cd82050dbf96901d412524ac2c7e401e67b53fcc1814b` だった。
