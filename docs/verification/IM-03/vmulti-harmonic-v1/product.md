# 型付き Case の二翼で有限直線後流の調和応答を照合する

上反角0°・5°と各翼 M2/N8・M4/N16・M8/N32 の六条件を、Ubuntu x86-64 の GCC 14.2.0 Release で照合した。主翼の上下動とピッチ、尾翼の上下動とピッチの四駆動を、各条件で一つの作用と因数分解に渡した。射影前の全係数行、全循環、翼別の力とモーメントは、同じ有限コアの独立参照に対して相対 $10^{-9}$ 以内だった。四駆動すべてで元行の相対 $L_2$ 残差は $10^{-11}$ 未満だった。 [独立問題の規約](../../../../rd/aero/multi-surface/records/multiple-surface-harmonic-protocol.md)・[製品結果](results.json)

## Case の元後縁で独立した二翼の後流を作る

型付き Case の静的基準配置から、Case の前後縁対応を通して元面と渦セルを発行する。有限後流の各帯は、同じ結果が発行した正準後縁から閉端まで延ばす。Case の規定する一定移流を保持時間だけ積分した変位を、各端の保持座標に厳密に加える。翼根の正準頂点は各翼の中で共有する。二翼間に接続辺は作らない。 [Case 前後縁対応](../../../specs/aero/wing-edge-pairing-contract.md)・[大域幾何と後縁転写](../../../specs/aero/global-geometry-te-contract.md)

内部の検証入口 `solve_case_wing_harmonic` は、同じ Case の基準配置、密度、一定の既知流れと移流、後流保持距離・時間を使う。翼別の集計範囲と原点は、Case の力・モーメント観測の組から決める。周波数と四つの複素右辺は検証側が明示する。Case に調和解析タグは追加していない。

基準流速は20 m/s、密度は1.2 kg/m³、コア半径は $10^{-5}$ m、閉端は128 m、短縮周波数は0.2。時間依存は $\operatorname{Re}(\Gamma e^{i\omega t})$ とする。各翼の上下動の右辺は $i\omega n_z$。ピッチの右辺は $-[20+i\omega(x-x_a)]n_z$ で、軸の $x_a$ は主翼0.375 m、尾翼5.25 m。駆動していない翼の右辺は厳密にゼロとした。

## 射影前の全係数行と元頂点の荷重を比較する

独立参照は、翼の物理座標から閉輪郭を直接積分し、製品の係数、後縁転写、荷重を呼ばずに算出した。現在の正のコアを両方式で同じ値にそろえ、製品は基準の年齢区間を二分した Gauss16、参照は四分した Gauss48 を用いた。係数と荷重の比較には、製品の中心と包含半径を含めた。

各係数行は独立参照の当該行の $L_2$ ノルムで正規化する。循環と元行残差は駆動ごとに独立循環または右辺の全行 $L_2$ ノルムで割る。翼別荷重は同じ翼・駆動における独立な力またはモーメントの3成分ノルムを使う。ゼロ成分にもこれらの尺度を用いるため、成分ごとのゼロ除算は生じない。

| 上反角 | 各翼の格子 | 係数行 | 循環 | 元行残差 | 翼別の力 | 翼別のモーメント |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0° | M2/N8 | 2.21e−14 | 6.72e−15 | 1.34e−13 | 8.78e−16 | 3.75e−15 |
| 0° | M4/N16 | 4.61e−14 | 8.44e−15 | 2.25e−13 | 2.23e−15 | 1.52e−14 |
| 0° | M8/N32 | 1.95e−13 | 2.05e−14 | 4.99e−13 | 3.12e−15 | 1.73e−14 |
| 5° | M2/N8 | 3.44e−14 | 9.98e−15 | 1.70e−13 | 3.53e−14 | 3.57e−14 |
| 5° | M4/N16 | 4.68e−14 | 1.08e−14 | 2.68e−13 | 1.91e−15 | 1.48e−14 |
| 5° | M8/N32 | 1.95e−13 | 2.05e−14 | 5.61e−13 | 3.95e−15 | 1.56e−14 |

各欄は六条件の行、駆動、翼について該当する最大値を上に丸めた値。各駆動・各翼の値は [結果JSON](results.json) にある。元面ごとの有向辺寄与と三角形面積項を、元頂点の複素力として保存した。各元面の後縁と対角線の辺荷重が除外されることを確認し、元頂点の仕事を独立な翼別の合力・モーメントと照合した。主翼の両半翼を傾けた対照では、共有根に生じる非零の元面寄与が逆符号で打ち消し合った。

相手翼の誘導ブロックを省いた対照の係数差は、全行列ノルム比で最小でも $7.86\times10^{-3}$ だった。別のコアを宣言した Case、誤った元面順、異なる翼の観測を混ぜた組、異なる閉端を持つ準備結果を拒否した。 [製品試験](../../../../tests/coupling/case_wing_geometry_test.cpp)

## 再現

六件の Case は [前後縁対応の検証](../wing-edge-pairing/product.md) と共通のバイト列を読む。元の理想核参照は `rd/evidence-index.json` に記録した `uvlas-rd-evidence-2026-11-part02.tar` から取り出す。生成器は tar 全体、引用元ファイル、六件の個別参照の SHA-256 を照合する。 [独立生成器](reference.py)・[凍結結果](../../../../rd/aero/multi-surface/records/multiple-surface-harmonic-refined-result.json)

| 証拠 | SHA-256 |
| --- | --- |
| 凍結した Gauss48 理想核結果の索引 | `72427831e87dd0a25f98ebccd2d0e3dac6f3bcbeb3ffa4c7181f6ba5a5b537bf` |
| 元参照の tar | `9404de04c6d33bb2195fcf4b82582c16864907ac3d5fd5031b9ef1b02ee6639f` |
| 有限コア独立生成器 | `86a5226879f6a39e93dbc285d1ece4ba80023b37ad453f120ea0061f62f6c332` |
| 有限コア独立参照の全六条件 JSON | `a0d44fd777be6a42909bdad846b068e0766930fff523dec95adf0eb4e211346f` |
| 既定試験用 M2/N8・0° JSON | `5f9b8603a48de1548c8a2e56145cbf5f51ca80772f9f2dd1d1808f080a040790` |

[ビルド手順](../../../build.md) に従い、依存ライブラリとコンパイラを用意する。独立生成器には Python 3.12 と NumPy 2.3.5 を使う。既定試験の1件は、元の tar から別の出力先に再生成し、チェックインしたバイト列と一致した。リポジトリのルートで実行する。

```bash
mkdir -p out/vmulti-harmonic
gh release download rd-evidence-2026-11 \
  --pattern uvlas-rd-evidence-2026-11-part02.tar --dir out/vmulti-harmonic
OPENBLAS_NUM_THREADS=1 python3 -B docs/verification/IM-03/vmulti-harmonic-v1/reference.py \
  --root . --bundle out/vmulti-harmonic/uvlas-rd-evidence-2026-11-part02.tar \
  --output out/vmulti-harmonic/all-six.json --fixture
cmake --preset linux-release -DUVLAS_VMULTI_FORMAL_FIXTURE=out/vmulti-harmonic/all-six.json
cmake --build --preset linux-release --target uvlas_test_coupling_case_wing_geometry uvlas_test_coupling_finite_wing_steady_fw_v01
ctest --preset linux-release --output-on-failure -V \
  -R '^(coupling[.](case_wing_geometry|vmulti_harmonic[.].*|finite_wing_steady.fw_v01))$'
```

対象は有限直線後流の幾何、元係数、循環と元荷重の照合。理想核に対する物理誤差の配分と自由後流の判定は、親パッケージ IM03-0 の受入で扱う。
