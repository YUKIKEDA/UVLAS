# 型付き Case から渦格子と後縁転写を準備する

上反角0°・5°と各翼の格子 M2/N8・M4/N16・M8/N32 の六条件を Ubuntu の GCC 14 Release で照合した。渦頂点1,704点と非貫通点1,344点の全位置が対照と一致した。全非貫通点の法線、元セルの順、翼根の正準頂点、224後縁帯の所有と転写符号も確認した。位置と単位法線の比較許容差は各成分でそれぞれ $5\times10^{-14}\,\mathrm m$ と $5\times10^{-14}$。 [採用仕様](../../../specs/aero/wing-edge-pairing-contract.md)・[独立対照](verification.md)

## 型付き Case の出所と照合

六条件の [Case fixture](../../../../tests/aero/fixtures/wing-edge-pairing-v1/) は、対応則の審査で照合した入力をバイト列のまま保管している。六件の SHA-256、凍結R&D結果と元ソースの SHA-256 は [独立検証コード](reference.py) に記載した。凍結結果の SHA-256 は `72427831e87dd0a25f98ebccd2d0e3dac6f3bcbeb3ffa4c7181f6ba5a5b537bf`。独立検証は Case の chart を有理数で評価し、元の物理翼幅を使う翼平面式と全点を厳密に比較する。製品試験は同じ Case を型検証と SI 変換に通し、製品の chart を使わずに翼平面式から各点の期待値を求める。

| 上反角 | 各翼の格子 | 渦頂点 | 非貫通点 | 後縁帯 |
| --- | --- | ---: | ---: | ---: |
| 0°、5°の各条件 | M2/N8 | 60 | 32 | 16 |
| 0°、5°の各条件 | M4/N16 | 180 | 128 | 32 |
| 0°、5°の各条件 | M8/N32 | 612 | 512 | 64 |
| 合計 | 6件 | 1,704 | 1,344 | 224 |

所有結果は `CaseWingGeometry::prepare` から取得する。結果から参照 chart、元面と局所添字、位置・速度・接線の包含半径、非貫通点の `ResolvedQuery`、大域幾何、正準後縁を参照できる。共有点の速度は既存の大域運動学照合に渡し、照合済みの速度も保持する。V-MULTI の全後縁帯で転写係数は $+1$、後流の向きは正準辺の向きと同じだった。向きを反転した別の試験では法線と後流の向きが反転し、転写係数は $+1$ を保った。試験は Case を介した幾何準備と後縁転写までを対象にする。親パッケージの V-MULTI 検証では、二翼の調和求解と荷重および後流全体を扱う。 [製品試験](../../../../tests/coupling/case_wing_geometry_test.cpp)

## 対応の切替と異常入力を同じ入口で確認する

共通翼幅の端が異なる面は暗黙対応で拒否し、明示した有理数の対応で位置と接線が得られた。前縁と後縁を結ぶ線分が交差する対応は拒否した。外周の切り欠き、離れた区分の重複、非貫通点で左右の接線が異なる節点でも準備は失敗した。共有の試験では、元面と共有方向の入替え、離散辺数の不一致、両端と渦頂点だけが一致する異なる曲線を検出した。内部に共有した後縁からは後流を放出しない。指定運動では時刻と片側を固定し、回転後の位置・法線・接線と並進速度を解析値で確認した。同じ位置で速度が異なる共有は拒否した。取消と作業量・容量の上限に達した後も、以前の準備結果を参照できた。

指定運動の共有は、静止を証明できる場合と、準備済み Frame および厳密な運動係数の一致を確認できる場合に認定する。異なる式・パラメータ化を同一の軌道として認定する追加の証明は、この入口には含まれない。

Ubuntu 24.04.5 の GCC 14.2.0 Release で `io.case_schema`、`coupling.case_motion_time_series`、`coupling.case_wing_geometry` の3件が通過した。幾何試験は8ケース、68,496アサーション。独立対照の通常実行と `-O` の出力は一致し、SHA-256 は `f9408673ee6db03a5f661cfd08b2de0f02df137ce9fe86a5960f743aa6d3bdbc` だった。

## 再現

Ubuntu で [ビルド手順](../../../build.md) に従い GCC 14 と依存ライブラリを用意し、リポジトリのルートで実行する。

```bash
python3 -B docs/verification/IM-03/wing-edge-pairing/reference.py \
  --cases tests/aero/fixtures/wing-edge-pairing-v1
cmake --preset linux-release
cmake --build --preset linux-release --target uvlas_test_coupling_case_wing_geometry uvlas_test_io_case_schema uvlas_test_coupling_case_motion_time_series
ctest --preset linux-release --output-on-failure \
  -R '^(io[.]case_schema|coupling[.]case_(wing_geometry|motion_time_series))$'
```
