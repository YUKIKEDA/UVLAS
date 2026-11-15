# 開発用Python拡張

`_uvlas` は、密行列の保持と状態・資格・確定の製品処理を呼ぶ低水準の開発APIです。物理ソルバーとケースファイルの解析入口はまだありません。互換性は固定していません。

```python
import numpy as np
import _uvlas

plan = _uvlas.ResourcePlan.create(4096).value
result = _uvlas.DenseReal.copy(np.asfortranarray([[1.0, 2.0], [3.0, 4.0]]), plan)
if result.ok:
    matrix = result.value
    array = matrix.view()  # 所有者を保持する読取り専用の配列
else:
    print(result.failure.code, result.failure.details)
```

`copy` はnative float64/complex128の2次元Fortran連続配列だけを受け、GIL保持中にコピーします。非連続、別のdtype、非有限値は失敗として返します。上の `asfortranarray` のように、必要な変換は呼出し側が明示してください。返すviewは、`matrix` や `result` の破棄後も有効です。WRITEABLEは再有効化できません。

生成・検査・状態の変更を行う処理は `Result` を返します。`view()`、`definition()`、確定状態の `payload()` は読取り専用のNumPy配列を直接返します。`ok` が偽なら `failure`、真なら `value` を使います。失敗は `FailureCode` の名前で区別でき、`details` に分かっているIDと資源量が入ります。未知情報は `None`。呼出し引数のPython型や構文の間違い、Python/NumPy自身のオブジェクト確保失敗には、Pythonの例外が使われます。

`Model.create` の定義と `TrialState` のpayloadは、1次元の連続uint8配列です。これらは呼出し側が用意するcanonicalなprotocol bytesで、変位や後流をfloat配列へ変換する入口ではありません。モデルと確定状態の配列は不変で、trialの取得は `payload_copy(plan)` による独立snapshotだけです。snapshotにもpayload予算が要ります。

状態の評価では、`trial.target().value` を計算前に取得し、そのtargetと計算した4種類の上界を `EvaluationReport.create` へ渡します。上界の順序は残差・幾何・時間・精度。`None` は未評価です。reportは上界の数値・世代を検査しますが、物理評価を代行しません。古いreportや別候補のreportは資格に使えません。

所有済み配列の検査、候補のhash、資格判定、確定の間はGILを解放します。可変trial、plan、packetへの同時操作はgetterも含めて即座に `InvalidInput` で拒否します。Python callbackは呼びません。通常GILを持つCPython 3.12の単一interpreterだけが対象です。入力を同期コピーしている間に、他のnativeコードで同じ入力メモリを書き換えないでください。

## 構築

CPython 3.12のinterpreterと開発ファイル、NumPy 2.xを用意します。NumPyの検証用固定版は `tests/python/requirements.txt` にあります。vcpkgのpython featureには、通常構成と別の導入prefixを指定してください。

```powershell
uv venv --python 3.12 out/python-env
uv pip install --python out/python-env/Scripts/python.exe -r tests/python/requirements.txt
pwsh tools/dev/msvc.ps1 cmake -S . -B out/build/python -G 'Ninja Multi-Config' -DUVLAS_BUILD_PYTHON=ON -DUVLAS_WARNINGS_AS_ERRORS=ON -DBUILD_TESTING=ON '-DVCPKG_MANIFEST_FEATURES=tests;python' -D "Python3_EXECUTABLE=$((Resolve-Path out/python-env/Scripts/python.exe).Path)" -D "CMAKE_TOOLCHAIN_FILE=$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" -D "VCPKG_INSTALLED_DIR=$PWD/out/vcpkg_python_installed"
pwsh tools/dev/msvc.ps1 cmake --build out/build/python --config Release --target uvlas_python
pwsh tools/dev/msvc.ps1 ctest --test-dir out/build/python -C Release -R '^python.ownership$' --output-on-failure
```

生成された拡張のディレクトリを `PYTHONPATH` に追加してimportします。Linuxでは同じoptionと外部interpreter、専用prefixを指定し、GCC 14/Ninjaで構築します。通常のlintは任意機能のソースを除外するため、この構成のcompile databaseを使って `clang-tidy -p out/build/python src/python/module.cpp` も実行してください。wheel配布はまだありません。
