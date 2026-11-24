# 構成、ビルド、試験

現在の製品コードはビルド情報、共通の失敗・誤差認定、checked寸法・資源計画、工程時間・メモリ計測とCLIの骨組みで、解析コマンドは無い。公開APIの互換性保証は始めていない。

## 必要なもの

| 道具 | 条件 |
|---|---|
| CMake | 3.25以上 |
| Ninja | Windows、Linuxとも |
| vcpkg | manifestの `builtin-baseline` のコミットを持つ版。環境変数 `VCPKG_ROOT` を設定する |
| Windows | Visual StudioのC++ツール（MSVC 19.44以上）。2026-11-11の確認はVisual Studio 2026（MSVC 19.51） |
| Linux | GCC 14と、対応する標準ライブラリ。認定済みは GCC 14.2（`linux-x64-gnu14.2-sse2`）。認定していないコンパイラの版では、基準結果の試験が落ちる |
| just | 任意。レシピは [justfile](../justfile) |
| Python | 3.10以上。試験（基準結果の照合）が使う。R&Dツリーの検査と追跡の道具は、`uv` が用意する3.12で走る |
| uv | `just` のレシピが、clang-format、clang-tidy、Python 3.12を呼ぶのに使う |

構成時に、`std::expected` のmonadic operations（`__cpp_lib_expected >= 202211L`）を検査する。Ubuntu 22.04標準のGCC 11は対象外。MSVCでは `cmake/CompilerPolicy.cmake` が `/std:c++23preview` を明示する。

## 依存

第三者ライブラリはvcpkgのmanifest（[vcpkg.json](../vcpkg.json)）が固定する（[ADR 0001](adr/0001-vcpkg-manifest-for-dependencies.md)）。CMakeの構成時に、Presetが渡すvcpkgのtoolchainが `out/vcpkg_installed/` へ導入する。初回はBoostなどの導入に時間がかかる（2026-11-11の開発機で約25分）。

| ライブラリ | 版 | 用途 |
|---|---|---|
| Eigen | 5.0.1 | 密・疎の線形代数（MPL-2.0部分のみ） |
| Spectra | 1.2.0（overlay port） | 非対称・複素固有値 |
| Boost Multiprecision、Hash2 | 1.92.0 | 精密参照、正確な線形拘束、SHA-256 |
| nlohmann/json | 3.12.0（UVLAS例外補修1、port #3） | ケースと状態のJSON |
| Catch2 | 3.16.1 | 試験（feature `tests`） |
| Google Benchmark | 1.9.5 | kernel計測専用（feature `benchmarks`、通常Presetは導入しない） |
| pybind11 | 3.1.0 | Python結合（feature `python`。現在のPresetでは導入しない） |

## Windows

Presetのジェネレータは Ninja Multi-Config で、コンパイラは環境の `cl` を使う（[ADR 0002](adr/0002-ninja-presets-with-the-msvc-environment.md)）。コマンドは [tools/dev/msvc.ps1](../tools/dev/msvc.ps1) を通して、MSVCの環境の中で実行する。このスクリプトは、Visual Studioの環境が自分の同梱版へ書き換える `VCPKG_ROOT` を、呼出し前の値へ戻す。

Windowsの製品CIでは、入力パーサーの確保量モデルを認定したMSVC 14.51（STL 145、202604）を使う。[CIの選択スクリプト](../tools/dev/ci-msvc.ps1)でコンパイラとSTL版をビルド前に照合し、製品とvcpkgに同じVisual Studioを指定する。runnerの既定が14.52へ切り替わった際は、入力パーサーが未認定の環境を拒否する分岐でMSVCの到達不能コード警告が発生した。認定済みの版がrunnerに無い場合は、固定したVisual Studio Build Toolsの配布物を使う。

```powershell
pwsh tools/dev/msvc.ps1 cmake --preset windows-msvc
pwsh tools/dev/msvc.ps1 cmake --build --preset windows-debug
pwsh tools/dev/msvc.ps1 ctest --preset windows-debug
pwsh tools/dev/msvc.ps1 cmake --build --preset windows-release
pwsh tools/dev/msvc.ps1 ctest --preset windows-release
.\out\build\windows-msvc\Release\uvlas.exe --version
```

同じことを `just gate` が行う。

## Linux

```sh
cmake --preset linux-debug   && cmake --build --preset linux-debug   && ctest --preset linux-debug
cmake --preset linux-release && cmake --build --preset linux-release && ctest --preset linux-release
cmake --preset linux-sanitize && cmake --build --preset linux-sanitize && ctest --preset linux-sanitize
```

`linux-sanitize` はAddressSanitizerとUndefinedBehaviorSanitizerを付けたDebug。コンパイラの指定など手元だけの設定は、追跡対象外の `CMakeUserPresets.json` に置く。共有Presetに開発機の絶対パスを書かない。

## ビルドの方針

- 生成物は `out/build/<preset>/`。`build/` はR&D時代の作業領域で、証拠が入っているため製品のビルドに使わない。
- 警告はエラー（`UVLAS_WARNINGS_AS_ERRORS`、Presetで有効）。
- 浮動小数点は、MSVCが `/fp:strict`、GCC／Clangが `-ffp-contract=off -fno-fast-math`。two-productと補償加算（N4-06）の前提であり、設計作業4の数値試作と同じ設定。
- Eigenは `EIGEN_MPL2_ONLY` と `EIGEN_DONT_PARALLELIZE`（N4-01）。

## 整形と静的検査

| コマンド | 内容 |
|---|---|
| `just format-check` | clang-format 21.1.2。`src/`、`include/`、`tests/` |
| `just lint` | clang-tidy 21.1.1。構成済みPresetのコンパイル情報（`compile_commands.json`）を使い、`src/` と `tests/` を検査する。検査項目と、対応する規則は [.clang-tidy](../.clang-tidy) に書いてある |

| `just lint-python` | `tools/` のPython。ruff 0.17.4（検査と整形）、pyright 1.1.408のstrict |

どれもゲートで、指摘はエラーになる。先に `cmake --preset …` で構成しておく。

## 試験

試験はCatch2で書き、`uvlas_add_test()`（[tests/CMakeLists.txt](../tests/CMakeLists.txt)）で登録する。CTestのラベルは `unit`、`property`、`integration`、`verification`、`negative` と、対応する受入ケース。Presetは、試験が0件のときもエラーにする。

| 試験 | 内容 |
|---|---|
| `cli.*`（5件） | CLIを別プロセスで起動し、ビルド情報、不正入力時の終了コード、標準出力と標準エラーの使い分けを確かめる |
| `toolchain.pinned_dependencies` | コンパイラへ届く依存の版がmanifestと一致すること、数値のビルド方針が効いていること、`fma` が積の誤差を正確に返すこと |
| `toolchain.linear_algebra` | 実・複素のpivot付きLU、64bit CSCのSparseLU＋COLAMD、密・疎LLT、固定前処理GMRESを既知解と原残差で照合する |
| `toolchain.arnoldi` | 実数作用素の複素共役固有値と、複素非Hermitian作用素の固有値を、行列を渡さないArnoldiで求める |
| `toolchain.precise_arithmetic` | Multiprecisionの精度、有理数の厳密算術、Hash2のSHA-256既知値を照合する |
| `toolchain.third_party_notices` | vcpkgの著作権・ライセンス通知が構成先へバイト一致で複製されていること。扱いは[第三者通知](../THIRD_PARTY_NOTICES.md) |
| `toolchain.notice_configuration` | 実ソースの複製でCatch2通知の6構成と、Benchmarkの既定OFF・明示OFF・ON・同一ビルド先のON→OFFを確認する。後者は登録済みの通知CTestも走らせ、ONの複製・上流通知欠落を拒否する。構成専用package fixtureを使い、依存原本の変更・C++ビルドは行わない |
| `core.failure` | 17種類の失敗、所有する診断、識別子の混用拒否、例外をResultへ変換する入口の契約 |
| `core.error_certificate` | 未評価・達成・未達の区別、予算境界、負値・非有限・未対応表現の拒否 |
| `core.checked_dimensions` | 符号付き64bit寸法・バイト数の検査、ゼロ次元、32bit外部APIへの検査付き変換 |
| `core.resource_plan` | 累積payloadの確保前照合、超過時の確保禁止、失敗時の予約解放と内側成功の計上保持、所有配列の寿命 |
| `core.identifier_compilation` | 正のコンパイル対照と、8種類のIDを異種IDから受け取るコンパイル失敗の確認 |
| `core.instrumentation` | 工程集計、無効時の結果不変、OSメモリ、条件と結果のJSON、入力拒否 |
| `core.allocation_audit` | C++の全new形式、workerでの確保、未対応監査の拒否、準備済みEigen演算とthread数別の結果比較 |
| `core.eigen_allocation_guard` | 同じprobeの正常対照と、Eigenの動的確保・reallocの期待診断による拒否 |
| `core.audit_configuration` | 未監査の製品ライブラリの直接・interface越しリンクを構成時に拒否 |

| `core.state` | 不変Model、親・候補・資格の所有とhash照合、一回確定、未達・未評価・変更後の拒否、全cache世代、Workspaceの再入・並行排他 |
| `numerics.storage` | 密・complex・CSC疎・PointSoA・DofLayoutの配置、所有、ID・写像・尺度・非有限の検査 |
| `numerics.eigen_views` | 非公開Eigen境界のpointer一致、padding、scalar・配置・stride・shape・非有限の拒否 |
| `numerics.storage_compilation` | 正のコンパイル対照、未対応scalarとtemporaryからの借用のコンパイル失敗 |
| `qa.numeric_canary` | 基準ケースの結果が、このプラットフォームの基準結果とビット一致すること（QA-04）。認定していないプラットフォームでは落ちる |
| `cli.build_info` | `uvlas --build-info` が来歴を返すこと |
| `cli.validate` | ケース検証の成功、構文と参照の失敗、資源上限を実行ファイルの終了コードと診断で確認する |
| `cli.inspect` | 保存済みの状態と来歴を読み、指紋の不一致、欠落、資源上限を実行ファイルで確かめる |

ソルバーの精度を確かめる試験は、まだ無い。

確保監査では通常版のuvlas_coreを混ぜず、同じ製品sourceを再コンパイルしたuvlas_core_auditを使う。Eigenのruntime確保検査とassertを全TUで揃え、Releaseでも監査用targetはNDEBUGを解除する。監視backendとthreadの範囲は [計測契約](specs/foundation/instrumentation-contract.md) に従い、未対応の経路を確保0として通さない。

kernelの計測は [benchmarks/README.md](../benchmarks/README.md) の独立した構成で行う。Benchmarkと通常構成の導入prefixを分ける。CPUプロファイルの手順は [kernel-profile](../.agents/skills/kernel-profile/SKILL.md) を使う。

## 開発用CLIの終了規約

| 入力 | 標準出力 | 標準エラー | 終了コード |
|---|---|---|---|
| `--help` | 利用方法 | 空 | 0 |
| `--version` | 開発版、コンパイラ、プラットフォーム、ビルド構成、命令セット | 空 | 0 |
| `--build-info` | 来歴のJSON：版、commit、未commitの変更の有無、コンパイラ、浮動小数点の設定、命令セット、プラットフォームの鍵、依存の版 | 空 | 0 |
| `validate case.json` | 構文・スキーマ・ケース内参照の検査完了を表示 | 失敗した工程と、取得できた元入力上の位置 | 成功は0、入力不正・未対応は2、資源上限は5 |
| `inspect result-dir` | 照合済みの `status.json` | 取得や照合に失敗したファイル、または資源上限の診断 | 照合成功は0、引数不正は2、資源上限は5、保存結果の不整合は6 |
| 引数なし・余分な引数 | 空 | 引数個数の診断と利用方法 | 2 |
| 想定外の例外 | 不定 | 不定 | 1 |
| 未知のオプション | 空 | 未知のオプションの診断と利用方法 | 2 |

`validate` は読取り専用のコマンドである。全スキーマとケース内の参照を検査し、外部参照があればマニフェスト直下のファイルを宣言順に処理する。ファイルごとに元バイト列の SHA-256 を照合する。MSH 4.1 ASCII は、照合したバイト列から節とタグを検査し、メッシュ本体に子参照がないことを確かめてから次のファイルへ進む。Gmsh の sidecar など、形式別の依存関係を検査できない宣言は、その位置を診断して終了コード 2 を返す。物理準備、実行、再開、結果の確認と実行終了時の状態記録は IM-06 の後続工程で扱う。

開発用の入力受入上限は、マニフェスト 8 MiB、外部参照 1 ファイル 32 MiB、外部ファイル全体 64 MiB、所有ペイロードの累積確保 256 MiB。パーサーの外部確保には別に 256 MiB の見積もり上限を適用する。ケースの入力制限も型検査時に適用する。

`inspect` は一つの結果ディレクトリから `status.json` と `provenance.json` を読み、状態の各項目と来歴の SHA-256 を照合する。来歴が `resolved-case.json` の指紋を持つ場合は、そのファイルのバイト列も照合する。照合できた状態は、元の JSON をそのまま標準出力へ返す。出力の `exit_code` は保存された解析の終了コードで、コマンド自身の終了コードは読取りと照合の結果を表す。状態・来歴は各 1 MiB、解決済みケースは 64 MiB、所有ペイロードの累積確保は 128 MiB を上限とする。

## 確認の記録

| 日付 | 環境 | 確認したこと |
|---|---|---|
| 2026-10-29 | Windows／MSVC 14.44、Ubuntu 22.04／GCC 14.4 | 依存なしの構成で、Debug・ReleaseのビルドとCLI試験5件。[記録](../rd/foundation/build/records/cpp23-2026-10-29.md) |
| 2026-11-11 | Windows／MSVC 19.51（Visual Studio 2026）、CMake 3.30.4、vcpkg 2026-10-29（baseline `2ec0f174`） | vcpkg manifestでの構成、Debug・Releaseのビルド、試験8件、clang-format 21.1.2とclang-tidy 21.1.1の検査 |
| 2026-11-12 | GitHub Actions：ubuntu-24.04／GCC 14.2.0、windows-2025／MSVC | PR #28のCI。Linuxは `linux-debug`、`linux-release`、`linux-sanitize` のビルドと試験8件、Windowsは同じ8件、フォーマット、Skillの複製、追跡、R&Dツリーの検査 |

開発機にLinux（GCC 14）の環境は無い。Linuxの確認は、CIの実行による。

## Python結合の任意構成

`UVLAS_BUILD_PYTHON=ON` とvcpkgの `tests;python` featureを、通常構成とは別の導入prefixへ指定する。外部CPython 3.12のInterpreter・Development.ModuleとNumPy 2.xを使う。構築例とAPIの使い方は [Python拡張](../src/python/README.md)。通常の `just lint` は任意moduleを除き、その構成のcompile databaseを使う専用のclang-tidyと `python.ownership` をWindows/LinuxのCIで実行する。ケースのrun/resumeとwheel配布はIM-06。

## 非公開のランタイム用ディレクトリ

ReleaseのCLI、CPython拡張、ライセンスを一つの場所へ配置し、ビルド先と異なる場所で実行できる。CPython 3.12と実行環境向けNumPy wheelを用意し、[構成・作成・移設の手順](verification/IM-06/private-runtime-image-producer-review.md)に従う。検証には新規venvを使う。

## 疎求解の対象検証

`numerics.sparse_spd` は外向きSPD bound、`numerics.sparse_triangular` は抽出・両置換・非重複span、`numerics.sparse_solve` は元CSCの求解と認定、`numerics.sparse_allocation` は準備済み三角kernelのC++/Eigen確保を検査する。`numerics.sparse_capacity` は10万未知数の合成三重対角で、準備・求解時間とprocess lifetime peak RSSを出す。容量試験を実務graphの性能やメモリprofileの認定に転用しない。

確保監査は製品のnumerics全体を監査用policyで再compileし、依存するcoreも同じpolicyの版を使う。通常版のEigen inline定義と混ぜない。

## 凍結参照fixtureの対象検証

`numerics.foundation_reference` は、既存PG24状態行列の代数adapterを製品DensePreparedで処理し、独立なDecimal参照と比較する。`tools.foundation_reference` は取込み・来歴・上書き拒否・全複素成分の負対照、`tools.foundation_reference_verify` は凍結fixtureの検証である。Pythonはrepo rootからmoduleとして起動し、CTestも同じ作業directoryを明示する。

```text
python -m tools.qa.foundation_reference verify tests/numerics/fixtures/j1-pg-free-history-24-v1
ctest --test-dir <configured-build> -C Debug -R "^(numerics.foundation_reference|tools.foundation_reference.*)$" --output-on-failure
```

この先行整備はIssue #15全体のgateではない。製品固有値・Arnoldi、GMRESを含む未達は [fixture契約](specs/foundation/reference-fixture-contract.md) に記す。取込み・別directoryでの再生成はreference-fixture Skillに従う。
