# 構成、ビルド、試験

現在の製品コードはビルド情報とCLIの骨組みだけで、解析コマンドは無い。公開APIの互換性保証は始めていない。

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
| nlohmann/json | 3.12.0 | ケースと状態のJSON |
| Catch2 | 3.16.1 | 試験（feature `tests`） |
| pybind11 | 3.1.0 | Python結合（feature `python`。現在のPresetでは導入しない） |

## Windows

Presetのジェネレータは Ninja Multi-Config で、コンパイラは環境の `cl` を使う（[ADR 0002](adr/0002-ninja-presets-with-the-msvc-environment.md)）。コマンドは [tools/dev/msvc.ps1](../tools/dev/msvc.ps1) を通して、MSVCの環境の中で実行する。このスクリプトは、Visual Studioの環境が自分の同梱版へ書き換える `VCPKG_ROOT` を、呼出し前の値へ戻す。

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
| `qa.numeric_canary` | 基準ケースの結果が、このプラットフォームの基準結果とビット一致すること（QA-04）。認定していないプラットフォームでは落ちる |
| `cli.build_info` | `uvlas --build-info` が来歴を返すこと |

ソルバーの精度を確かめる試験は、まだ無い。

## 開発用CLIの終了規約

| 入力 | 標準出力 | 標準エラー | 終了コード |
|---|---|---|---|
| `--help` | 利用方法と開発状態 | 空 | 0 |
| `--version` | 開発版、コンパイラ、プラットフォーム、ビルド構成 | 空 | 0 |
| `--build-info` | 来歴のJSON：版、commit、未commitの変更の有無、コンパイラ、浮動小数点の設定、命令セット、プラットフォームの鍵、依存の版 | 空 | 0 |
| 引数なし・余分な引数 | 空 | 引数個数の診断と利用方法 | 2 |
| 想定外の例外 | 不定 | 不定 | 1 |
| 未知のオプション | 空 | 未知のオプションの診断と利用方法 | 2 |

解析CLIの終了コードはIO4-01が定める。実装はIM-06。

## 確認の記録

| 日付 | 環境 | 確認したこと |
|---|---|---|
| 2026-10-29 | Windows／MSVC 14.44、Ubuntu 22.04／GCC 14.4 | 依存なしの構成で、Debug・ReleaseのビルドとCLI試験5件。[記録](../rd/foundation/build/records/cpp23-2026-10-29.md) |
| 2026-11-11 | Windows／MSVC 19.51（Visual Studio 2026）、CMake 3.30.4、vcpkg 2026-10-29（baseline `2ec0f174`） | vcpkg manifestでの構成、Debug・Releaseのビルド、試験8件、clang-format 21.1.2とclang-tidy 21.1.1の検査 |
| 2026-11-12 | GitHub Actions：ubuntu-24.04／GCC 14.2.0、windows-2025／MSVC | PR #28のCI。Linuxは `linux-debug`、`linux-release`、`linux-sanitize` のビルドと試験8件、Windowsは同じ8件、フォーマット、Skillの複製、追跡、R&Dツリーの検査 |

開発機にLinux（GCC 14）の環境は無い。Linuxの確認は、CIの実行による。
