# ビルドと開発用CLI

現在はビルド基盤のみを実装している。CLIには `--help` と `--version` があり、解析コマンドは未実装。公開APIの互換性保証は開始していない。

## 必要な環境

- CMake 3.25以上。
- C++23を扱えるコンパイラと標準ライブラリ。構成時に `std::expected` のmonadic operations（`__cpp_lib_expected >= 202211L`）を検査する。
- Windowsの標準Preset：Visual Studio 2022のx64 C++ツール（MSVC 19.44以上）とWindows SDK。通常のPowerShellから実行可能。
- Linuxの標準Preset：GCCまたはClangとNinja。

Ubuntu 22.04標準のGCC 11は必要なライブラリ機能に未対応のため対象外。新しいコンパイラを選ぶ際は新規のビルドディレクトリを使うか、`CMAKE_CXX_COMPILER` を明示して再構成する。

MSVCではCMake 3.30・4.1のC++23指定が `/std:c++latest` に対応するため、`cmake/CompilerPolicy.cmake` で `/std:c++23preview` を明示する。後続規格の機能へ暗黙に依存しない方針とする。`std::mdspan` は現時点で必須条件にしておらず、必要な機能の採用時に対応条件を追加する。

現在の計算コードに外部数値ライブラリはない。数値依存は技術試作後に追加・版固定する。

## Windows

プロジェクトルートで実行する。

```powershell
cmake --preset windows-msvc
cmake --build --preset windows-debug
ctest --preset windows-debug
cmake --build --preset windows-release
ctest --preset windows-release
.\build\windows-msvc\Release\uvlas.exe --version
```

Windowsはマルチ構成のビルドディレクトリを使い、DebugとReleaseの生成物を分ける。手動実行時も構成のディレクトリを明示する。

## Linux

プロジェクトルートで実行する。

```sh
cmake --preset linux-debug
cmake --build --preset linux-debug
ctest --preset linux-debug
cmake --preset linux-release
cmake --build --preset linux-release
ctest --preset linux-release
./build/linux-release/uvlas --version
```

ローカルのコンパイラ指定などは、必要に応じて追跡対象外の `CMakeUserPresets.json` へ置く。共有Presetに開発機の絶対パスを保存しない。

## 開発用CLIの終了規約

| 入力 | 標準出力 | 標準エラー | 終了コード |
|---|---|---|---|
| `--help` | 利用方法と開発状態 | 空 | 0 |
| `--version` | 開発版、コンパイラ、プラットフォーム、ビルド構成 | 空 | 0 |
| 引数なし・余分な引数 | 空 | 引数個数の診断と利用方法 | 2 |
| 未知のオプション | 空 | 未知のオプションの診断と利用方法 | 2 |

数値計算のエラー分類・終了コードは、解析CLIを追加する段階で定義する。

## 試験

CTestはCLIを別プロセスで起動し、成功時のビルド情報、不正入力時の終了コード、標準出力と標準エラーの使い分けを確認する。Presetは試験が0件の場合もエラーにする。

初期C++20の結果は [T01検証記録](verification/build-2026-10-29.md)、C++23移行後の結果は [移行検証記録](verification/cpp23-2026-10-29.md) に記録した。ソルバーの精度を検証する試験はまだ含まない。

今回検証した隔離ツールチェーンと再現用パッケージ仕様は [ツールチェーン手順](../tools/toolchains/README.md) を参照。
