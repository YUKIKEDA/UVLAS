# 0004: 製品の試験はCatch2で書く

- Status: Accepted
- Date: 2026-11-11
- Issue: #1

## Context

採用仕様は試験フレームワークを定めていない。これまでの製品試験は、CLIを別プロセスで起動するCMakeスクリプトだけだった。IM-00からは、型・解法・失敗契約を対象に、許容差付きの数値比較と負対照を多数書く。

## Decision

Catch2 v3（3.16.1）を使う。vcpkgのfeature `tests` で導入し、試験の実行ファイルにだけリンクする。試験は `uvlas_add_test()`（[tests/CMakeLists.txt](../../tests/CMakeLists.txt)）で登録し、CTestのラベルで種類と受入ケースを示す。

## Consequences

- 許容差付きの比較（`WithinAbs`、`WithinRel`）と、入力を並べた試験（generator）が使える。
- 1実行ファイルが1つのCTest試験になる。細かい失敗箇所はCatch2の出力で読む。
- Catch2は製品のライブラリ、CLI、wheelへリンクしない。配布物のライセンス通知に含めない。
