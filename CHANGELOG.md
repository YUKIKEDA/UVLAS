# 変更履歴

公開した版を、新しい順に、[Keep a Changelog](https://keepachangelog.com/ja/1.1.0/) の形で記す。版の意味は [品質保証計画 QA-05](docs/quality/quality-plan.md) が定める：PATCHは結果を変えない。MINORは結果を変えてよく、変えた変更は全て「結果に影響する変更」に、[結果変更記録](docs/quality/result-changes/README.md) の番号と一緒に載る。

公開した版は、まだ無い。

## [Unreleased]

### 結果に影響する変更

解析の能力が、まだ無いため、解析結果への影響はない。

- [0001](docs/quality/result-changes/0001-linux-gcc14-platform.md)：プラットフォーム `linux-x64-gnu14.2-sse2` を認定した。標準数学ライブラリの2関数が、Windowsと1ulp異なる。

### 追加

- ビルドの来歴を返す `uvlas --build-info`。

### 既知の問題

- 解析コマンドは実装されていない。
