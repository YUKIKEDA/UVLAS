# 0005: kernel計測は独立したGoogle Benchmark実行ファイルで行う

- Status: Accepted
- Date: 2026-11-15
- Issue: #25

## Context

計測基盤は最初の数値kernelより前に使える必要がある。Google Benchmarkを使うことは開発規約とIssue #25で採用済みであり、このADRは具体的な版と導入境界を記録する。

## Decision

Google Benchmark 1.9.5をvcpkgのbenchmarks featureとoverrideで固定する。実物のportはApache-2.0を宣言し、ソース取得のSHA512を固定している。導入されたcopyrightと公開ヘッダを確認した。

Benchmarkは独立したuvlas_kernel_benchにだけリンクする。製品ライブラリ・CLI・wheelにリンクせず、通常の試験にも時間の閾値を設けない。featureの異なる構成は専用の導入prefixを使い、通常構成の導入物を変更しない。

## Consequences

最初の対象は製品のchecked_extent_productである。BenchmarkのJSONに機械条件・設定・build情報を含める。小kernelの時間から全ソルバーの実務性能を認定しない。CIではBenchmarkを実行せず、OSごとの性能採取は開発者が条件を固定して行う。
