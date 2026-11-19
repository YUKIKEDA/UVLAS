# PR #186：訂正した対照の観測値を公開前に照合する

2026-11-19。指摘 `discussion_r4629082700`、対象head `d832d3ecc72015905cbfa8df7d8419a3f658459c`。現行runnerの公開前検査だけを訂正し、固定生成核と48系列の数値証拠は変更しない。

旧discovery.exeと訂正版は、どちらもcontrolsが終了コード0で `bad_prime_rank_loss_rejected=true` を返した。実際の両出力を `review-186-controls-fixture.json` に保存した。旧版には有理rank2/剰余rank1の観測項目がない。数値transportをstub化して旧入口へ渡すと、出力が既存の場合も未作成の場合も、旧形式のまま48行の結果を公開することを先に再現した。

現行 `run-discovery.py.txt` は、controlsがobjectであり、`bad_prime_rational_rank` がintの2、`bad_prime_modular_rank` がintの1、拒否flagがboolのTrueであることを検査する。boolをintとして受け入れる検査にはしない。条件は明示的な失敗判定であり、`python -O` でも有効である。合わなければprimitiveや数値処理へ進まず、完成結果を置換しない。

公開処理のprobeは、通常と最適化の各38系列で次を確認した。

- 正しい実controls形式だけが48のstub行を公開する。
- 旧形式、各rankの欠落・誤値・float/bool/string/null、拒否flagのFalse/整数/null/欠落、object以外を拒否する。
- 拒否時はprimitiveと数値transportの呼出しが0回で、既存byteまたは未作成状態を保持し、tempを残さない。

これは公開条件の検証であり、数値生成の再実行ではない。通常版 `review-186-controls-gate-result.json`、最適化版 `review-186-controls-gate-optimized-result.json` を保存した。既存のtesting規約に従って対照の名前でなく観測値を確認し、evidence規約に従って全件合格時だけ公開する。同趣旨の規約は増やさない。

旧runnerと旧manifestは `review-186-original-runner.py.txt` と `review-186-original-manifest.json.txt` にbyte保全した。既存manifestのrunner引用先だけをその原本コピーへ結び直し、hash値を保持した。新runnerの合格を旧48実行の証拠へ転記していない。生成核・48数値結果・元検証記録のbyteは元headと一致する。

```text
python -B review-186-controls-gate-probe.py.txt --legacy review-186-original-runner.py.txt --runner run-discovery.py.txt --fixture review-186-controls-fixture.json --output <新規結果>
python -B -O review-186-controls-gate-probe.py.txt --legacy review-186-original-runner.py.txt --runner run-discovery.py.txt --fixture review-186-controls-fixture.json --output <新規最適化結果>
```

独立レビュー担当が軽量対照の結果を確認した。CIや広い48数値系列は重複実行せず、PRとthreadは利用者の確認まで開いたままにする。親#140の未達条件は不変である。
