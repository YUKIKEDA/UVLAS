---
name: reference-fixture
description: UVLASの凍結R&D参照を、原本hashと独立oracleを保って製品fixtureへ取り込み、別directoryで再生成・検証する。数値方式の選定やR&D結果の上書きには使わない。
---

# 凍結参照を取り込む

製品試験が実際に読むfixtureを扱う。参照を製品の計算結果で作らず、元の記録・入力・生成sourceと変換の境界を保つ。契約は [fixture契約](../../../docs/specs/foundation/reference-fixture-contract.md)。

## 参照を探す

引用名がbuild/でも、直接のbuildだけで現物の有無を判定しない。`rd/evidence-index.json` のSHAとbundleを調べ、`tools.rd.data_store.cached(name)` の正規cacheを確認する。未取得なら既存evidence工具で必要なbundleを取得する。元のcacheはread-onlyで扱う。

この流れの実例はPG24参照の3artifactである。元audit・索引・実体hashを照合する。生成ソースも元audit・現在file・fixtureの3者で同定する。外部原本が見つからない場合、似た合成問題を旧id/hashへ差し替えない。

## 別directoryへ作る

repo rootからmoduleとして実行する。

```text
python -m tools.qa.foundation_reference import --output out/new-pg-reference
python -m tools.qa.foundation_reference verify out/new-pg-reference
```

別cacheを使うときだけ `--cache-root <cache-root>` をimportへ指定する。既存outputは拒否される。3つの原本のバイト列は保持し、派生JSONはUTF-8/LFとする。引用済みの原本に改行変換をかけない。作業ファイルのバイト列と、gitのblobから読み出した内容のSHA-256が一致することを、初回commit前に確認する。

このimportは物理R&Dを再実行せず、原Pを読んだ代数adapterと独立Decimal参照だけを作る。原本の物理生成をやり直す必要がある場合は、別のscratch copyと新しい参照idを使い、凍結した旧参照と比較する。元R&D sourceやcacheの上書きはしない。

## 比較して渡す

- 新fixtureを採用する前に、既存fixtureとの差をbytesとJSON値で示す。値が動いたら原因と許容差への影響を判断し、結果変更の規則へ従う。
- verifyと、fixtureを読む対象の製品試験を実行する。工具は標準libraryで、ruff/pyright strictを通す。無関係な数値試験を再実行しない。
- hash不一致、shape/order/finite、boolの数値化、複素pairの各成分、元のソースとfixture hashの同時付替え、既存outputへの上書きを負対照にする。
- 元のid/sha/hintとmap型sources/artifactsを記録する。絶対pathや利用者名を公開せず、必要なら既存redactorで別copyを作り前後hashを残す。
- fixtureの一致と製品solverの合格を分ける。未実装の固有値・Arnoldiや未実行のOSを達成扱いにせず、Issue #15を準備PRで閉じない。

Skillを修正したら `.claude/skills` の鏡像も同じcommitで同期する。
