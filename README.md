# UVLAS

実務で利用できる品質のUVLM–ANCFシェル連成ソルバーを開発するプロジェクト。

## 現在地

2026-11-11：方式と契約の採用（構造、空力、連成・時間・安定性、数値基盤、入出力）と受入計画、最初の実装パッケージIM-00の開始審査を終え、リポジトリの基盤を整えた段階。**製品のソルバー本体は未着手。** `src/` にあるのはビルド情報とCLIの骨組みだけである。S10実験参照の認定、製品規模の性能、全実務受入は未完了。

作業の順序と状態は [ロードマップ](docs/roadmap.md)、実装の順序と残リスクは [開始審査](docs/plans/implementation-start-review.md) にある。

## 目標

低速・非圧縮・付着流を対象に、複数揚力面と後流の相互干渉、等方材・積層複合材のANCFシェルと梁を組み合わせた構造、集中質量・偏心慣性・剛体結合、静的双方向連成、非線形過渡応答、大変形平衡周りのフラッタ予測を扱う。CLIとPython APIを提供し、WindowsとLinuxに対応する。

代表性能目標は空力約5,000パネル・構造約10万自由度、単一PCでフラッタ速度探索全体12時間以内。主測定プロファイルは64GB搭載・プロセス48GiB。製品での達成実績はまだない。精度・適用範囲・測定条件は [受入計画](docs/requirements/practical-acceptance-plan.md) に定める。

外部格子はNastran系の限定サブセットとGmsh。途中確定状態と探索からの再開は必須。対応形式と拒否の条件は [入出力仕様](docs/specs/foundation/io-execution-adoption-specification.md) に定める。

内部マイルストーンと公開リリースを区別する。初回の実務適用に必要な品質と機能を満たすまで、部分公開しない。配布はソース、両OSのCLI、Python wheel。自身のソースはApache-2.0で、第三者コードとデータの条件は別に管理する。

## 文書

| 知りたいこと | 場所 |
|---|---|
| 文書の索引（要求、採用仕様、計画） | [docs/README.md](docs/README.md) |
| 作業の順序と状態 | [docs/roadmap.md](docs/roadmap.md) |
| 構成、ビルド、試験 | [docs/build.md](docs/build.md)、[docs/conventions.md](docs/conventions.md) |
| 開発の手順 | [CONTRIBUTING.md](CONTRIBUTING.md)、エージェント向けは [AGENTS.md](AGENTS.md) |
| R&Dの経緯、採用しなかった候補、証拠 | [rd/README.md](rd/README.md) |
| Rules、知見と根拠 | [docs/engineering/README.md](docs/engineering/README.md) |

[PG-ANCF技術報告書・第1稿](docs/papers/pg-ancf-technical-report-v1.md)（[PDF](docs/papers/pdf/pg-ancf-technical-report-v1.pdf)）に、採用した参照係数空間の定式化、恒等式の導出、保存結果の図表、先行研究との関係と未解決事項をまとめている。未査読で、新規性や一般の厚さ一様収束を確定するものではない。

## ビルド

C++23、CMake、Ninja、vcpkg。手順は [docs/build.md](docs/build.md)。

```text
just gate          # 構成、ビルド、試験
just format-check  # 整形の検査
```
