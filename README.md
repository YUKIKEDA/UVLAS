# UVLAS

実務で利用できる品質のUVLM–ANCFシェル連成ソルバーを開発するプロジェクト。

## 現在地

数値基盤として、密・疎行列、密・疎の基準解法、GMRES、固有値抽出、精度保持と明示投影、状態・所有権の契約、Python結合、計測基盤を実装した。凍結した設計試作の参照を製品経路へ接続し、[IM-00の製品検証の終了条件を満たした](docs/plans/im00-completion-review.md)。次はPG・疎拘束、有限渦kernel、入力・配布基盤の個別契約を並列に具体化する。**物理kernelは未実装。** UVLM–ANCFの物理系の組立・連成、S10実験参照の認定、製品規模の性能、全実務受入は未完了。

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
