# UVLAS

実務で利用できる品質のUVLM–ANCFシェル連成ソルバーを開発するプロジェクト。

<!-- active-uvlm-status:start -->
2026-11-06・最新：**実装前1〜4の採用設計、5の受入計画、6の最初の実装パッケージ審査を完了。次はIM-00：数値表現・求解・状態契約の製品実装。** 両OSの数値／Python接続を検証し、19ケース・37製品要求を対応付け、1072ソース・687成果物を監査。**S10実験参照の認定、実PG全規模・全自由後流の製品性能、全実務受入・リリースは未完了。本体コードは未着手。**

[開始判定・実装順と残リスク](docs/implementation-start-review.md)、[受入計画](docs/practical-acceptance-plan.md)、[ケース台帳](docs/verification/practical-acceptance-case-catalog.json)。
<!-- active-uvlm-status:end -->

## 現在地

方式・離散式・主要I/F・不変条件・参照・受入計画を整理し、最初の製品実装IM-00の開始条件を満たした。現在のsrcはビルド情報／CLI基盤のみ。設計実験のPython／C++コードは製品ソルバーではない。

[全体ロードマップ](docs/roadmap.md)、[開始審査とIM-00〜08](docs/implementation-start-review.md)から作業と終了条件を確認できる。実験資料の不足、性能、精度消費、外部格子変換等の残リスクは開始審査R-01〜07で管理する。

## 目標

低速・非圧縮・付着流を対象に、複数揚力面と後流の相互干渉、等方材・積層複合材のANCFシェルと梁を組み合わせた構造、集中質量・偏心慣性・剛体結合、静的双方向連成、非線形過渡応答、大変形平衡周りのフラッタ予測を扱う。CLIとPython APIを提供し、WindowsとLinuxに対応する。

代表性能目標は空力約5,000パネル・構造約10万自由度、単一PCでフラッタ速度探索全体12時間以内。主測定プロファイルは64GB搭載・プロセス48GiB。製品での達成実績はまだない。精度・適用範囲・測定条件は[受入計画](docs/practical-acceptance-plan.md)に定める。

外部格子はNastran系の限定サブセットとGmsh。対応形式・カード・変換・拒否・再開は[入出力仕様](docs/io-execution-adoption-specification.md)に定める。途中確定状態と探索からの再開は必須。

内部マイルストーンと公開リリースを区別する。初回の実務適用に必要な品質・機能を満たすまで部分公開しない。ソース、両OSのCLI、Python wheelを提供する計画。自身のソースはApache-2.0、第三者コード・データの条件は別に管理する。

## 採用仕様と根拠

- [製品要求](docs/product-requirements.md)、[精度・証拠条件](docs/acceptance-criteria.md)
- [構造統合](docs/structural-adopted-specification.md)、[主シェルPG-ANCF](docs/main-shell-adopted-specification.md)
- [UVLM統合](docs/uvlm-integrated-specification.md)、[連成・時間・安定性](docs/coupled-adopted-specification.md)
- [初期・荷重微分の最終規約](docs/j1-initialization-output-derivative-addendum.md)
- [数値基盤・型・I/F](docs/numerical-foundation-adoption-specification.md)、[入出力・実行・再開](docs/io-execution-adoption-specification.md)
- [19ケースの受入計画](docs/practical-acceptance-plan.md)、[実装開始判定](docs/implementation-start-review.md)
- [採用・不採用の経緯](docs/structural-model-adoption-review.md)、[設計台帳](docs/design-decisions.md)、[実装前計画](docs/preimplementation-plan.md)

## ビルド基盤

C++23。現在のビルド／CLI基盤の実行方法は[ビルド手順](docs/build.md)、両OSの確認は[C++23検証](docs/verification/cpp23-2026-10-29.md)。数値依存の限定試作は[4の終了判断](docs/verification/numerical-foundation-adoption-completion.md)。採用証拠の参照元は凍結し、製品へ移植した結果は別に検証する。

## 開発・レビュー用の再利用知識

[PG-ANCF技術報告書・第1稿](docs/papers/pg-ancf-technical-report-v1.md)（[PDF](output/pdf/pg-ancf-technical-report-v1.pdf)）に、採用した参照係数空間の定式化、恒等式の導出、保存結果の図表、先行研究との関係と未解決事項を統合した。未査読で、新規性や一般の厚さ一様収束を確定するものではない。

[Skills・Rulesと知見の入口](docs/engineering/README.md)。R&D、仕様・契約設計、実装レビューの3つのSkillと、共通Rules、UVLAS固有のレビュー規約、失敗を含む知見の根拠をまとめている。本体実装を別の担当者へ引き継ぎ、レビューに利用する。
