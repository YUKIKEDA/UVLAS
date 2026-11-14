# 文書の索引

製品開発で読む文書はここにある。R&Dの経緯と採用しなかった候補は [rd/](../rd/README.md) にある。

## まず読むもの

| 文書 | 内容 |
|---|---|
| [roadmap.md](roadmap.md) | 作業の順序と現在地。行のIDとIssueの対応 |
| [conventions.md](conventions.md) | ディレクトリ構成とモジュール |
| [build.md](build.md) | 構成、ビルド、試験の手順 |
| [engineering/](engineering/README.md) | 共通Rules（SE-01〜14）、知見と根拠、レビュー・プロファイル |
| [adr/](adr/) | 決定の記録 |
| [quality/quality-plan.md](quality/quality-plan.md) | 品質保証計画（QA-01〜10）：追跡、基準結果、結果の変更管理、版の意味、来歴、環境の認定、公開の審査 |
| [quality/traceability.md](quality/traceability.md) | 追跡表：各IDの定義、試験、検証記録（道具が書く） |
| [quality/case-status.md](quality/case-status.md) | 受入ケース19件の状態 |
| [quality/result-changes/](quality/result-changes/README.md) | 結果変更記録 |

## 要求と受入：何を作るか、何をもって完了とするか

| 文書 | 内容 |
|---|---|
| [product-requirements.md](requirements/product-requirements.md) | 製品要求。対象、範囲、利用経路 |
| [project-start.md](requirements/project-start.md) | 利用者決定の台帳（Q1〜Q27、D00〜D25） |
| [acceptance-criteria.md](requirements/acceptance-criteria.md) | 受入精度と証拠の条件 |
| [practical-acceptance-plan.md](requirements/practical-acceptance-plan.md) | 受入計画A5-01〜07。19ケース、S10、P-CURVED-BOX、性能 |
| [performance-plan.md](requirements/performance-plan.md) | 性能試験の方針。代表ケースの具体条件は受入計画A5-06 |
| [verification-resources.md](requirements/verification-resources.md) | 公開参照の候補 |

## 採用仕様：実装の正本

優先関係は [レビュー・プロファイル](engineering/uvlas-review-profile.md) に従う。初期微分と受理後の荷重は、J1の追加規約が優先する。

| 領域 | 文書 | 決定ID |
|---|---|---|
| 系全体 | [coupled-problem-definition.md](specs/coupled-problem-definition.md)、[numerical-conventions.md](specs/numerical-conventions.md) | |
| 構造 | [structural-adopted-specification.md](specs/structure/structural-adopted-specification.md) | ST-01〜10 |
| 主シェル | [main-shell-adopted-specification.md](specs/structure/main-shell-adopted-specification.md) | SH-D01〜08 |
| 材料・積層 | [material-shell-adopted-specification.md](specs/structure/material-shell-adopted-specification.md)、[構成契約](specs/structure/material-layer-adoption-contract.md)、[初期ひずみ契約](specs/structure/material-initial-strain-contract.md) | MAT-D01〜13 |
| 梁 | [beam-adopted-specification.md](specs/structure/beam-adopted-specification.md) | |
| 接続・支持・慣性 | [connection-adopted-specification.md](specs/structure/connection-adopted-specification.md) | |
| 空力 | [uvlm-integrated-specification.md](specs/aero/uvlm-integrated-specification.md)、[誘導](specs/aero/uvlm-induction-adopted-specification.md)、[後流・荷重](specs/aero/uvlm-wake-adopted-specification.md) | AE-01〜07 |
| 連成・時間・安定性 | [coupled-adopted-specification.md](specs/coupled/coupled-adopted-specification.md)、[追加規約](specs/coupled/j1-initialization-output-derivative-addendum.md) | J1-D01〜09 |
| 連成の補助契約 | [離散平衡とchart](specs/coupled/j1-discrete-equilibrium-and-chart-contract.md)、[静的継続](specs/coupled/j1-static-continuation-error-contract.md)、[時間・安定性](specs/coupled/coupled-time-stability-specification.md)、[初期化残差](specs/coupled/coupled-initialization-residual-contract.md)、[状態と探索](specs/coupled/coupled-state-search-contract.md) | |
| 決定の台帳 | [design-decisions.md](specs/design-decisions.md) | AE-01〜07ほか、各領域の決定の定義 |
| 数値基盤・型・I/F | [numerical-foundation-adoption-specification.md](specs/foundation/numerical-foundation-adoption-specification.md)、[線形精度契約](specs/foundation/linear-solve-accuracy-contract.md) | N4-01〜10 |
| 入出力・実行・再開 | [io-execution-adoption-specification.md](specs/foundation/io-execution-adoption-specification.md) | IO4-01〜07 |

採用仕様が参照する導出と根拠（PG-ANCFの定義、UVLMの個別契約など）は、`rd/` の該当系統にある。仕様の文面を直すときは、追加規約か新しい版として足す。依存管理はN4-01を [ADR 0001](adr/0001-vcpkg-manifest-for-dependencies.md) が読み替えている。

## 計画

| 文書 | 内容 |
|---|---|
| [implementation-start-review.md](plans/implementation-start-review.md) | 開始審査。実装パッケージIM-00〜08、残リスクR-01〜07 |

実装前の計画、要求と証拠の対応表（当時のもの。現在は追跡表）、当時の現在地の記録は、履歴として [rd/planning/](../rd/planning/) にある。

## 報告書

[PG-ANCF技術報告書・第1稿](papers/pg-ancf-technical-report-v1.md)（[PDF](papers/pdf/pg-ancf-technical-report-v1.pdf)）。未査読。
