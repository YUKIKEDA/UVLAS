# 追跡表

`python tools/qa/trace.py matrix` が書く。手で直さない。IDの定義は文書から、試験はCatch2のタグから、検証記録は `docs/verification/` の `ids` から読む。一つのIDを調べるには `python tools/qa/trace.py show <ID>`。

「試験」「記録」の数が0であることは、そのIDがまだ実装も検証もされていないことを示す。数があることは、合格を意味しない。合否は検証記録と [品質保証計画](quality-plan.md) の状態で言う。

## 製品要求

| ID | 定義 | 試験 | 記録 | コード |
|---|---|---:|---:|---:|
| API-01 | [product-requirements.md:38](../requirements/product-requirements.md) | 0 | 0 | 0 |
| API-02 | [product-requirements.md:39](../requirements/product-requirements.md) | 0 | 0 | 0 |
| BC-01 | [product-requirements.md:26](../requirements/product-requirements.md) | 0 | 0 | 0 |
| DIST-01 | [product-requirements.md:45](../requirements/product-requirements.md) | 0 | 0 | 0 |
| DIST-02 | [product-requirements.md:46](../requirements/product-requirements.md) | 0 | 0 | 0 |
| ENV-01 | [product-requirements.md:40](../requirements/product-requirements.md) | 0 | 0 | 0 |
| ENV-02 | [product-requirements.md:41](../requirements/product-requirements.md) | 0 | 0 | 0 |
| GEO-01 | [product-requirements.md:19](../requirements/product-requirements.md) | 0 | 0 | 0 |
| IO-01 | [product-requirements.md:33](../requirements/product-requirements.md) | 0 | 0 | 0 |
| IO-02 | [product-requirements.md:34](../requirements/product-requirements.md) | 0 | 0 | 0 |
| IO-03 | [product-requirements.md:35](../requirements/product-requirements.md) | 0 | 0 | 0 |
| LOAD-01 | [product-requirements.md:30](../requirements/product-requirements.md) | 0 | 0 | 0 |
| MAT-01 | [product-requirements.md:24](../requirements/product-requirements.md) | 0 | 0 | 0 |
| MAT-02 | [product-requirements.md:25](../requirements/product-requirements.md) | 0 | 0 | 0 |
| OUT-01 | [product-requirements.md:31](../requirements/product-requirements.md) | 0 | 0 | 0 |
| OUT-02 | [product-requirements.md:32](../requirements/product-requirements.md) | 0 | 0 | 0 |
| PERF-01 | [product-requirements.md:42](../requirements/product-requirements.md) | 0 | 0 | 0 |
| PERF-02 | [product-requirements.md:43](../requirements/product-requirements.md) | 0 | 0 | 0 |
| PERF-03 | [product-requirements.md:44](../requirements/product-requirements.md) | 0 | 0 | 0 |
| PHY-01 | [product-requirements.md:17](../requirements/product-requirements.md) | 0 | 0 | 0 |
| PHY-02 | [product-requirements.md:18](../requirements/product-requirements.md) | 0 | 0 | 0 |
| QUAL-01 | [product-requirements.md:60](../requirements/product-requirements.md) | 0 | 0 | 0 |
| QUAL-02 | [product-requirements.md:61](../requirements/product-requirements.md) | 0 | 0 | 0 |
| QUAL-03 | [product-requirements.md:62](../requirements/product-requirements.md) | 0 | 0 | 0 |
| QUAL-04 | [product-requirements.md:63](../requirements/product-requirements.md) | 0 | 0 | 0 |
| QUAL-05 | [product-requirements.md:64](../requirements/product-requirements.md) | 0 | 0 | 0 |
| QUAL-06 | [product-requirements.md:65](../requirements/product-requirements.md) | 0 | 0 | 0 |
| QUAL-07 | [product-requirements.md:66](../requirements/product-requirements.md) | 0 | 0 | 0 |
| RUN-01 | [product-requirements.md:36](../requirements/product-requirements.md) | 0 | 0 | 0 |
| RUN-02 | [product-requirements.md:37](../requirements/product-requirements.md) | 0 | 0 | 0 |
| SOL-01 | [product-requirements.md:27](../requirements/product-requirements.md) | 0 | 0 | 0 |
| SOL-02 | [product-requirements.md:28](../requirements/product-requirements.md) | 0 | 0 | 0 |
| SOL-03 | [product-requirements.md:29](../requirements/product-requirements.md) | 0 | 0 | 0 |
| STR-01 | [product-requirements.md:20](../requirements/product-requirements.md) | 0 | 0 | 0 |
| STR-02 | [product-requirements.md:21](../requirements/product-requirements.md) | 0 | 0 | 0 |
| STR-03 | [product-requirements.md:22](../requirements/product-requirements.md) | 0 | 0 | 0 |
| STR-04 | [product-requirements.md:23](../requirements/product-requirements.md) | 0 | 0 | 0 |

## 採用仕様の決定

| ID | 定義 | 試験 | 記録 | コード |
|---|---|---:|---:|---:|
| AE-01 | [uvlm-integrated-specification.md:1](../specs/aero/uvlm-integrated-specification.md) | 0 | 0 | 0 |
| AE-02 | [design-decisions.md:55](../specs/design-decisions.md) | 0 | 0 | 0 |
| AE-03 | [design-decisions.md:56](../specs/design-decisions.md) | 0 | 0 | 0 |
| AE-04 | [design-decisions.md:57](../specs/design-decisions.md) | 0 | 0 | 0 |
| AE-05 | [design-decisions.md:58](../specs/design-decisions.md) | 0 | 0 | 0 |
| AE-06 | [design-decisions.md:59](../specs/design-decisions.md) | 0 | 0 | 0 |
| AE-07 | [uvlm-integrated-specification.md:73](../specs/aero/uvlm-integrated-specification.md) | 11 | 5 | 0 |
| BM-D01 | [beam-adopted-specification.md:9](../specs/structure/beam-adopted-specification.md) | 0 | 0 | 0 |
| BM-D02 | [beam-adopted-specification.md:10](../specs/structure/beam-adopted-specification.md) | 0 | 0 | 0 |
| BM-D03 | [beam-adopted-specification.md:11](../specs/structure/beam-adopted-specification.md) | 0 | 0 | 0 |
| BM-D04 | [beam-adopted-specification.md:12](../specs/structure/beam-adopted-specification.md) | 0 | 0 | 0 |
| BM-D05 | [beam-adopted-specification.md:13](../specs/structure/beam-adopted-specification.md) | 0 | 0 | 0 |
| BM-D06 | [beam-adopted-specification.md:14](../specs/structure/beam-adopted-specification.md) | 0 | 0 | 0 |
| BM-D07 | [beam-adopted-specification.md:15](../specs/structure/beam-adopted-specification.md) | 0 | 0 | 0 |
| BM-D08 | [beam-adopted-specification.md:16](../specs/structure/beam-adopted-specification.md) | 0 | 0 | 0 |
| BM-D09 | [beam-adopted-specification.md:17](../specs/structure/beam-adopted-specification.md) | 0 | 0 | 0 |
| CN-D01 | [connection-adopted-specification.md:9](../specs/structure/connection-adopted-specification.md) | 0 | 0 | 0 |
| CN-D02 | [connection-adopted-specification.md:10](../specs/structure/connection-adopted-specification.md) | 0 | 0 | 0 |
| CN-D03 | [connection-adopted-specification.md:11](../specs/structure/connection-adopted-specification.md) | 0 | 0 | 0 |
| CN-D04 | [connection-adopted-specification.md:12](../specs/structure/connection-adopted-specification.md) | 0 | 0 | 0 |
| CN-D05 | [connection-adopted-specification.md:13](../specs/structure/connection-adopted-specification.md) | 0 | 0 | 0 |
| CN-D06 | [connection-adopted-specification.md:14](../specs/structure/connection-adopted-specification.md) | 0 | 0 | 0 |
| CN-D07 | [connection-adopted-specification.md:15](../specs/structure/connection-adopted-specification.md) | 0 | 0 | 0 |
| CN-D08 | [connection-adopted-specification.md:16](../specs/structure/connection-adopted-specification.md) | 0 | 0 | 0 |
| CN-D09 | [connection-adopted-specification.md:17](../specs/structure/connection-adopted-specification.md) | 0 | 0 | 0 |
| EIG-01 | [eigen-extraction-contract.md:7](../specs/foundation/eigen-extraction-contract.md) | 18 | 3 | 0 |
| EIG-02 | [eigen-extraction-contract.md:28](../specs/foundation/eigen-extraction-contract.md) | 12 | 4 | 0 |
| EIG-03 | [eigen-extraction-contract.md:49](../specs/foundation/eigen-extraction-contract.md) | 20 | 5 | 1 |
| EIG-04 | [eigen-extraction-contract.md:81](../specs/foundation/eigen-extraction-contract.md) | 10 | 3 | 0 |
| EIG-05 | [eigen-extraction-contract.md:95](../specs/foundation/eigen-extraction-contract.md) | 11 | 2 | 0 |
| EIG-06 | [eigen-extraction-contract.md:117](../specs/foundation/eigen-extraction-contract.md) | 9 | 1 | 0 |
| EIG-07 | [eigen-extraction-contract.md:132](../specs/foundation/eigen-extraction-contract.md) | 7 | 3 | 0 |
| EIG-08 | [eigen-extraction-contract.md:142](../specs/foundation/eigen-extraction-contract.md) | 10 | 2 | 0 |
| EIG-09 | [eigen-extraction-contract.md:158](../specs/foundation/eigen-extraction-contract.md) | 4 | 1 | 0 |
| EIG-10 | [eigen-extraction-contract.md:166](../specs/foundation/eigen-extraction-contract.md) | 15 | 3 | 0 |
| EIG-11 | [eigen-extraction-contract.md:176](../specs/foundation/eigen-extraction-contract.md) | 16 | 7 | 0 |
| EIG-12 | [eigen-extraction-contract.md:192](../specs/foundation/eigen-extraction-contract.md) | 25 | 3 | 0 |
| EIG-13 | [eigen-extraction-contract.md:215](../specs/foundation/eigen-extraction-contract.md) | 18 | 4 | 0 |
| EIG-14 | [eigen-extraction-contract.md:235](../specs/foundation/eigen-extraction-contract.md) | 0 | 3 | 0 |
| IO4-01 | [case-schema-contract.md:93](../specs/foundation/case-schema-contract.md) | 111 | 15 | 0 |
| IO4-02 | [io-execution-adoption-specification.md:27](../specs/foundation/io-execution-adoption-specification.md) | 0 | 5 | 0 |
| IO4-03 | [io-execution-adoption-specification.md:35](../specs/foundation/io-execution-adoption-specification.md) | 1 | 1 | 0 |
| IO4-04 | [io-execution-adoption-specification.md:54](../specs/foundation/io-execution-adoption-specification.md) | 0 | 1 | 0 |
| IO4-05 | [io-execution-adoption-specification.md:66](../specs/foundation/io-execution-adoption-specification.md) | 0 | 1 | 0 |
| IO4-06 | [io-execution-adoption-specification.md:72](../specs/foundation/io-execution-adoption-specification.md) | 0 | 1 | 0 |
| IO4-07 | [io-execution-adoption-specification.md:82](../specs/foundation/io-execution-adoption-specification.md) | 0 | 1 | 0 |
| J1-D01 | [j1-discrete-equilibrium-and-chart-contract.md:5](../specs/coupled/j1-discrete-equilibrium-and-chart-contract.md) | 0 | 0 | 0 |
| J1-D02 | [j1-discrete-equilibrium-and-chart-contract.md:10](../specs/coupled/j1-discrete-equilibrium-and-chart-contract.md) | 0 | 0 | 0 |
| J1-D03 | [j1-discrete-equilibrium-and-chart-contract.md:11](../specs/coupled/j1-discrete-equilibrium-and-chart-contract.md) | 0 | 0 | 0 |
| J1-D04 | [j1-discrete-equilibrium-and-chart-contract.md:12](../specs/coupled/j1-discrete-equilibrium-and-chart-contract.md) | 0 | 0 | 0 |
| J1-D05 | [j1-discrete-equilibrium-and-chart-contract.md:13](../specs/coupled/j1-discrete-equilibrium-and-chart-contract.md) | 0 | 0 | 0 |
| J1-D06 | [j1-discrete-equilibrium-and-chart-contract.md:14](../specs/coupled/j1-discrete-equilibrium-and-chart-contract.md) | 0 | 0 | 0 |
| J1-D07 | [j1-discrete-equilibrium-and-chart-contract.md:15](../specs/coupled/j1-discrete-equilibrium-and-chart-contract.md) | 0 | 0 | 0 |
| J1-D08 | [j1-discrete-equilibrium-and-chart-contract.md:16](../specs/coupled/j1-discrete-equilibrium-and-chart-contract.md) | 0 | 0 | 0 |
| J1-D09 | [j1-discrete-equilibrium-and-chart-contract.md:17](../specs/coupled/j1-discrete-equilibrium-and-chart-contract.md) | 0 | 0 | 0 |
| MAT-D01 | [material-initial-strain-contract.md:129](../specs/structure/material-initial-strain-contract.md) | 0 | 0 | 0 |
| MAT-D02 | [material-initial-strain-contract.md:129](../specs/structure/material-initial-strain-contract.md) | 0 | 0 | 0 |
| MAT-D03 | [material-initial-strain-contract.md:129](../specs/structure/material-initial-strain-contract.md) | 0 | 0 | 0 |
| MAT-D04 | [material-layer-adoption-contract.md:9](../specs/structure/material-layer-adoption-contract.md) | 0 | 0 | 0 |
| MAT-D05 | [material-layer-adoption-contract.md:10](../specs/structure/material-layer-adoption-contract.md) | 0 | 0 | 0 |
| MAT-D06 | [material-layer-adoption-contract.md:11](../specs/structure/material-layer-adoption-contract.md) | 0 | 0 | 0 |
| MAT-D07 | [material-layer-adoption-contract.md:12](../specs/structure/material-layer-adoption-contract.md) | 0 | 0 | 0 |
| MAT-D08 | [material-shell-adopted-specification.md:9](../specs/structure/material-shell-adopted-specification.md) | 0 | 0 | 0 |
| MAT-D09 | [material-shell-adopted-specification.md:10](../specs/structure/material-shell-adopted-specification.md) | 0 | 0 | 0 |
| MAT-D10 | [material-shell-adopted-specification.md:11](../specs/structure/material-shell-adopted-specification.md) | 0 | 0 | 0 |
| MAT-D11 | [material-shell-adopted-specification.md:12](../specs/structure/material-shell-adopted-specification.md) | 0 | 0 | 0 |
| MAT-D12 | [material-shell-adopted-specification.md:13](../specs/structure/material-shell-adopted-specification.md) | 0 | 0 | 0 |
| MAT-D13 | [material-shell-adopted-specification.md:14](../specs/structure/material-shell-adopted-specification.md) | 0 | 0 | 0 |
| N4-01 | [numerical-foundation-adoption-specification.md:6](../specs/foundation/numerical-foundation-adoption-specification.md) | 14 | 4 | 0 |
| N4-02 | [numerical-foundation-adoption-specification.md:19](../specs/foundation/numerical-foundation-adoption-specification.md) | 162 | 41 | 4 |
| N4-03 | [numerical-foundation-adoption-specification.md:37](../specs/foundation/numerical-foundation-adoption-specification.md) | 99 | 39 | 2 |
| N4-04 | [numerical-foundation-adoption-specification.md:47](../specs/foundation/numerical-foundation-adoption-specification.md) | 0 | 0 | 0 |
| N4-05 | [numerical-foundation-adoption-specification.md:57](../specs/foundation/numerical-foundation-adoption-specification.md) | 104 | 17 | 0 |
| N4-06 | [numerical-foundation-adoption-specification.md:75](../specs/foundation/numerical-foundation-adoption-specification.md) | 157 | 33 | 0 |
| N4-07 | [numerical-foundation-adoption-specification.md:85](../specs/foundation/numerical-foundation-adoption-specification.md) | 98 | 7 | 1 |
| N4-08 | [numerical-foundation-adoption-specification.md:100](../specs/foundation/numerical-foundation-adoption-specification.md) | 121 | 42 | 1 |
| N4-09 | [numerical-foundation-adoption-specification.md:120](../specs/foundation/numerical-foundation-adoption-specification.md) | 65 | 33 | 1 |
| N4-10 | [numerical-foundation-adoption-specification.md:130](../specs/foundation/numerical-foundation-adoption-specification.md) | 11 | 7 | 0 |
| QA-01 | [quality-plan.md:13](../quality/quality-plan.md) | 0 | 0 | 0 |
| QA-02 | [quality-plan.md:23](../quality/quality-plan.md) | 0 | 0 | 0 |
| QA-03 | [quality-plan.md:31](../quality/quality-plan.md) | 0 | 0 | 0 |
| QA-04 | [quality-plan.md:39](../quality/quality-plan.md) | 0 | 0 | 0 |
| QA-05 | [quality-plan.md:51](../quality/quality-plan.md) | 0 | 0 | 0 |
| QA-06 | [quality-plan.md:65](../quality/quality-plan.md) | 0 | 0 | 0 |
| QA-07 | [quality-plan.md:71](../quality/quality-plan.md) | 0 | 0 | 1 |
| QA-08 | [quality-plan.md:82](../quality/quality-plan.md) | 0 | 0 | 0 |
| QA-09 | [quality-plan.md:88](../quality/quality-plan.md) | 0 | 0 | 0 |
| QA-10 | [quality-plan.md:100](../quality/quality-plan.md) | 0 | 0 | 0 |
| SH-D01 | [main-shell-adopted-specification.md:11](../specs/structure/main-shell-adopted-specification.md) | 5 | 1 | 0 |
| SH-D02 | [main-shell-adopted-specification.md:12](../specs/structure/main-shell-adopted-specification.md) | 0 | 0 | 0 |
| SH-D03 | [main-shell-adopted-specification.md:13](../specs/structure/main-shell-adopted-specification.md) | 32 | 20 | 0 |
| SH-D04 | [main-shell-adopted-specification.md:14](../specs/structure/main-shell-adopted-specification.md) | 0 | 0 | 0 |
| SH-D05 | [main-shell-adopted-specification.md:15](../specs/structure/main-shell-adopted-specification.md) | 5 | 3 | 0 |
| SH-D06 | [main-shell-adopted-specification.md:16](../specs/structure/main-shell-adopted-specification.md) | 0 | 0 | 0 |
| SH-D07 | [main-shell-adopted-specification.md:17](../specs/structure/main-shell-adopted-specification.md) | 20 | 8 | 0 |
| SH-D08 | [main-shell-adopted-specification.md:18](../specs/structure/main-shell-adopted-specification.md) | 0 | 0 | 0 |
| SI-D01 | [structural-adopted-specification.md:9](../specs/structure/structural-adopted-specification.md) | 0 | 0 | 0 |
| SI-D02 | [structural-adopted-specification.md:10](../specs/structure/structural-adopted-specification.md) | 0 | 0 | 0 |
| SI-D03 | [structural-adopted-specification.md:11](../specs/structure/structural-adopted-specification.md) | 0 | 0 | 0 |
| SI-D04 | [structural-adopted-specification.md:12](../specs/structure/structural-adopted-specification.md) | 0 | 0 | 0 |
| SI-D05 | [structural-adopted-specification.md:13](../specs/structure/structural-adopted-specification.md) | 0 | 0 | 0 |
| SI-D06 | [structural-adopted-specification.md:14](../specs/structure/structural-adopted-specification.md) | 0 | 0 | 0 |
| SI-D07 | [structural-adopted-specification.md:15](../specs/structure/structural-adopted-specification.md) | 0 | 0 | 0 |
| SI-D08 | [structural-adopted-specification.md:16](../specs/structure/structural-adopted-specification.md) | 0 | 0 | 0 |
| SI-D09 | [structural-adopted-specification.md:17](../specs/structure/structural-adopted-specification.md) | 0 | 0 | 0 |
| ST-01 | [design-decisions.md:39](../specs/design-decisions.md) | 0 | 0 | 0 |
| ST-02 | [design-decisions.md:40](../specs/design-decisions.md) | 0 | 0 | 0 |
| ST-03 | [design-decisions.md:41](../specs/design-decisions.md) | 0 | 0 | 0 |
| ST-04 | [design-decisions.md:42](../specs/design-decisions.md) | 0 | 0 | 0 |
| ST-05 | [design-decisions.md:43](../specs/design-decisions.md) | 0 | 0 | 0 |
| ST-06 | [design-decisions.md:44](../specs/design-decisions.md) | 0 | 0 | 0 |
| ST-07 | [design-decisions.md:45](../specs/design-decisions.md) | 0 | 0 | 0 |
| ST-08 | [design-decisions.md:46](../specs/design-decisions.md) | 0 | 0 | 0 |
| ST-09 | [design-decisions.md:47](../specs/design-decisions.md) | 0 | 0 | 0 |
| ST-10 | [design-decisions.md:48](../specs/design-decisions.md) | 0 | 0 | 0 |
| UI-D01 | [uvlm-induction-adopted-specification.md:9](../specs/aero/uvlm-induction-adopted-specification.md) | 17 | 8 | 0 |
| UI-D02 | [uvlm-induction-adopted-specification.md:10](../specs/aero/uvlm-induction-adopted-specification.md) | 10 | 2 | 0 |
| UI-D03 | [uvlm-induction-adopted-specification.md:11](../specs/aero/uvlm-induction-adopted-specification.md) | 27 | 10 | 0 |
| UI-D04 | [uvlm-induction-adopted-specification.md:12](../specs/aero/uvlm-induction-adopted-specification.md) | 16 | 3 | 0 |
| UI-D05 | [uvlm-induction-adopted-specification.md:13](../specs/aero/uvlm-induction-adopted-specification.md) | 25 | 10 | 0 |
| UI-D06 | [uvlm-induction-adopted-specification.md:14](../specs/aero/uvlm-induction-adopted-specification.md) | 53 | 10 | 0 |
| UW-D01 | [uvlm-wake-adopted-specification.md:9](../specs/aero/uvlm-wake-adopted-specification.md) | 0 | 0 | 0 |
| UW-D02 | [uvlm-wake-adopted-specification.md:10](../specs/aero/uvlm-wake-adopted-specification.md) | 0 | 0 | 0 |
| UW-D03 | [uvlm-wake-adopted-specification.md:11](../specs/aero/uvlm-wake-adopted-specification.md) | 0 | 0 | 0 |
| UW-D04 | [uvlm-wake-adopted-specification.md:12](../specs/aero/uvlm-wake-adopted-specification.md) | 3 | 2 | 0 |
| UW-D05 | [uvlm-wake-adopted-specification.md:13](../specs/aero/uvlm-wake-adopted-specification.md) | 0 | 0 | 0 |
| UW-D06 | [uvlm-wake-adopted-specification.md:14](../specs/aero/uvlm-wake-adopted-specification.md) | 0 | 0 | 0 |
| UW-D07 | [uvlm-wake-adopted-specification.md:15](../specs/aero/uvlm-wake-adopted-specification.md) | 17 | 5 | 0 |
| UW-D08 | [uvlm-wake-adopted-specification.md:16](../specs/aero/uvlm-wake-adopted-specification.md) | 3 | 2 | 0 |

## 受入計画とケース

| ID | 定義 | 試験 | 記録 | コード |
|---|---|---:|---:|---:|
| A5-01 | [practical-acceptance-plan.md:5](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| A5-02 | [practical-acceptance-plan.md:24](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| A5-03 | [practical-acceptance-plan.md:45](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| A5-04 | [practical-acceptance-plan.md:65](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| A5-05 | [practical-acceptance-plan.md:87](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| A5-06 | [practical-acceptance-plan.md:97](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| A5-07 | [practical-acceptance-plan.md:111](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| P-CURVED-BOX | [practical-acceptance-plan.md:65](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| PERF-FL-01 | [case-status.md:29](../quality/case-status.md) | 0 | 0 | 0 |
| V-AERO | [practical-acceptance-plan.md:36](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| V-BEAM | [practical-acceptance-plan.md:32](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| V-J1 | [practical-acceptance-plan.md:38](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| V-JOINT | [practical-acceptance-plan.md:33](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| V-LOAD | [practical-acceptance-plan.md:34](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| V-MAT | [practical-acceptance-plan.md:31](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| V-MULTI | [practical-acceptance-plan.md:37](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| V-PG | [practical-acceptance-plan.md:30](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| V-RECOVER | [practical-acceptance-plan.md:35](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| V-SEARCH | [practical-acceptance-plan.md:39](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| V-STATE | [practical-acceptance-plan.md:40](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| V-USE | [practical-acceptance-plan.md:41](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| X-S10-DECAY | [practical-acceptance-plan.md:55](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| X-S10-FLUTTER | [practical-acceptance-plan.md:54](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| X-S10-GVT-LE | [practical-acceptance-plan.md:51](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| X-S10-GVT-TE | [practical-acceptance-plan.md:52](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |
| X-S10-STATIC | [practical-acceptance-plan.md:53](../requirements/practical-acceptance-plan.md) | 0 | 0 | 0 |

## 実装パッケージ

| ID | 定義 | 試験 | 記録 | コード |
|---|---|---:|---:|---:|
| IM-00 | [design-decisions.md:4](../specs/design-decisions.md) | 0 | 0 | 0 |
| IM-01 | [implementation-start-review.md:40](../plans/implementation-start-review.md) | 0 | 0 | 0 |
| IM-02 | [implementation-start-review.md:41](../plans/implementation-start-review.md) | 0 | 0 | 0 |
| IM-03 | [implementation-start-review.md:42](../plans/implementation-start-review.md) | 0 | 0 | 0 |
| IM-04 | [implementation-start-review.md:43](../plans/implementation-start-review.md) | 0 | 0 | 0 |
| IM-05 | [implementation-start-review.md:44](../plans/implementation-start-review.md) | 0 | 0 | 0 |
| IM-06 | [implementation-start-review.md:45](../plans/implementation-start-review.md) | 0 | 0 | 0 |
| IM-07 | [implementation-start-review.md:46](../plans/implementation-start-review.md) | 0 | 0 | 0 |
| IM-08 | [implementation-start-review.md:47](../plans/implementation-start-review.md) | 0 | 0 | 0 |

## 残リスク

| ID | 定義 | 試験 | 記録 | コード |
|---|---|---:|---:|---:|
| R-01 | [implementation-start-review.md:69](../plans/implementation-start-review.md) | 0 | 0 | 0 |
| R-02 | [implementation-start-review.md:70](../plans/implementation-start-review.md) | 1 | 1 | 0 |
| R-03 | [implementation-start-review.md:71](../plans/implementation-start-review.md) | 0 | 0 | 0 |
| R-04 | [implementation-start-review.md:72](../plans/implementation-start-review.md) | 0 | 0 | 0 |
| R-05 | [implementation-start-review.md:73](../plans/implementation-start-review.md) | 0 | 0 | 0 |
| R-06 | [implementation-start-review.md:74](../plans/implementation-start-review.md) | 0 | 0 | 0 |
| R-07 | [implementation-start-review.md:75](../plans/implementation-start-review.md) | 0 | 0 | 0 |
