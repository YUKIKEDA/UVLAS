# 構成とモジュール

リポジトリの配置の正本。モジュールやディレクトリを足すPRが、同じPRでこのページの行を足す。

## ディレクトリ

| パス | 内容 |
|---|---|
| `include/uvlas/<module>/` | 製品の公開ヘッダ。現状は内部APIで、互換性は固定していない |
| `src/<module>/` | 製品の実装。`src/app/` はCLI |
| `tests/` | 製品の試験（Catch2、CTest）。`tests/<module>/` に分ける |
| `tests/support/` | 確保経路の監査とthread数別のバイト比較を行う試験専用の枠 |
| `benchmarks/` | Google Benchmarkによる製品kernelの計測。製品やCTestへリンクしない |
| `cmake/` | コンパイラ方針、共通オプション、依存、vcpkgのoverlay port |
| `docs/requirements/` | 製品要求、受入精度、受入計画、利用者決定の履歴 |
| `docs/specs/<area>/` | 採用仕様。実装の正本 |
| `docs/plans/` | 開始審査（実装パッケージ、順序、残リスク） |
| `docs/engineering/` | 共通Rules、知見と根拠、レビュー・プロファイル |
| `docs/adr/` | 決定の記録。1決定1ファイル |
| `docs/quality/` | 品質保証計画、追跡表、受入ケースの状態、結果変更記録 |
| `qa/baselines/<case>/` | 基準結果。プラットフォームごとに1ファイル。`tools/qa/results.py` だけが書く |
| `qa/<case>/` | 基準ケースのプログラム |
| `docs/papers/` | 技術報告書 |
| `docs/verification/<IM-ID>/` | 製品の検証記録。最初の記録を作るPRで作る |
| `rd/<domain>/<family>/` | R&Dの記録。系統ごとに設計メモ、`records/`、`experiments/`、`fixtures/`。`rd/planning/` は実装前の計画と台帳 |
| `rd/evidence-index.json` | 記録がハッシュで引用する大きい証拠ファイルの索引。実体はGitHubのリリースのアセット |
| `tools/rd/` | R&Dツリーの再配置と検査 |
| `tools/qa/` | 追跡（`trace.py`）、基準結果と結果の変更管理（`results.py`） |
| `tools/dev/` | 開発用の補助（MSVC環境の起動、Skillの複製、CIの変更範囲判定と集約の検査） |
| `.agents/skills/`、`.claude/skills/` | Skillの正本と、その複製 |
| `.cursor/rules/` | 規則 |
| `ruff.toml`、`pyrightconfig.json`、`.clang-format`、`.clang-tidy` | 静的検査の設定。どれもゲートが使う |
| `out/` | ビルド生成物、vcpkgの導入物。git対象外 |
| `build/` | R&Dのスクリプトを走らせるための一時領域。git対象外。消してよい |
| `.dev/` | 一時置き場。git対象外 |

## 製品モジュール

依存の向きは上から下。空力は構造の内部を知らず、構造は渦や後流の内部を知らない。連成が両者を組み合わせる。入出力からコアを呼び、要素や誘導速度の関数からファイル出力やCLIへ依存させない。

| モジュール | 責任 | 導入するパッケージ |
|---|---|---|
| `core` | ビルド情報、強いID、checked寸法、ResourcePlan、`Failure`、`ErrorCertificate`、状態・資格・cache・Workspace、工程時間・メモリ計測 | IM00-3の[失敗・認定の型](specs/foundation/failure-certificate-contract.md)、IM00-2の[寸法・資源契約](specs/foundation/dimensions-resource-contract.md)、IM00-9の[状態契約](specs/foundation/state-contract.md)、[計測契約](specs/foundation/instrumentation-contract.md) |
| `numerics` | 密・疎の保持、PointSoA、DofLayout、精度保持・補償原残差・明示投影、疎LU/正定値資格付き疎LLT。密基準解法はIM00-6。固有抽出の共通型・原残差／Gram認定・作用予算はIM00-11 | IM00-4の[保持・境界契約](specs/foundation/numerical-storage-contract.md)、IM00-5の[精度契約](specs/foundation/precision-contract.md)、IM00-7の[疎求解契約](specs/foundation/sparse-solve-contract.md)、[固有抽出契約](specs/foundation/eigen-extraction-contract.md)。一般密行列の完全スペクトル、選択pair、独立した左右clusterは `eigen_full.hpp`。所有packetは固定尺度・元係数・精度世代・認定を保持する。原M資格付き対称一般化全スペクトルは [eigen_generalized.hpp](../include/uvlas/numerics/eigen_generalized.hpp)。一般作用Arnoldiは [eigen_arnoldi.hpp](../include/uvlas/numerics/eigen_arnoldi.hpp) の所有結果・Workspace・元作用再認定を含む。検証範囲は対応する実装Issueと `docs/verification/IM-00` の記録で確認する |
| `structure` | PG参照準備用のexact有理数・二変数多項式、rational graph・半セル・矩形活動map、全原辞書・full Vと必須場のexact構成map。材料、梁、接続、荷重、回復は後続 | IM-01、IM-02。[exact代数の資源契約](specs/structure/pg-exact-algebra-resource-contract.md)と[graph・cell契約](specs/structure/pg-graph-cell-contract.md)に従う準備処理を `exact_polynomial.hpp` と `pg_graph.hpp` が公開する。[全原辞書の契約](specs/structure/pg-full-dictionary-contract.md)と `pg_coefficients.hpp` はanchor補完と原列由来を保持する。活動rank・条件付け後の物理基底・求積・物理系の完成を表さない。多角形と全積分分割の実装開始条件は[分割契約](specs/structure/pg-integration-partition-contract.md)を参照する |
| `aero` | [有限三角形の厳密述語](specs/aero/exact-triangle-predicates-contract.md)と[所有形状・厳密query位置式](specs/aero/owned-geometry-query-contract.md)を準備部品として実装。元の接続・alias・世代を保持する。[局所packet](specs/aero/local-geometry-packet-contract.md)が導出量の包絡とIDに基づく局所incidence・接触政策を扱う。[P1核の追加契約](specs/aero/p1-induction-kernel-contract.md)に従い、局所P1基底・全13変数微分・丸め包絡・固定遠方則の照合を実装した。大域幾何、TE所有、後流、荷重、物理全体系の実務受入は後続 | IM-03 |
| `coupling` | 初期化、強連成、確定、安定性と探索 | IM-04、IM-05 |
| `io` | 復号済みASCII IDとnative整数の厳格変換、元byte列とSHA-256・木・raw tokenを所有する入力文書。ケースのnested schema・物理Model、外部格子、結果、checkpointは後続 | IM06-1の[厳格scalar契約](specs/foundation/strict-input-scalars-contract.md)と [input_scalar.hpp](../include/uvlas/io/input_scalar.hpp)、IM06-2の[入力文書・資源契約](specs/foundation/input-document-resource-contract.md)と [input_document.hpp](../include/uvlas/io/input_document.hpp)。文書全体を検査してroot envelopeを確認するが、物理Modelは生成しない。全nested schemaの設計は[全case契約](specs/foundation/case-schema-contract.md)と[フィールド辞書](specs/foundation/case-schema-fields.md)に定める。FD-08の初期摂動には[完全な物理初期状態の追加規約](specs/foundation/case-initial-physical-state-contract.md)を優先する。schema validatorの実装完了を示すものではない。[型付きschemaの実装契約](specs/foundation/typed-case-schema-implementation-contract.md)が所有型・資源・全検査の開始条件を具体化する |
| `python` | 任意機能のCPython拡張。所有する数値保持と状態protocol、NumPy境界 | IM00-10の[Python境界契約](specs/foundation/python-binding-contract.md) |
| `app` | CLI | 現在は `--help` と `--version`。解析コマンドはIM-06 |

表のモジュールのうち、いまディレクトリがあるのは `core`、`numerics`、`structure`、`aero`、`io`、`app`、`python`。他は、最初のソースを足すPRで作る。

## 名前

- CMakeターゲット：`uvlas_<module>`。試験は `uvlas_add_test(<module>.<topic> …)`。
- 名前空間：`uvlas::<module>`。
- ロードマップの行ID：基盤は `F-n`、実装パッケージは `IM00-n` 〜 `IM08-n`。
- 採用仕様の決定ID（SH-D、ST、AE、J1-D、N4、IO4、A5）と受入ケースID（V-*、X-*、P-*）は、仕様の表記のまま使う。
