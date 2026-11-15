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
| `numerics` | 密・疎の保持、PointSoA、DofLayout、精度保持・補償原残差・明示投影、疎LU/正定値資格付き疎LLT。密基準解法はIM00-6 | IM00-4の[保持・境界契約](specs/foundation/numerical-storage-contract.md)、IM00-5の[精度契約](specs/foundation/precision-contract.md)、IM00-7の[疎求解契約](specs/foundation/sparse-solve-contract.md) |
| `structure` | PG参照係数、材料、梁、接続、荷重、回復 | IM-01、IM-02 |
| `aero` | 渦核、後流、荷重、TE所有 | IM-03 |
| `coupling` | 初期化、強連成、確定、安定性と探索 | IM-04、IM-05 |
| `io` | ケース、外部格子、結果、checkpoint | IM-06 |
| `python` | 任意機能のCPython拡張。所有する数値保持と状態protocol、NumPy境界 | IM00-10の[Python境界契約](specs/foundation/python-binding-contract.md) |
| `app` | CLI | 現在は `--help` と `--version`。解析コマンドはIM-06 |

表のモジュールのうち、いまディレクトリがあるのは `core`、`numerics`、`app`、`python`。他は、最初のソースを足すPRで作る。

## 名前

- CMakeターゲット：`uvlas_<module>`。試験は `uvlas_add_test(<module>.<topic> …)`。
- 名前空間：`uvlas::<module>`。
- ロードマップの行ID：基盤は `F-n`、実装パッケージは `IM00-n` 〜 `IM08-n`。
- 採用仕様の決定ID（SH-D、ST、AE、J1-D、N4、IO4、A5）と受入ケースID（V-*、X-*、P-*）は、仕様の表記のまま使う。
