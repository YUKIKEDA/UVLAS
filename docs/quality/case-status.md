# 受入ケースの状態

受入計画（[A5](../requirements/practical-acceptance-plan.md)）の19ケースの、現在の状態。状態を進めるPRが、この表の行を書き換え、根拠の検証記録を名指しする（[QA-02](quality-plan.md)）。根拠の記録のない行を進めない。

状態は、種類の違う証拠を別々に持つ：`Specified`（入力・参照・評価・閾値を定めた）、`SourceAvailable`（資料を入手した）、`ReferenceQualified`（参照を比較に使えると認定した）、`ImplementationVerified`（製品の実装が仕様に適合する）、`NumericallyQualified`（数値予算を満たす）、`ExperimentQualified`（実験参照と不確かさを認定した）、`Accepted`（受入）。一つの合格へ潰さない。

2026-11-11：製品のソルバーは未実装で、`ImplementationVerified` 以上のケースは無い。

| ケース | 対象 | 状態 | 根拠の記録 |
|---|---|---|---|
| V-PG | 主シェル：曲面、厚さ、歪み、支持、独立荷重 | Specified | なし |
| V-MAT | 材料、積層、初期ひずみ、層別回復 | Specified | なし |
| V-BEAM | 梁、偏心質量と慣性 | Specified | なし |
| V-JOINT | 接続、支持、全厚さトレース、反力 | Specified | なし |
| V-LOAD | 既知入力、指定運動、突風、支持の仕事 | Specified | なし |
| V-RECOVER | 応力、断面力、層の両側の回復 | Specified | なし |
| V-AERO | UVLM：有限核、20帯、合力、抗力 | Specified | なし |
| V-MULTI | 複数揚力面、共有とTE所有 | Specified | なし |
| V-J1 | 強連成、全未知量の線形化、初期境界微分 | Specified | なし |
| V-SEARCH | 静的継続、mode追跡、速度探索 | Specified | なし |
| V-STATE | 再開：同一環境のビット一致、破損と非互換の拒否 | Specified | なし |
| V-USE | 入力経路、CLIとPython、配布物の導入 | Specified | なし |
| X-S10-GVT-LE | 公開実験S10：地上振動試験（前縁側） | SourceAvailable | [資料の診断](../../rd/acceptance/s10/records/s10-reference-metadata-result.json)。不確かさと校正履歴は未認定 |
| X-S10-GVT-TE | 公開実験S10：地上振動試験（後縁側） | SourceAvailable | 同上 |
| X-S10-STATIC | 公開実験S10：静的変形 | SourceAvailable | 同上。支持座標、marker番号、未変形の基準が未認定 |
| X-S10-FLUTTER | 公開実験S10：フラッタ | SourceAvailable | 同上 |
| X-S10-DECAY | 公開実験S10：解放応答 | Specified | 原ファイルと構成が未選定 |
| P-CURVED-BOX | 実務統合入力：後退テーパ曲面の箱構造と尾翼 | Specified | なし |
| PERF-FL-01 | 性能受入：P-CURVED-BOX、約10万自由度＋5,000パネル、48GiB／12時間。条件はA5-06が本文で定める | Specified | なし |
