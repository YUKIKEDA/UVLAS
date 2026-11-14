# 構造モデルの採用審査

<!-- active-uvlm-status:start -->
2026-11-06・最新：**実装前1〜4の採用設計、5の受入計画、6の最初の実装パッケージ審査を完了。次はIM-00：数値表現・求解・状態契約の製品実装。** 両OSの数値／Python接続を検証し、19ケース・37製品要求を対応付け、1072ソース・687成果物を監査。**S10実験参照の認定、実PG全規模・全自由後流の製品性能、全実務受入・リリースは未完了。本体コードは未着手。**

[開始判定・実装順と残リスク](../../docs/plans/implementation-start-review.md)、[受入計画](../../docs/requirements/practical-acceptance-plan.md)、[ケース台帳](../acceptance/plan/records/practical-acceptance-case-catalog.json)。
<!-- active-uvlm-status:end -->

2026-10-31に開始した構造モデルの採用審査。**2026-11-05に1.1〜1.5を統合して終了。** 現在の採用式・適用条件は[構造統合仕様](../../docs/specs/structure/structural-adopted-specification.md)、根拠は[終了判断](../structure/integration/records/structural-adoption-results.md)を参照。以下の初期候補・未達の記録は履歴として保持する。

今回の到達点：[歪み格子の非一様薄肉曲げ](../structure/main-shell/q9-plate/records/q9-distorted-plate-results.md) を事前仕様化して8ケースを実施し、別経路の観測・成果物監査を完了した。数値処理の必須条件と最細精度・最後の次数は成立したが、薄肉一様性の基準は未達。補正なし完全積分Q9の採用保留を維持する。[無せん断適合性の解析反例](../structure/main-shell/q9-plate/q9-distorted-shear-compatibility.md) も整理し、次に閉じる対象をST-04の対策成分・具体方式とした。構造審査全体の終了ではない。

2026-11-06現在：実装前1〜4の採用設計、5の受入計画、6の最初の依存閉包審査を終了。IM-00の製品実装に着手可能。本体は未着手で、S10実験認定・製品規模性能・全実務受入は未完了。

## 既存成果から維持する判断

- 補正なし3833を主シェルとして採用しない。台形の純曲げ補間で人工せん断が生じ、求積点の増加だけでは解消しなかった。[比較経緯](../structure/main-shell/q9-plate/q9-thickness-quadratic-design.md) を採用審査の起点とする。初期の3423／3833調査順を最新の採用順位と読み違えない。
- Q9・厚さ二次・81座標は比較基盤として保持。現在は共有Q3・厚さ二次・144座標も比較し、第一段の表現性が成立した。位置・厚さ微分の意味を保持し、純曲げの厳密再現だけで採用しない。
- 同じ3D補間空間のHex27との座標変換・独立組立照合は、質量・接線・境界仮想仕事の根拠として利用する。薄肉構造としての適格性は別に確認する。
- 薄肉時には求積、組立、残差評価、解・補正の保持、返値の精度を分離する。初回の未達を後続の高精度試験で上書きせず、[精度診断](../structure/main-shell/q9-plate/records/q9-precision-diagnostic-experiment.md) と合わせて読む。Decimalを製品依存として採用したことにはしない。
- 同じ参照配置・材料軸・指定初期ひずみ・エネルギーに基づく内力と一貫接線を用い、初期応力の幾何接線を保持する。限定的な材料・予応力試験の成功を全ての材料・変形へ一般化しない。

## 採用を閉じるための統合表

| 対象 | 既存の証拠と契約 | 残る採用判断・証拠 | 完了時の成果物 |
|---|---|---|---|
| ST-01/04 主シェル・求積・安定性 | [PG-ANCF採用仕様](../../docs/specs/structure/main-shell-adopted-specification.md)、[R4整合](../structure/main-shell/paired-geometry/records/paired-reference-r4-results.md)、[R5採用](../structure/main-shell/paired-geometry/records/paired-shell-r5-results.md) | **1.1の採用判断は完了**。製品全実務受入・具体数値依存は後続 | 採用式、必要微分、適用／拒否、全トレース、求積と保存証拠を固定済み |
| ST-02 参照配置・幾何 | [主シェル仕様](../../docs/specs/structure/main-shell-adopted-specification.md)、[構造統合仕様](../../docs/specs/structure/structural-adopted-specification.md)、接合参照面・後退テーパ | **構造モデル採用は完了**。外部格子の具体変換／誤差予算は4／5 | 参照写像・物理厚さ／法線・正則性・面対応を固定済み |
| ST-03/CV-03 材料・積層・初期ひずみ | [材料採用仕様](../../docs/specs/structure/material-shell-adopted-specification.md)、[1.2検証・採否](../structure/material/records/material-shell-adoption-results.md) | **1.2の数値仕様・主シェル結合は完了**。梁は別の一般化ひずみ／断面則へ縮約。実材料の用途別受入は5 | 全3D入力・等価則、層／軸／初期ひずみ、全接線・質量・両側回復の境界を固定済み |
| ST-05/07 質量・偏心慣性 | [接続仕様](../../docs/specs/structure/connection-adopted-specification.md)、[同一平衡のM/K/C](../structure/integration/records/structural-combined-result.json)、質量台帳 | **構造モデル採用は完了**。動的高階微分・具体型／性能は3／4 | 同じ運動学のM・慣性・重力・所有を統合済み |
| ST-06 梁 | [SE(3)梁採用仕様](../../docs/specs/structure/beam-adopted-specification.md)、[1.3検証・採否](../structure/beam/records/beam-adoption-results.md) | **1.3の採用判断は完了**。実接合は1.4、動的方式／具体精度は3／4、実務範囲は5 | 軸・せん断・ねじり・曲げ、断面／反り境界、DOF・接線・質量、シェル受渡しを固定済み |
| ST-08 界面・支持・剛体結合 | [接続採用仕様](../../docs/specs/structure/connection-adopted-specification.md)、[1.4の実平衡・モード・有限面・非アフィン比較](../structure/connection/records/connection-adoption-results.md) | **1.4の採用判断は完了**。任意接合の収束・局所応力受入は5、動的拘束は3 | C/J/H・時間微分、独立化と拒否、原反力、面内／厚さ収束条件を固定済み |
| ST-09 荷重・支持・物理減衰 | [統合仕様](../../docs/specs/structure/structural-adopted-specification.md)、[荷重／減衰](../structure/integration/records/structural-load-damping-result.json)、[時系列と独立支持参照](../structure/integration/records/structural-history-state-result.json) | **1.5採用完了**。時間ステージと全空力依存は2／3、係数の実材料同定は5 | 全荷重接線、客観粘性、動く支持の仕事、イベント／片側／外挿拒否を固定済み |
| ST-10 内力・接線・結果回復 | [統合仕様](../../docs/specs/structure/structural-adopted-specification.md)、[応力／切断面](../structure/integration/records/structural-recovery-result.json)、[全残差の接線](../structure/integration/records/structural-combined-result.json) | **1.5採用完了**。製品生成方式・具体型は4、用途別回復受入は5 | 同じ状態・測度・軸・層側・切断面・基準点と収支、原反力を固定済み |

## 審査の進行順と依存

1. ST-01〜04について既存証拠の適用範囲を整理し、主シェルの最終採用に不足する比較ケースを、薄さ・歪み・曲率・積層・荷重・評価量・独立参照・許容差・停止条件まで具体化する。既存試験と同じケースを無条件に再実行しない。
2. その結果と並行して梁候補を審査し、シェルと梁が共有する位置・方向・仕事の境界を閉じる。未選定の梁を剛体線やシェルの細帯で代替した扱いにしない。
3. 共通の運動学から質量・偏心慣性、接合・支持、減衰・荷重を導く。接合を自由度の単純共有で済ませず、厚さ座標変換と物理的な接続条件を区別する。
4. 内力・接線・出力回復を統合し、採用対象ごとの不変条件・検証根拠と台帳IDを対応させる。未達は未達として残し、成立した候補だけを採用仕様へ昇格させる。

後続の2.UVLM、3.J1、4.数値処理・I/F、5.実務ケース・受入根拠、6.着手審査との接続条件も記録する。構造の質量・接線・荷重・支持・減衰が確定しても、時間積分や連成・線形安定性の採用を単独で決定したことにはしない。

## この段階の終了条件

上表の対象に採用式・適用範囲・必要な接続契約・独立検証が対応し、ST-01〜10とCV-03の採用状態が更新されていること。製品I/Fへ渡す量・単位・微分・反力・履歴を列挙し、J1と数値依存へ残す共同決定を明示する。未解決事項を実装者の裁量で補う仕様にはしない。

2026-11-05：**2.2採用完了、現在地2.3。** [後流・全荷重仕様](../../docs/specs/aero/uvlm-wake-adopted-specification.md)、[16群の判定](../aero/wake-loads/records/uvlm-wake-adoption-results.md)。製品の時間方式・UVLM全体・本体実装開始は未完了。

2026-11-05・2.3更新：[独立参照の実数値・未達・次の終了条件](../aero/unsteady-reference/records/uvlm-independent-reference-progress.md)。定常の総量・感度は限定成立、非定常の空間収束を追加中。UVLM全体と本体実装開始は未完了。

2026-11-05・第1参照キャンペーン終了：**2.2完了、現在地2.3。** 定常2翼・20帯、非定常上下運動、主翼／尾翼の独立比較が限定成立。[成立範囲・未達・残条件・監査](../aero/unsteady-reference/records/uvlm-reference-campaign-summary-2026-11-05.md)。2.3全体・2.4・本体実装開始は未完了。
