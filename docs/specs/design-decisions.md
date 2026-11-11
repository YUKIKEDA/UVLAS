# 実装前の設計決定台帳

<!-- active-uvlm-status:start -->
2026-11-06・最新：**実装前1〜4の採用設計、5の受入計画、6の最初の実装パッケージ審査を完了。次はIM-00：数値表現・求解・状態契約の製品実装。** 両OSの数値／Python接続を検証し、19ケース・37製品要求を対応付け、1072ソース・687成果物を監査。**S10実験参照の認定、実PG全規模・全自由後流の製品性能、全実務受入・リリースは未完了。本体コードは未着手。**

[開始判定・実装順と残リスク](../plans/implementation-start-review.md)、[受入計画](../requirements/practical-acceptance-plan.md)、[ケース台帳](../../rd/acceptance/plan/records/practical-acceptance-case-catalog.json)。
<!-- active-uvlm-status:end -->

版：1.1、2026-10-29。対象はQ1〜Q27の要求を満たすために必要な技術決定。RF-04の製品範囲はQ27 Bで確定。最新の採用式・具体規約は構造／UVLM／J1／N4／IO4、受入計画はA5に基準化。最初のIM-00開始審査は合格。以下の初期棚卸し・追記はその経緯として保持する。

全体の現在地・6作業の終了条件は[ロードマップ](../../rd/planning/preimplementation-history.md)。進め方・成熟段階・開始条件は [実装前設計計画](../../rd/planning/preimplementation-plan.md)、I/Fと不変条件の作成要件は [設計契約・不変条件](../../rd/planning/design-contracts-and-invariants.md) を参照。

## 台帳の読み方

2026-10-31の現在位置：[AVL参照診断の終了判定](../../rd/aero/avl-reference/records/avl-reference-diagnostic-disposition.md) により、直近の診断を終了した。独立配分対照4件と監査は成立したが、M24の合算誤差予算は未評価で、この経路を数値精度の受入参照に採用しない。RF-03／E-AERO／VV-03とJ1は閉じない。以後の主作業は [構造モデルの採用審査](../../rd/planning/structural-model-adoption-review.md)。本表の初期状態と後続記録を統合し、採用が成立した対象だけを基準化する。

- P0〜P4は実装前設計計画の段階。期限は日付ではなく、どの段階を閉じるまでに決めるかを表す。
- 「閉じる証拠」は作成必須の成果物であり、現在存在する証拠ではない。出典・式・契約・ケースのIDと成果物の位置を追加してから基準化する。
- 状態は初期棚卸しの結果。候補調査や草案があっても「採用・整合確認済み」とはしない。
- RF系の用途・範囲・品質の変更はユーザー判断。それ以外はアシスタントが根拠付きの技術案を作る。要求変更を伴う場合はユーザーへ選択肢を示す。
- 先行決定は入力となる設計を示す。J1共同決定群はAE-04、AE-05、TM-02、TM-04、FL-01。候補を反復比較し、荷重・後流・積分・連成・線形化が同じ離散系として成立してから群全体を基準化する。無理に個別の先行順へ分解しない。

## A. 解析問題と共通規約

| ID | 決める内容 | 先行決定 | 閉じる証拠 | 期限・状態 |
|---|---|---|---|---|
| RF-01 | 代表的な静的・過渡・フラッタ問題、入力・観測量・比較条件 | 合意済み要求 | 寸法・材料・支持・流れ・入力履歴・QoIを持つ問題定義と要求対応。ゼロ近傍の基準量も特定 | A5でSI入力の実務統合問題、公開S10の具体run／条件、観測・閾値・不足時の停止を指定。19ケース／37製品要求を対応。製品受入は未完了 |
| RF-02 | 定量的な適用範囲と、計算可能性／物理的適用性の境界 | RF-01 | 幾何・厚さ・変形・ひずみ・流れ・支持・入力に関する仮定と検出／警告／拒否の区分 | P1・理論条件と診断の区分を整理。定量値は未確定 |
| RF-03 | 受入証拠の範囲、参照の採否、校正と評価の分離 | RF-01 | 参照の版・条件・利用条件・不確かさ・取得手順と、不足が影響する要求の一覧 | A5・S10診断で取得物・版・校正／評価分離・不確かさ不足を記録。S10実験参照は未認定、実務受入は保留。凍結した独立数値参照は移植に使用可能 |
| RF-04 | 初期応力・初期ひずみの初回対応範囲 | RF-01 | 外力による予荷重、指定初期ひずみ、既存応力場の入力を区別したユーザー決定と受入範囲 | Q27 Bの範囲を構造統合／J1に採用。指定初期ひずみ・平衡・幾何剛性・履歴を保持。A5のV-MAT／V-J1で製品移植を検証 |
| CV-01 | 座標・単位・符号・幾何向き・係数・基準点 | RF-01 | 現草案を採用式へ対応付けた記号辞書、変換式、符号・単位の独立確認例 | P1・草案 |
| CV-02 | 自由度・局所大域インデックス・配列配置・スケール | ST-01、ST-06、ST-08、AE-01 | 全成分の量・単位・並び・接続、拘束前後の対応、残差・未知量の尺度と配列辞書 | P3・未定義 |
| CV-03 | 材料軸・積層順・応力ひずみ測度・せん断・出力配置 | ST-03、ST-06、CV-01 | 基底と配置、層・評価点、テンソル変換、成分順、仕事共役性と出力定義 | P3・採用材料・梁・構造統合仕様で測度・基底・層／片側・Mandel／工学・断面共役量と物理wrenchを固定。構造1完了。製品出力型とファイル形式は4（[構造統合仕様](structure/structural-adopted-specification.md)、[終了判断](../../rd/structure/integration/records/structural-adoption-results.md)） |

## B. 構造の定式化と離散化

| ID | 決める内容 | 先行決定 | 閉じる証拠 | 期限・状態 |
|---|---|---|---|---|
| ST-01 | ANCFシェルの具体的定式化・節点DOF・形状関数・厚さ表現 | RF-02、CV-01 | 原式、補間と全必要微分、要素座標・単位、許容形状、三角形の受入／再格子化方針、候補比較 | P2/P3・SH-D01〜08の参照係数PG-ANCFを採用し構造1を完了。一般曲面・薄肉・物理PDE・後退テーパ・全位置トレースの成立根拠を保持。具体疎解法／精度は4、実務全形状は5（[構造統合仕様](structure/structural-adopted-specification.md)、[終了判断](../../rd/structure/integration/records/structural-adoption-results.md)） |
| ST-02 | 初期曲面・参照配置・幾何写像・歪みの扱い | ST-01 | 参照ヤコビアン、勾配の構成、零ひずみ条件、反転・退化・品質条件と例 | P2/P3・固定参照写像、実法線／厚さ、正則性、後退テーパ領域、接合物理面の対応を採用。構造1完了。外部格子の具体的変換と幾何誤差予算は4／5（[構造統合仕様](structure/structural-adopted-specification.md)、[終了判断](../../rd/structure/integration/records/structural-adoption-results.md)） |
| ST-03 | 等方・積層弾性の材料則と変形・ひずみの仮定 | ST-01、RF-02/04 | 初期ひずみを含むエネルギー密度／応力関係、材料軸変換、厚さ条件、入力と許容パラメータ、接線 | P2/P3・MAT採用仕様の全3D Green／第二PK・等価入力・積層・材料軸・初期ひずみ・両側回復を統合。構造1完了。実材料の用途別受入は5（[構造統合仕様](structure/structural-adopted-specification.md)、[終了判断](../../rd/structure/integration/records/structural-adoption-results.md)） |
| ST-04 | 求積・ロッキング対策・内部変数と静的縮約 | ST-01〜03 | 積分点・重み、ANS/EAS等の採否と具体式、内部反復・凝縮・接線整合、異常モードの確認方法 | P2/P3・採用PG空間・全3Dエネルギー・求積・参照条件付けの根拠を維持し構造1完了。未達の旧方式や人工補正を既定へ戻さない。具体精度／疎求解は4（[構造統合仕様](structure/structural-adopted-specification.md)、[終了判断](../../rd/structure/integration/records/structural-adoption-results.md)） |
| ST-05 | シェル・梁の質量、慣性項、整合／集中化 | ST-01、ST-06 | 質量式、構成依存性、速度依存慣性項、離散運動エネルギー、剛体質量特性 | P2/P3・シェル一定M、梁／剛体M(q)と慣性、重力、同一平衡M/K/Cを統合し構造1完了。動的接線M,qqと必要三階運動学は3／4へ明示して引継ぎ（[構造統合仕様](structure/structural-adopted-specification.md)、[終了判断](../../rd/structure/integration/records/structural-adoption-results.md)） |
| ST-06 | 梁理論、回転／勾配DOF、断面・材料、偏心 | RF-02、CV-01 | 軸・断面・せん断・ねじり・回転表現、剛性・質量と必要な接線、大回転の運動学、シェルとの境界 | P2/P3・BM-D01〜09のSE(3)梁を採用し構造1完了。軸・せん断・ねじり・曲げ、連成断面と初期ひずみ、自由反り・拘束反り拒否・局所回復・精度条件を維持（[構造統合仕様](structure/structural-adopted-specification.md)、[終了判断](../../rd/structure/integration/records/structural-adoption-results.md)） |
| ST-07 | 集中質量・重心偏心・慣性テンソルの組込み | ST-01、ST-06 | 基準点と軸、運動エネルギー・一般化慣性項、平行軸変換、要素DOFへの写像 | P2/P3・CN-D05/06の粒子／偏心剛体・M／必要微分・慣性／重力、質量の所有／出自／追加置換と二重計上拒否を統合し構造1完了（[構造統合仕様](structure/structural-adopted-specification.md)、[終了判断](../../rd/structure/integration/records/structural-adoption-results.md)） |
| ST-08 | 固定・剛体結合・交差接続・指定運動と拘束解法 | ST-01、ST-06、ST-07 | C(q,t)、階数、C_qと必要微分、速度・加速度整合、消去／乗数等の採否、反力、許容仮想変位 | P2/P3・CN-D01〜09の全位置接着・明示剛体取付・指定支持・冗長／不整合・反力と仕事を採用し構造1完了。動的拘束と初期化は3、疎構造は4、接合用途別収束は5（[構造統合仕様](structure/structural-adopted-specification.md)、[終了判断](../../rd/structure/integration/records/structural-adoption-results.md)） |
| ST-09 | 物理減衰、弾性・減衰支持、重力・外力・追従荷重 | ST-03、ST-08 | 荷重と仕事、支持ポテンシャル・散逸、構成依存荷重の接線、物理／数値減衰の区分 | P2/P3・SI-D04/05/07として荷重測度・物体／世界成分・圧力／全接線、客観粘性、弾性／減衰支持・基準仕事、明示時系列を採用し1.5完了。時刻／空力依存は2／3、実材料係数は5（[構造統合仕様](structure/structural-adopted-specification.md)、[終了判断](../../rd/structure/integration/records/structural-adoption-results.md)） |
| ST-10 | 内力・接線の導出・生成方法と応力／断面力回復 | ST-02〜09、CV-03 | エネルギーから内力・接線への対応、手導出／AD等の採否、内部変数を含む微分、応力・反力・断面積分の式 | P2/P3・SI-D06/08/09として全残差・全接線、測度／軸／層側の応力、物理切断面・原反力・収支、状態規約を統合し1.5完了。必要動的高階微分は3／4、製品生成・型は4、実務回復受入は5（[構造統合仕様](structure/structural-adopted-specification.md)、[終了判断](../../rd/structure/integration/records/structural-adoption-results.md)） |

## C. 空力の定式化と離散化

| ID | 決める内容 | 先行決定 | 閉じる証拠 | 期限・状態 |
|---|---|---|---|---|
| AE-01 | 渦の表現、幾何格子・渦格子・荷重点・循環未知量 | RF-02、CV-01 | 点・辺・面のIDと接続、具体配置、法線・有向辺・循環の符号、曲面と後縁の定義 | P3・空間方式／入力／主要契約を2.1〜2.4で統合採用。定常・分布・抗力・非定常12観測と誤差規則、全相互作用・全荷重／PG写像・状態／再開を監査。AE-04/05の製品時間・初期整合もJ1共同決定3で採用済み（[統合仕様](aero/uvlm-integrated-specification.md)、[終了判断](../../rd/aero/adoption/records/uvlm-adoption-completion.md)） |
| AE-02 | 非貫通・Kutta条件、影響係数・右辺と評価 | AE-01、AE-03 | 定式化と離散式、行列の次元、境界条件の組立、翼端・対称面の扱い、線形残差 | P3・空間方式／入力／主要契約を2.1〜2.4で統合採用。定常・分布・抗力・非定常12観測と誤差規則、全相互作用・全荷重／PG写像・状態／再開を監査。AE-04/05の製品時間・初期整合もJ1共同決定3で採用済み（[統合仕様](aero/uvlm-integrated-specification.md)、[終了判断](../../rd/aero/adoption/records/uvlm-adoption-completion.md)） |
| AE-03 | 有限渦線の誘導速度、コア・近接・自己誘導 | AE-01、CV-01 | 採用式、特異・退化条件、コアモデルと尺度、除外／停止規則、独立参照 | P3・空間方式／入力／主要契約を2.1〜2.4で統合採用。定常・分布・抗力・非定常12観測と誤差規則、全相互作用・全荷重／PG写像・状態／再開を監査。AE-04/05の製品時間・初期整合もJ1共同決定3で採用済み（[統合仕様](aero/uvlm-integrated-specification.md)、[終了判断](../../rd/aero/adoption/records/uvlm-adoption-completion.md)） |
| AE-04 | 循環・非循環荷重、圧力差・力・作用点・モーメント | AE-01/02、J1 | 採用する荷重式、各項の速度・循環・履歴・時刻、自己誘導の扱い、係数と作用点、二重計上防止 | P3・空間方式／入力／主要契約を2.1〜2.4で統合採用。定常・分布・抗力・非定常12観測と誤差規則、全相互作用・全荷重／PG写像・状態／再開を監査。AE-04/05の製品時間・初期整合もJ1共同決定3で採用済み（[統合仕様](aero/uvlm-integrated-specification.md)、[終了判断](../../rd/aero/adoption/records/uvlm-adoption-completion.md)） |
| AE-05 | 後流初期化、放出、移流、反復内の更新、保持・遠方処理 | AE-01〜03、J1 | 初期平衡後流、接続と循環受渡し、離散更新、放出時刻、打切り・近似の条件と誤差評価 | P3・空間方式／入力／主要契約を2.1〜2.4で統合採用。定常・分布・抗力・非定常12観測と誤差規則、全相互作用・全荷重／PG写像・状態／再開を監査。AE-04/05の製品時間・初期整合もJ1共同決定3で採用済み（[統合仕様](aero/uvlm-integrated-specification.md)、[終了判断](../../rd/aero/adoption/records/uvlm-adoption-completion.md)） |
| AE-06 | 複数揚力面・後流の相互干渉、近接・交差、対称性 | AE-01〜05 | 干渉の全組合せ、共有辺・重複・自己項、幾何異常とモデル適用性の区分 | P3・空間方式／入力／主要契約を2.1〜2.4で統合採用。定常・分布・抗力・非定常12観測と誤差規則、全相互作用・全荷重／PG写像・状態／再開を監査。AE-04/05の製品時間・初期整合もJ1共同決定3で採用済み（[統合仕様](aero/uvlm-integrated-specification.md)、[終了判断](../../rd/aero/adoption/records/uvlm-adoption-completion.md)） |
| AE-07 | 一様流・突風・指定運動・始動条件・入力補間 | RF-01、CV-01、AE-02 | 空間・時間・座標の定義、既知速度と翼面速度の区分、始動履歴、欠測・範囲外・イベントの扱い | P3・空間方式／入力／主要契約を2.1〜2.4で統合採用。定常・分布・抗力・非定常12観測と誤差規則、全相互作用・全荷重／PG写像・状態／再開を監査。AE-04/05の製品時間・初期整合もJ1共同決定3で採用済み（[統合仕様](aero/uvlm-integrated-specification.md)、[終了判断](../../rd/aero/adoption/records/uvlm-adoption-completion.md)） |

## D. 写像・連成・時間・安定性

| ID | 決める内容 | 先行決定 | 閉じる証拠 | 期限・状態 |
|---|---|---|---|---|
| TR-01 | 空力面と構造の対応、埋込み、オフセット、探索 | ST-01/02/06/08、AE-01 | x=g(q,t)、物質座標、位置・速度・法線、Jと明示的時間依存、写像不能の扱い | P1/P2・固定物質／厚さ座標の曲面Q9で頂点・辺中点・三角重心の写像と法線を照合。多チャート・一般オフセット・探索の採用は未完了 |
| TR-02 | 一般化荷重・仕事・合力モーメントと写像の微分 | TR-01、AE-04、ST-08 | Q=J^T f等の適用条件、独立モーメント・偏心、荷重の構成依存性と接線、保存する量の成立条件 | P3・J1-D01〜09の数式・状態・必要微分・拒否／精度規約を共同採用。実PG／全自由後流の変形平衡・全接線・全履歴線形化・時系列、初期微分・姿勢・再開・枝／探索を限定検証。具体数値解法／抽出器・型／I/FはN4／IO4で採用済み。受入計画はA5、製品移植・受入結果はIM-00〜08（[仕様](coupled/coupled-adopted-specification.md)、[終了判断](../../rd/coupled/j1/records/j1-adoption-completion.md)） |
| TM-01 | 静的／動的な初期化と整合した初期履歴 | ST-03/08/09、AE-05/07 | 初期ひずみを保持した平衡／非平衡の区分、初期q,v,a,lambdaと循環・後流、拘束整合、始動と履歴の生成順 | P3・J1-D01〜09の数式・状態・必要微分・拒否／精度規約を共同採用。実PG／全自由後流の変形平衡・全接線・全履歴線形化・時系列、初期微分・姿勢・再開・枝／探索を限定検証。具体数値解法／抽出器・型／I/FはN4／IO4で採用済み。受入計画はA5、製品移植・受入結果はIM-00〜08（[仕様](coupled/coupled-adopted-specification.md)、[終了判断](../../rd/coupled/j1/records/j1-adoption-completion.md)） |
| TM-02 | 構造・後流の時間積分、刻み、評価時刻と履歴 | ST-05/08/10、J1 | 採用積分式、パラメータ、ステージ時刻、荷重・循環差分の順、固定／可変刻み、数値減衰、誤差次数 | P3・J1-D01〜09の数式・状態・必要微分・拒否／精度規約を共同採用。実PG／全自由後流の変形平衡・全接線・全履歴線形化・時系列、初期微分・姿勢・再開・枝／探索を限定検証。具体数値解法／抽出器・型／I/FはN4／IO4で採用済み。受入計画はA5、製品移植・受入結果はIM-00〜08（[仕様](coupled/coupled-adopted-specification.md)、[終了判断](../../rd/coupled/j1/records/j1-adoption-completion.md)） |
| TM-03 | 非線形静的平衡、荷重増分・分岐・継続と停止 | ST-08〜10、AE-02/04、TR-02 | 離散残差・接線・初期推定、緩和／線探索等の採否、平衡分岐、特異点・未収束・反力の判定 | P3・J1-D01〜09の数式・状態・必要微分・拒否／精度規約を共同採用。実PG／全自由後流の変形平衡・全接線・全履歴線形化・時系列、初期微分・姿勢・再開・枝／探索を限定検証。具体数値解法／抽出器・型／I/FはN4／IO4で採用済み。受入計画はA5、製品移植・受入結果はIM-00〜08（[仕様](coupled/coupled-adopted-specification.md)、[終了判断](../../rd/coupled/j1/records/j1-adoption-completion.md)） |
| TM-04 | 強連成の離散問題、反復、緩和・接線・収束 | TR-02、J1 | 全未知量と残差、予測→反復→確定の擬似コード、後流・履歴の扱い、正規化と許容差、近似接線の省略項 | P3・J1-D01〜09の数式・状態・必要微分・拒否／精度規約を共同採用。実PG／全自由後流の変形平衡・全接線・全履歴線形化・時系列、初期微分・姿勢・再開・枝／探索を限定検証。具体数値解法／抽出器・型／I/FはN4／IO4で採用済み。受入計画はA5、製品移植・受入結果はIM-00〜08（[仕様](coupled/coupled-adopted-specification.md)、[終了判断](../../rd/coupled/j1/records/j1-adoption-completion.md)） |
| TM-05 | 試行・棄却・刻み変更・停止・確定の状態遷移 | TM-01〜04 | 変更可能な状態、成功／失敗ごとの遷移、履歴・後流の更新回数、キャッシュ世代、再試行時の復元 | P3・J1-D01〜09の数式・状態・必要微分・拒否／精度規約を共同採用。実PG／全自由後流の変形平衡・全接線・全履歴線形化・時系列、初期微分・姿勢・再開・枝／探索を限定検証。具体数値解法／抽出器・型／I/FはN4／IO4で採用済み。受入計画はA5、製品移植・受入結果はIM-00〜08（[仕様](coupled/coupled-adopted-specification.md)、[終了判断](../../rd/coupled/j1/records/j1-adoption-completion.md)） |
| FL-01 | 平衡周りの線形化と安定性問題 | TM-03、ST-10、TR-02、J1 | 拘束・空力・後流を含む連続／離散系、摂動変数と演算子、数値減衰、過渡系との対応と比較方法 | P3・J1-D01〜09の数式・状態・必要微分・拒否／精度規約を共同採用。実PG／全自由後流の変形平衡・全接線・全履歴線形化・時系列、初期微分・姿勢・再開・枝／探索を限定検証。具体数値解法／抽出器・型／I/FはN4／IO4で採用済み。受入計画はA5、製品移植・受入結果はIM-00〜08（[仕様](coupled/coupled-adopted-specification.md)、[終了判断](../../rd/coupled/j1/records/j1-adoption-completion.md)） |
| FL-02 | 固有値解法、構造・空力の縮約とモード対応 | FL-01、NU-01 | 残差、対象固有対、縮約・モード数の収束、近接／重複モードの追跡、質量内積等の定義 | P3・J1-D01〜09の数式・状態・必要微分・拒否／精度規約を共同採用。実PG／全自由後流の変形平衡・全接線・全履歴線形化・時系列、初期微分・姿勢・再開・枝／探索を限定検証。具体数値解法／抽出器・型／I/FはN4／IO4で採用済み。受入計画はA5、製品移植・受入結果はIM-00〜08（[仕様](coupled/coupled-adopted-specification.md)、[終了判断](../../rd/coupled/j1/records/j1-adoption-completion.md)） |
| FL-03 | 速度探索、臨界値精度、発散・フラッタ・不明の分類 | FL-01/02、RF-01 | 探索区間・追加点・停止、平衡分岐とモードの追跡、成長率の不確かさ、再開データ | P3・J1-D01〜09の数式・状態・必要微分・拒否／精度規約を共同採用。実PG／全自由後流の変形平衡・全接線・全履歴線形化・時系列、初期微分・姿勢・再開・枝／探索を限定検証。具体数値解法／抽出器・型／I/FはN4／IO4で採用済み。受入計画はA5、製品移植・受入結果はIM-00〜08（[仕様](coupled/coupled-adopted-specification.md)、[終了判断](../../rd/coupled/j1/records/j1-adoption-completion.md)） |

## E. 数値基盤・契約・外部形式・検証

| ID | 決める内容 | 先行決定 | 閉じる証拠 | 期限・状態 |
|---|---|---|---|---|
| NU-01 | 線形・非線形・固有値解法と外部数値実装の利用 | ST-08/10、AE-02、TM-04、FL-01 | 行列の対称性・正定値性等の条件、直接／反復・前処理、悪条件・特異系診断、利用／自作の境界と版・利用条件 | N4-01〜10／IO4-01〜07で型・所有権・寿命・疎密解法・精度・失敗・再開・容量政策を採用。両OSの限定試作を監査。製品移植・全実規模性能・実務受入は未完了（[開始審査](../plans/implementation-start-review.md)） |
| NU-02 | 浮動小数点、スケーリング、並列化、再現性 | CV-02、TM-02/04、NU-01 | 精度型、残差尺度、演算順、最適化フラグ、決定性の範囲、許容差・診断 | N4-01〜10／IO4-01〜07で型・所有権・寿命・疎密解法・精度・失敗・再開・容量政策を採用。両OSの限定試作を監査。製品移植・全実規模性能・実務受入は未完了（[開始審査](../plans/implementation-start-review.md)） |
| NU-03 | 行列・演算子・キャッシュ・状態の記憶と容量 | CV-02、AE-05、FL-01/02、IF-01 | 実務規模のメモリと計算量、疎パターン・分解・後流、行列明示／作用形式、更新条件と支配項の見積もり | N4-01〜10／IO4-01〜07で型・所有権・寿命・疎密解法・精度・失敗・再開・容量政策を採用。両OSの限定試作を監査。製品移植・全実規模性能・実務受入は未完了（[開始審査](../plans/implementation-start-review.md)） |
| IF-01 | Model/State/Trial/Workspace/Resultとデータ辞書 | CV-02/03、TM-05、FL-01 | フィールド、型・次元・配置・所有者・寿命、依存・可変性、シリアライズ対象、キャッシュ世代 | N4-01〜10／IO4-01〜07で型・所有権・寿命・疎密解法・精度・失敗・再開・容量政策を採用。両OSの限定試作を監査。製品移植・全実規模性能・実務受入は未完了（[開始審査](../plans/implementation-start-review.md)） |
| IF-02 | 主要処理のI/Fと擬似コード | IF-01、ST/AE/TR/TM/FLの対象決定 | 契約・不変条件文書の対象処理に事前事後条件、時間、引数・返値、変更範囲、失敗保証、状態遷移と計算順を記載 | N4-01〜10／IO4-01〜07で型・所有権・寿命・疎密解法・精度・失敗・再開・容量政策を採用。両OSの限定試作を監査。製品移植・全実規模性能・実務受入は未完了（[開始審査](../plans/implementation-start-review.md)） |
| IF-03 | 入力検査・診断・失敗分類・許容差の解決 | RF-02、CV-02、IF-02 | 未知・重複・非有限値・不整合・数値失敗の検出箇所、原因ID、終了状態、失敗時の保証 | N4-01〜10／IO4-01〜07で型・所有権・寿命・疎密解法・精度・失敗・再開・容量政策を採用。両OSの限定試作を監査。製品移植・全実規模性能・実務受入は未完了（[開始審査](../plans/implementation-start-review.md)） |
| IF-04 | 結果・チェックポイントの形式と互換性 | CV-03、TM-05、FL-03、IF-01 | 保存項目の具体型、版・ハッシュ、確定単位、出力継続、書込み失敗、CLI/Python共有契約 | N4-01〜10／IO4-01〜07で型・所有権・寿命・疎密解法・精度・失敗・再開・容量政策を採用。両OSの限定試作を監査。製品移植・全実規模性能・実務受入は未完了（[開始審査](../plans/implementation-start-review.md)） |
| IF-05 | Nastran/Gmshと内部モデルの変換 | ST-01/06〜08、CV-01/02、IF-03 | 対応版・カード・要素、座標・単位・材料・接続変換、非対応の拒否、変換前後の比較ケース | N4-01〜10／IO4-01〜07で型・所有権・寿命・疎密解法・精度・失敗・再開・容量政策を採用。両OSの限定試作を監査。製品移植・全実規模性能・実務受入は未完了（[開始審査](../plans/implementation-start-review.md)） |
| PK-01 | C++・Python境界、実行時依存、配布・対応版 | IF-01〜04、NU-01/02 | ABIと配列所有権、例外・中断、対象OS/CPU/Python、wheel/CLIの構成、依存・通知・導入試験計画 | 対象境界の実装前・方針のみ |
| VV-01 | 要求・式・処理・ケースの対応と独立参照 | RF-03、対象のST/AE/TR/TM/FL | ケースごとの数値入力、参照または生成法、出典、適用条件、評価区間、許容差・合否・追跡 | P3・未完成 |
| VV-02 | 不変条件の成立条件、検査方法、尺度・閾値 | 対象の理論決定、IF-02 | 不変条件カタログを式と検証点へ落とし、比較誤差と適用除外を事前指定 | P3・作成要件を定義した段階 |
| VV-03 | 空間・時間・後流・反復・縮約の誤差予算 | TM-02/04、FL-02、NU-02、VV-01 | 独立に変えるパラメータ、細分化列、基準量、推定法、数値誤差予算、非漸近時の扱い | P3・20%方針のみ |
| VV-04 | 異常系・再開・配布・性能の受入ケース | IF-03〜05、PK-01、NU-03、VV-01 | 破損・未収束・再開・clean installの具体的入力と操作、期待状態・診断、性能ケースの全設定 | 対象実装前に仕様化・ケース名のみ |

## 最初に閉じる範囲

RF-01〜03とCV-01を起点に、全体系の問題定義を作る。構造・空力候補を調べた後、J1を含む離散系を閉じる。主要I/Fはその結果から導出する。IF-02だけ先に関数名・クラス名を確定しない。

実装開始判定では、対象パッケージについてこの台帳から必要な先行・共同決定をたどり、未基準化の項目が残らないことを確認する。例えば誘導速度カーネルにも、CV-01、AE-01/03、必要なIF・VV契約と全体系での位置付けが必要。ライブラリがビルドできたことだけでは解除しない。

この台帳は決めるべきことの初期棚卸しであり、採用仕様そのものではない。調査で新しい決定が必要と分かった場合は、ID、依存、影響、期限を追加し、既存の基準化済み項目への影響を点検する。

2026-10-29追記：[全体系の問題定義](coupled-problem-definition.md) と [要求・証拠対応](../../rd/planning/requirements-evidence-map.md) を作成。全体系の式、仕事・状態の接続、平衡／過渡／安定性の関係まで整理したが、具体離散化・実務ケースの数値条件・参照採用は未完了。どの項目も本追記だけで基準化済みにはしない。

Q27 Bに基づく [初期ひずみ設計](../../rd/structure/material/initial-strain-design.md) と [ANCFシェル候補の一次比較](../../rd/structure/main-shell/q9-plate/ancf-shell-candidate-audit.md) を作成。前者は要求・候補式・解析参照、後者は一次資料の照合状況。いずれも要素・材料の採用完了ではない。

[3423数式仕様案](../../rd/structure/main-shell/ancf-3423/ancf-3423-formulation.md) を追加。ANCFの勾配を含む24座標、初期曲面の参照写像、ANS/EAS、初期ひずみ、質量、内力・縮約接線、回復・呼出契約を具体化した。固定版ソースの懸案と解析参照は対応文書へ分け、数値試験や実装開始の合格とは扱わない。

[層別EAS整合審査](../../rd/structure/main-shell/ancf-3423/eas-layer-consistency.md) でEAS-D01〜04、ANS-D01を記録。層ごとの零平均条件・基底・局所解法は技術判断を進めたが、3423主シェル採用やST-04全体の基準化とは判定しない。

[曲面対応シェル比較](../../rd/structure/main-shell/q9-plate/curved-shell-comparison.md) を追加。3833型の厚さ二次補間を固定版CPPと照合し、48座標Q8型の論文とは別候補として扱う。数学仕様11文書の主要式・添字・行列もLaTeXへ整備し、数式のコード表記を解消した。表記変更を理論・数値検証の完了証拠には数えない。

3833-D01〜08を [数式仕様](../../rd/structure/main-shell/q9-plate/ancf-3833-formulation.md) に記載。72座標・必要微分・参照写像・材料基底・初期ひずみ・内力・接線・質量・求積を具体化。[事前仕様付き代数試験](../../rd/structure/main-shell/q9-plate/records/3833-algebraic-experiment.md) の11判定が成立。3×3×2弾性経路は平坦な無応力1要素でrank66、同じ求積の質量はrank48。厚さ2点質量と安定化なしの2×2×2弾性を除外し、3×3×2弾性を比較候補に残す。ST-01〜05/10全体の完了ではない。

[薄肉・歪み試験](../../rd/structure/main-shell/q9-plate/records/3833-bending-distortion-experiment.md) の76判定は、純曲げ・人工せん断の独立解析参照との一致を確認。8節点台形のひずみ汚染を踏まえて主シェルへの単純採用を見送り。[9節点設計](../../rd/structure/main-shell/q9-plate/q9-thickness-quadratic-design.md) で不足多項式の保持と27節点二次体積要素への固定座標変換を導出した。荷重平衡・格子収束・非線形・接続・独立FEM照合は未完了。

続く [外部FEM・荷重試験](../../rd/structure/main-shell/q9-plate/records/q9-independent-load-experiment.md) でST-01/04/05の証拠を追加。81座標の線形質量・接線・面荷重と単要素純曲げ平衡を独立実装で照合した。NU-02には薄肉時の組立後変換・求解の精度問題と、積分前変換・反復改良の改善結果を追加。製品の尺度・誤差評価・解法契約は未確定で、依存ライブラリ採用や各項目全体の基準化とはしない。次は格子収束の事前仕様。

[非一様板曲げ・格子収束](../../rd/structure/main-shell/q9-plate/records/q9-plate-convergence-experiment.md) で面内連続・厚さ二次の参照を導出し、ST-04/VV-03の限定証拠を追加。厚さ比0.1・0.01の最細格子は1%誤差基準を満たす。最薄の細格子と粗格子求積は未達。NU-02の残差評価・精度、ST-04の面荷重積分を切り分けるまで最薄の収束・ロッキングの結論を保留する。全体の基準化・実装開始にはしない。

[精度切り分け](../../rd/structure/main-shell/q9-plate/records/q9-precision-diagnostic-experiment.md) の別経路で平坦・均一等方・矩形格子の全12ケースが成立し、最薄16×16の収束判定を行える証拠を追加。NU-02には [精度契約案](foundation/linear-solve-accuracy-contract.md) を作成した。高精度保持解と倍精度投影の残差は区別する。ST-01/04/05・NU系全体の基準化、製品用Decimal採用、歪み・曲面への一般化とはしない。

[歪み・曲面荷重試験](../../rd/structure/main-shell/q9-curved/records/q9-curved-load-experiment.md) の20ケース・座標変換10組が成立。ST-02には格納幾何に対する有理数Bernstein係数の全域判定、荷重・ゲージ・剛体運動の限定証拠を追加し、[参照幾何契約案](../../rd/structure/conventions/reference-geometry-contract.md) を作成した。h・r/d/cと物理厚さ／中央面の区別、界面の厚さ座標再パラメータ化を明示。ST-08の界面変換の数値照合、一般曲面の収束、材料・非線形等の採用は未完了。

[厚さ座標界面契約](../../rd/structure/conventions/thickness-coordinate-interface.md) と [二要素試験](../../rd/structure/main-shell/q9-plate/records/q9-interface-reparameterization-experiment.md) で正常10ケース・負の対照6件が成立。IFACE-D01〜03として係数基底変換、力の双対／行列合同変換、右手系・物理層IDの保持を採用。実形状の厚さ変更・異なる法線・交差部材の接続、積層異方性・初期ひずみ・非線形の検証とはしない。次はST-03の材料軸・層・初期ひずみと内力・接線の規約を審査する。

[材料・初期ひずみ契約](structure/material-initial-strain-contract.md) にMAT-D01〜03を記録し、[材料点・Q9接線試験](../../rd/structure/main-shell/q9-dynamics/records/q9-material-tangent-experiment.md) で限定照合。参照軸・測度・工学せん断・幾何剛性、非ゼロ内力の回転恒等式、層場と厚さ座標変換を確認した。合成弾性定数の数式確認であり、実材料適用限界・任意有限初期ひずみ・非線形平衡・固有値・過渡・再開・製品採用は未完了。

[平衡・線形化契約](../../rd/structure/main-shell/q9-dynamics/prestress-equilibrium-linearization.md) のINIT-D01〜03と [4座標拘束のQ9試験](../../rd/structure/main-shell/q9-dynamics/records/q9-prestress-equilibrium-experiment.md) を追加。総ポテンシャル線探索の丸めによる未達を保持し、同じエネルギーの直接増分経路で平衡・固有値・微小応答・離散固有値・振幅極限を照合。非線形拘束の接線曲率と慣性を含めた。同じモデルの受渡しの限定証拠で、TM-02/J1、一般の平衡・連成安定性・再開・製品解法の採用ではない。

[全座標平衡契約](../../rd/structure/main-shell/q9-plate/q9-full-space-equilibrium-contract.md) のFULL-D01〜03と [単要素試験](../../rd/structure/main-shell/q9-plate/records/q9-full-space-equilibrium-experiment.md) を追加。全81座標へ12線形支持を課し、非アフィンな初期推定から69未知量で平衡・反力・予応力スペクトルを照合。外部Hex27基底の独立組立、69モードの部分空間、線形応答の刻み収束が成立した。モーメント基準点の実装不整合と修正前結果を保持。多要素・一般曲面・積層・非線形過渡・連成・再開の採用は未完了。

[多要素収束契約](../../rd/structure/main-shell/q9-dynamics/multielement-prestress-convergence-contract.md) にMULTI-D01〜03を記録。[非一様平衡の内部4格子・外部2格子試験](../../rd/structure/main-shell/q9-dynamics/records/q9-multielement-prestress-convergence.md) で、静的誤差・有限参照の変化・物理モード対応を分離して照合した。成分別参照組立は全225座標で元の式と一致。数値的な比較指標であり、真の固有値誤差の厳密上界・純薄板曲げ・曲面・積層・実翼・製品の組立方式の採用完了とはしない。

[非線形ステップ・再開契約](../../rd/structure/main-shell/q9-dynamics/nonlinear-step-restart-contract.md) にSTEP-D01〜03を追加。確定状態の不変性、親世代、一度だけの確定、初期ひずみ・時間荷重を含むモデル識別を具体化した。[全69自由度の限定試験](../../rd/structure/main-shell/q9-dynamics/records/q9-nonlinear-restart-experiment.md) では時間製造解の2次収束と棄却・再開一致が成立。一方、自然応答は2560刻みでも事前次数基準未達。TM-02/J1採用、自然応答の時間精度、後流・連成・製品再開の完了とはしない。

[時間誤差・刻み契約](../../rd/structure/main-shell/q9-dynamics/transient-error-and-step-contract.md) にTIME-D01〜03を追加。[自然応答診断](../../rd/structure/main-shell/q9-dynamics/records/q9-transient-accuracy-diagnostic.md) で、全69モードを保持した増分内力・接線の同値性、反復許容差感度、短窓の別積分法照合と2次収束を確認。初回1280刻みの速度誤差未達を残し、同じ閾値の2560刻み経路が成立した。半周期の高次速度成分の差は未解消であり、ステップ受理と時間精度適格性を区別する。一般の残差尺度・帯域・積分法J1の採用完了とはしない。

[保持する応答・初期化・積分候補](../../rd/structure/main-shell/q9-dynamics/transient-bandwidth-and-integrator-design.md) にBAND-D01・INIT-D04・INT-D01を追加。未処理出力を基本とし、4初期化型、一般化α法の係数規約と算法加速度履歴を定義。[限定比較](../../rd/structure/main-shell/q9-dynamics/records/q9-integrator-candidate-comparison.md) では線形増幅・2次精度・中点同値性・履歴保持再開が成立したが、同じ初期値問題の速度・物理加速度・応力増分は候補ごとに未達。数値散逸を既定採用せず、保持荷重平衡からの解放を別初期化として比較する。J1と製品方式の選定は未完了。

[保持荷重解放契約](../../rd/coupled/j1/load-release-initialization-contract.md) にRELEASE-D01〜03を追加。[全69自由度の限定試験](../../rd/structure/main-shell/q9-dynamics/records/q9-load-release-experiment.md) で保持体積力の独立積分、非線形保持平衡、解放後初期加速度、イベント初期状態の保存を確認した。T1/20の未処理観測量は中点則で成立。一般化α法は刻み二乗付き残差の固定許容差で初回未達となり、動的な力・速度尺度の追加経路で同じ誤差基準を満たした。初回未達は保存。一般拘束・任意時刻イベント・連成方式の採用は未完了。

[過渡イベント確定契約](../../rd/structure/main-shell/q9-dynamics/transient-event-transaction-contract.md) にEVENT-D01〜03を追加。[移動中の途中解放試験](../../rd/structure/main-shell/q9-dynamics/records/q9-event-transaction-experiment.md) で境界整列、179受理ステップと1イベントの分離、181出力キー、4保存位置からの全確定状態一致、13不正チェックポイント拒否を確認した。初回の検査集計不整合と修正を記録し、軌道ハッシュは不変。時刻型・一般イベント・障害耐久性・連成系の採用ではない。次はJ1の循環・荷重・後流時刻と状態集合を具体化する。

[UVLM時刻・循環・荷重・後流契約](../../rd/aero/wake-loads/uvlm-discrete-time-contract.md) にUVT-D01〜03を追加。B1一次比較基準の時刻表・旧循環の後流受渡し・構造との一括残差、B2微分の履歴要件を定義した。[2面の限定試験](../../rd/aero/wake-loads/records/uvlm-time-topology-experiment.md) で20判定が成立。誘導と力の辺集合、付着境界と物質列、反復と受理を区別。製造右辺・指定移流・処方荷重であり、物理揚抗力・自由後流・ANCF強連成・製品方式の採用とはしない。

[渦コア・自己項・荷重契約](../../rd/aero/kernel/vortex-core-self-and-load-contract.md) にCORE-D01〜03を追加。[限定試験](../../rd/aero/kernel/records/vortex-core-and-steady-load-experiment.md) で30線分・固定コア渦輪・定常平板の解析極限と幅感度を照合し、16判定が成立。自己直線と全渦輪、有限線分と延長線距離、幾何接触と有限な正則化速度、共有辺半径を区別。初回のJSON保存不整合を修正・記録した。製品のコア半径・成長、統一渦格子・曲面非定常荷重・J1は未選定。

[渦格子・物質対応仕様](../../rd/aero/kernel/aerogrid-material-load-geometry.md) にAGRID-D01〜03を追加。[曲面写像・キャンバ・瞬時付加質量試験](../../rd/aero/kernel/records/aerogrid-camber-addedmass-experiment.md) で同じ幾何候補の荷重点・面積・双対写像を照合した。初回で未検査だった曲面法線を修正し、初回記録を保持した上で法線差分・未提供拒否を追加。最小キャンバの揚力差0.01498%、瞬時非循環モーメント差0.006103%。全非定常応答・強連成・多チャート・製品方式の採用ではない。


[線形後流・調和応答契約](../../rd/aero/wake-loads/uvlm-linear-wake-response-contract.md) にLWAKE-D01〜03を追加。[調和上下運動試験](../../rd/aero/unsteady-reference/records/uvlm-harmonic-heave-experiment.md) で旧循環放出・有限記憶後流・BDF1荷重を共通幾何へ接続した。周期履歴の時間前進は同じ離散式の周波数解と一致したが、独立Theodorsen参照の初回14判定中5件、追加7判定中3件が未達。追加最細のk=1は複素差2.0312%、振幅差2.0279%で2%基準を満たさない。空間と時間の誤差相殺が疑われるため、rの事後選択・J1採用は行わない。次は後縁・新生後流・荷重求積の整合と刻み極限を切り分ける。


[後縁・後流区間・荷重の整合診断](../../rd/aero/wake-loads/uvlm-trailing-edge-consistency.md) にTE-D01〜03を追加。μの連続・移流と渦層勾配／後縁圧力差を区別し、W-P1節点履歴・先頭陰的項を導出した。[限定試験](../../rd/aero/wake-loads/records/uvlm-trailing-edge-diagnostic.md) は解析縮約と3D渦辺、左右重み、周期時間前進を照合。一様M=128は連続後流・厳密微分でもk=1の複素差2.85%が残り、W-P1だけでは精度未達。物理パネル分割をcosine配置へ変えた周波数診断はM=256で最大複素差0.272%、128→256差0.269%となり追加条件が成立した。初回未達とJSON保存不整合を保存。非一様格子の実際の時間前進、後縁圧力・モーメント・自由後流・強連成とJ1は未完了。


[非一様W-P1・弱圧力契約](../../rd/aero/wake-loads/uvlm-graded-wake-pressure-contract.md) にGWP-D01〜03を追加。[限定試験](../../rd/aero/wake-loads/records/uvlm-graded-wake-pressure-experiment.md) で非一様束縛格子・節点後流・先頭陰的項・BDF1荷重を同じ時間前進へ接続し、17判定が成立した。M=256の総揚力複素差は最大0.504%、k=.5の時間前進の独立理論差は最細0.578%。圧力平均の保存則・跳び・時間／周波数一致を照合したが、後縁近傍平均/p*は最細k=1でも約0.664であり、点wise Kutta圧力の合格とはしない。製品方式・局所圧力・ピッチ／モーメント・自由後流・強連成・J1は未採用。


[局所圧力・後縁極限・ピッチ参照](../../rd/aero/wake-loads/uvlm-local-pressure-pitch-reference.md) にLPP-D01〜03を追加。[限定試験](../../rd/aero/wake-loads/records/uvlm-local-pressure-pitch-experiment.md) は独立な加速度ポテンシャル圧力、上下／2回転軸ピッチのF/M、保存則・基準点移動・周期時間前進を照合。初回18判定中16成立・局所圧力2未達。別の1024/2048細分化で最大局所参照差は0.441%へ減少したが、隣接差2.319%、非貫通残差4.36e-11で追加2条件未達。初回・追加結果を保持し、次は求解残差の精度評価とP0トレースの空間誤差を切り分ける。局所圧力・点wise Kutta条件・一般3D・J1の採用とはしない。


[残差・トレース精度契約](../../rd/aero/wake-loads/uvlm-residual-and-trace-accuracy.md) にRTA-D01〜03を追加。[限定診断](../../rd/aero/wake-loads/records/uvlm-residual-trace-experiment.md) でTwoProductを含む残差、80桁Decimal選択行、同一A/Bの反復改良を照合。2048分割の初期補償残差も最大4.40e-11で未達だったが、float64の1補正で最大6.55e-12へ改善。圧力変化は最大2.95e-11で、空間変化2.319%は未達のまま。独立μ参照で区間積分と両端トレースを分解し、最悪粗格子区間では合成トレース項が大きいことを確認。求解器・数値型・圧力再構成・一般3D・J1の採用とはしない。


[実数基底・複素後縁列更新契約](../../rd/aero/wake-loads/uvlm-real-base-column-update.md) にRCU-D01〜03を追加。[固定圧力区間の追加収束](../../rd/aero/wake-loads/records/uvlm-column-update-convergence.md) は同じ有効A/Bを実数LU＋複素1列更新で解き、4096分割の最大局所平均誤差0.350%、2048→4096差0.365%で13条件が成立。4096分割k=.1はfloat64補正が停滞し、高低保持で成立したが投影後の2右辺は残差未達。型・消費表現を省略しない。圧力場・観測端点・参照は不変。次は有限幅3DのW-P1・誘導・荷重写像を具体化する。汎用ソルバー・製品数値型・点圧力・J1は未採用。

### 有限幅3次元W-P1の作用素・状態候補

[3D契約](../../rd/aero/wake-loads/uvlm-wp1-three-dimensional-contract.md) にW3D-D01〜03を追加。AE-01/02/05、TM-04の部分具体化として、幾何頂点とは別の列×帯μ、一次境界渦線＋面内渦層、同一コア、先頭後流の複数列陰的項を定義した。[照合](../../rd/aero/wake-loads/records/uvlm-wp1-three-dimensional-experiment.md) は独立面積分・P0細分化・共有辺・製造非貫通を確認。基礎17判定のうち幅一列還元だけ未達を保存し、事前追加した求積比較3判定は成立した。一般3Dの面外誘導候補までで、近接求積・自己誘導・自由後流・物理翼荷重・全接線・J1は未完了。

### 近接求積と有限翼の物理荷重比較

[近接契約](../../rd/aero/kernel/uvlm-near-field-quadrature-contract.md) にNQ-D01〜03を追加。距離・基底別誤差・最終監査・丸め床・失敗伝播を具体化した。[事前固定した診断](../../rd/aero/kernel/records/uvlm-near-field-screening-experiment.md) では96次でも距離d/c=.01のL誤差20.9%が残る。診断6条件は成立したが、適応評価器は未実装・未採用。[有限翼比較仕様](../../rd/aero/finite-wing/records/finite-wing-load-reference-specification.md) のFW-D01〜03は後退テーパ／上反角の2翼を対象とし、独立AVL 3.52の軸・コア・荷重オプションを原本監査した。20入力を準備したが参照計算は未実施。抗力二次量・非定常・実験・J1の完了とはしない。

### 適応求積候補とAVL出力の疎通

[適応求積試験](../../rd/aero/unsteady-reference/records/uvlm-adaptive-quadrature-experiment.md) は12判定成立。平面の近接14列、一般3Dの面積分参照12列、一定μの接近18評価、変換・反復・失敗処理を限定照合した。単一セルの固定幾何候補であり、一般3D近接の可変μ各列・バッチ・全接線・荷重予算は未完了。[AVL実行プロトコル](../../rd/aero/avl-reference/records/avl-reference-execution-protocol.md) では初回の粘性フラグ不適合を保持し、修正3試行の設定・出力・反復一致を確認。FEはdCp等で全循環・全要素力ではなく、非貫通残差も未取得。本比較へ入る前に参照ソースビルド・コア監査・計測出力を整備する。J1は未選定。

### AVL原本ビルド・計測とコア0の採用保留

[ソース計測試験](../../rd/aero/avl-reference/records/avl-source-measurement-experiment.md) でGNU Fortran/C 14.4.0による倍精度原本をビルドし、元の係数行列・循環・非貫通残差・渦線力と作用点を取得した。初回40計測条件はCOMMON宣言位置の誤りで力0件となり未達を保持。修正40条件と原本20条件の14判定は成立し、通常出力の不変性と力／モーメントの保存を確認した。ただし [全分量の観察](../../rd/aero/avl-reference/records/avl-core-symmetry-observation.json) ではコア0に不要な横力・ヨーモーメントが残る。一次揚力／ピッチ感度の0.05%条件だけで採用しない。次は直線上／近傍の速度核・鏡映対応を切り分ける。固定20帯・全格子・UVLAS荷重比較・J1は未完了。

### 直線近傍の桁落ちと参照修正候補（2026-10-31）

[切分け・修正回帰](../../rd/aero/avl-reference/records/avl-collinear-kernel-experiment.md) で、コア0の横力・ヨーを有限渦線延長近傍の数値的相殺に特定した。元の力評価点・影響係数を追加取得し、78,336ペアをDecimal80桁と照合、64ペアを160桁で監査した。元循環のまま有限辺外の安定等価式だけで対称性が回復。別Fortran候補では非圧縮・0コア・有限辺外だけを変更し、元の微小コア20条件の通常出力不変、0コア20条件の対称性・残差・荷重和・循環不変を確認した。修正候補IDはAVL3.52/DP/zero-core-exterior-stable-v1。未修正版は不適格のまま保存する。固定20帯・全格子・UVLAS比較・J1は未完了。

### 固定帯の保存と全格子の未達（2026-10-31）

[固定20帯・全格子試験](../../rd/aero/avl-reference/records/avl-fixed-band-convergence-experiment.md) で2翼・各10格子・5迎角の100条件を評価した。283,680寄与の合力／O周りモーメント保存は成立し、事前10条件は8成立・2未達。総量は弦系列の漸近性が不足し、分布はN48→96でも約21.2%の差が残るため参照受入を保留。製造楕円分布では元帯の合力が厳密でも点集中分類だけでN96に7.62%の差が生じた。一定密度の別再構成の良好な結果を元の観測へ置き換えない。次は一次荷重と有限迎角荷重、観測／再構成と保存則を分けて定義・検証する。UVLAS比較・J1は未完了。

### 一次荷重・保存再構成と帯境界適合（2026-10-31）

[分布観測契約](../../rd/aero/avl-reference/spanwise-load-observation-contract.md) のSW-D01〜03で点測度と支持辺＋補正偶力の再構成を区別し、全剛体仕事・合力／モーメントを保存する候補を定義した。一般のANCF荷重伝達へは転用しない。[一次荷重試験](../../rd/aero/avl-reference/records/avl-linear-load-wrench-experiment.md) は有限迎角分解と保存・独立線積分が成立したが、収束3条件が未達。物理翼・観測・合否を変えない [帯境界適合格子](../../rd/aero/avl-reference/records/avl-band-aligned-grid-experiment.md) の別100条件でも点観測と再構成は帯ごとに一致し、総量／分布の漸近性は未達だった。残る循環解・荷重評価の格子依存と求解感度を次に分離する。参照受入・UVLAS比較・J1は未完了。


### 非対称参照系の求解感度（2026-10-31）

[8系の求解感度診断](../../rd/aero/avl-reference/records/avl-solve-sensitivity-experiment.md) は、修正版11条件が成立し、元解／局所直接解／反復改良の荷重差が相対最大1.86e−15だった。検査した格納系の求解精度変更を、格子非漸近性の対策として採用しない。条件数は下からの推定で、上界として消費しない。[精度契約](foundation/linear-solve-accuracy-contract.md) へ対象演算子・推定種別・出力随伴・補償加算順序を追加した。次に係数生成と空間離散化を分けて診断する。参照の収束受入・製品解法の選定は保留。


### 係数生成と比較基底の責任境界（2026-10-31）

[係数照合](../../rd/aero/avl-reference/records/avl-aic-generation-results.md) は7条件成立。検査した係数の生成誤差は最大1.092e−13 m⁻¹、小系全係数の変更による荷重差は最大2.36e−16。大系の部分和を全出力誤差の保証としない。[比較場契約](../../rd/aero/avl-reference/steady-horseshoe-ring-comparison-contract.md) により馬蹄渦の係数と渦輪の一定ポテンシャル差を区別し、同じ添字の直接比較を採用しない。基底変更は作用素・観測重み・キャッシュ世代も含める。契約の代数試験と物理空間の共通細分試験は次段階で、参照受入は保留。


### 基底変換と場比較の実施後判断（2026-10-31）

[比較場試験](../../rd/aero/avl-reference/records/avl-field-overlay-results.md) で強さ・作用素・重みの整合変換と保存則を確認し、重みを変えない負の対照を検出した。共通細分は厳密有理数幾何の版を採用対象の参照診断とし、旧Decimal版を残す。40格子の元荷重は不変。後縁の弦方向4系列のみ漸近候補で、全域P0場と固定平均は未達。前線移動を含むP0差分への方向条件の意味を解析反例で限定し、条件緩和による参照受入は行わない。次に製造投影・符号付きトレース積分の切分けを行う。


### 製造場の既知誤差と符号付き観測（2026-10-31）

[80投影・36トレース比較](../../rd/aero/avl-reference/records/avl-projection-trace-results.md) で、P0方向条件の未達は物理誤差の存在証明にならず、総揚力の差はトレース差の正負相殺にも依存することを確認した。元の未達を取り消さず、精度目標も変えない。[保存観測案](../../rd/aero/finite-wing/finite-wing-observation-error-contract.md) でτ/Kと元点偶力、帯のT₀/T₁/K₀による荷重評価と数学的な上限の使用範囲を定義した。次は恒等式・モーメント保存を検査し、格子差と真の解への誤差評価を区別した参照適格性へ接続する。


### τ/K観測と参照容量の拡張（2026-10-31）

[保存観測](../../rd/aero/avl-reference/records/avl-conservative-observation-results.md) でT₁と元RV偶力を保つ必要性と荷重恒等式、射影後ノルムの場間差上限を確認した。真の解への上界とは区別する。スパン不足に対して [2翼N320](../../rd/aero/avl-reference/records/avl-span320-results.md) を別仕様で実施。NVX/NVMAX以外を変えず、4回帰がバイト一致し、最大約3.42GBで新規2計算が成立した。次の4段階でも収束条件は未達で参照受入を保留。次は対称条件に限定した縮約を全翼データと照合する設計段階であり、全翼解を置換したり製品へ対称性を要求したものではない。

### 対称静的参照の限定使用（2026-10-31）

[全翼／半翼の実施結果](../../rd/aero/avl-reference/records/avl-symmetry-reduction-results.md) により、対象2翼の静的対称部分空間に限って、検証したspacing-pi64版を使用する。初回の幾何不一致と旧版は保持する。旧／新分点の意図した移動には解析包絡を適用し、修正版の全翼／半翼一致しきい値を緩めない。[N640の結果](../../rd/aero/avl-reference/records/avl-symmetry-fine-observation-results.md) は揚力／帯分布の条件付き推定のみ予算内であり、モーメント未評価のため参照全体の採用は保留する。次の細分化は容量・資源仕様を別に定めてから行う。製品一般用途への対称制約や新たな公開APIは導入しない。

### N1280候補の容量だけを拡張する参照経路（2026-10-31）

[事前仕様と経路選択](../../rd/aero/avl-reference/records/avl-span1280-capacity-experiment.md) は外部BLAS・反復法・行列再計算へ変更せず、同じ逐次LUのままNVX／NVMAX16000・NSMAX1000とする。[別版の回帰](../../rd/aero/avl-reference/records/avl-span1280-capacity-results.md) は全41オブジェクト再ビルドと既存8件・66出力のバイト一致が成立した。容量版の割当13.35144GiBと実未知数での算術値を区別し、N1280は新たな資源ゲートと専用ランナーが成立してから実行する。最細格子の精度や参照全体の採用はまだ決定しない。

### N1280の結果による参照判定（2026-10-31）

[N1280評価](../../rd/aero/avl-reference/records/avl-span1280-results.md) では、局所KSMAX不足の初回失敗を保持し、別版の単一定数修正と8回帰後に8計算を完了した。保存・残差・資源・成果物監査が成立しても、スパン次数が約2.25〜2.39から3.22〜3.51へ変わり、両翼5観測量が全て未達。N640で得た揚力／帯分布の条件付き推定をN1280の採用根拠にせず、新基準点の合算予算を全てnullとした。しきい値を緩めず、既存データの次数変動と次の計算経路を設計する。参照全体の採用は保留。

### 次数変動の診断と次の配分感度試験（2026-10-31）

[N1280保存データ診断](../../rd/aero/avl-reference/records/avl-span-order-diagnostic-results.md) から、帯間正負相殺を今回の総量差の説明とせず、翼端の二乗寄与集中を感度試験の動機とする。4つの2項モデルの事後比較は尾部モデルの採用としない。次は総数を保つM12・N1280のtip-rich／tip-poorの両翼4件で配分への応答を調べる。 [事前仕様](../../rd/aero/avl-reference/records/avl-tip-redistribution-experiment.md) と入力は準備済み、計算は未実施。非相似格子からRichardson次数を求めず、M24参照の合算予算nullを維持する。

### 配分感度後の数値感度確認（2026-10-31）

[固定総数の4計算](../../rd/aero/avl-reference/records/avl-tip-redistribution-results.md) は実行・保存・原本とZIPの一致・6場の保存検査を満たした。応答は微小かつ配分増減で非対称であり、翼端だけの原因や精度改善と判断しない。元点モーメントと有限辺・補正偶力の項を別観測として維持する。次に [対称縮約下の求解感度](../../rd/aero/avl-reference/symmetric-reference-solve-sensitivity-design.md) を別計測版で設計する。既存の全翼用感度アダプターの対称拒否を外すだけでは転用せず、W_s=WEの観測写像と元作用素の不変を検証する。参照全体の採用は保留。

### 対称求解感度の結果と独立配分対照（2026-10-31）

[対称求解感度](../../rd/aero/avl-reference/records/avl-symmetric-solve-sensitivity-results.md) の小系・元状態・数値感度検査が成立した。観測重みはWEで半翼へ写し、元の計算結果を保持した局所診断に限定した。試した求解変更はM12配分応答を説明する大きさではなかったため、次は [2因子の独立対照](../../rd/aero/avl-reference/records/avl-allocation-factorial-experiment.md) とする。内側だけ／翼端だけを変更した両翼4入力は準備済みだが未実行。今回の感度比を前進誤差・係数生成誤差・M24尾部の保証にしない。参照全体の受入は保留。

2026-11-05：**2.2採用完了、現在地2.3。** [後流・全荷重仕様](aero/uvlm-wake-adopted-specification.md)、[16群の判定](../../rd/aero/wake-loads/records/uvlm-wake-adoption-results.md)。製品の時間方式・UVLM全体・本体実装開始は未完了。

2026-11-05・2.3更新：[独立参照の実数値・未達・次の終了条件](../../rd/aero/unsteady-reference/records/uvlm-independent-reference-progress.md)。定常の総量・感度は限定成立、非定常の空間収束を追加中。UVLM全体と本体実装開始は未完了。

2026-11-05・第1参照キャンペーン終了：**2.2完了、現在地2.3。** 定常2翼・20帯、非定常上下運動、主翼／尾翼の独立比較が限定成立。[成立範囲・未達・残条件・監査](../../rd/aero/unsteady-reference/records/uvlm-reference-campaign-summary-2026-11-05.md)。2.3全体・2.4・本体実装開始は未完了。
