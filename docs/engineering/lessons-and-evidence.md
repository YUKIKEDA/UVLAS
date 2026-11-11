# 今回のR&Dから抽出した知見と根拠

2026-11-06。対象はUVLASの定式化・Python／C++設計試作・独立参照・ビルド基盤・実装前審査。**製品本体を実務運用した経験の記録ではない。** 個別試験の成立範囲を保ち、経験則と厳密な保証を区別して再利用する。

## 進め方で改めるべきこと

| ID | 今回の問題・根拠 | 再利用する判断 |
|---|---|---|
| P-01 | ユーザーは「プロセスを学ぶ」ことではなく実務品質の製品を要求した。[目的](../requirements/product-requirements.md) | 内部で知識を残すことと製品目的を分ける。難しい問題の免除理由に学習・段階公開を使わない |
| P-02 | 要件QAだけでは式・主要I/F・不変条件がなく、実装開始を宣言できなかった。[開始条件](../../rd/planning/preimplementation-plan.md) | 合意事項を実装可能な決定へ変換する工程を独立に置く。開始判定は対象処理の依存閉包で行う |
| P-03 | 主シェル未解決のまま独立領域を先行する案は、核心を後回しにしただけと指摘された。[撤回とv2](../../rd/structure/main-shell/review/structural-progress-control-v2.md) | 独立作業の進捗と核心の解決を分離する。ボトルネックに決定的な識別実験と終了条件を置く |
| P-04 | 150座標案は歪みで共有曲げ空間が減ることを保存核から識別できた。[原因診断](../../rd/structure/main-shell/review/main-shell-root-cause.md) | 新候補は不採用理由を変える構成だけを数える。次数・重み・格子の変更を原因説明の代用にしない |
| P-05 | 多数の履歴文書と現在地表示に古い状態が残り得た。[ロードマップ](../../rd/planning/preimplementation-history.md) | 現在の採用・未達・次作業を一つの入口へ集約。監査ファイル数や文書数で製品進捗を示さない |
| P-06 | 作業4の限定試作は独立した事前仕様書の整備が不十分だったと記録した。[4の終了判断](../../rd/foundation/numerics/records/numerical-foundation-adoption-completion.md) | 実行前の判定と、実行後の解釈を区別。後付け仕様を事前固定済みと言わない。全探索メモの重い文書化は不要 |
| P-07 | 最新の開始判定もIM-00の開始許可であり、本体未着手・実験／性能未認定。[開始審査](../plans/implementation-start-review.md) | 一つのgateの合格を全体へ広げない。今後のレビューでは過去の開始宣言そのものの根拠も検査する |

## 数値・力学・比較の知見

| ID | 得られた知見 | 証拠と成立範囲／再利用時の注意 |
|---|---|---|
| N-01 | 局所の剛体零mode・単要素純曲げが成立しても、共有・支持後の空間が必要な曲率を失う場合がある | [主シェル原因](../../rd/structure/main-shell/review/main-shell-root-cause.md)。特定sectorの解析と特定配置の数値核を、任意格子の証明へ広げない |
| N-02 | 薄肉では残る膜・せん断の零次ひずみが曲げに対して増幅する。固定厚さの収束だけでは不十分 | [歪み薄肉](../../rd/structure/main-shell/q9-plate/records/q9-distorted-plate-results.md)、[原因のスケール](../../rd/structure/main-shell/review/main-shell-root-cause.md)。厚さと格子を独立に変え、誤差増幅を観測する |
| N-03 | 体積・界面・上下面・実法線厚さ・位置全トレースは一組で仕事を満たす必要がある | [局所／大域](../../rd/structure/main-shell/paired-geometry/paired-reference-local-global-contract.md)、[主シェル採用](../specs/structure/main-shell-adopted-specification.md)。材料や場の仮定を変えた案へ元の品質証拠を流用しない |
| N-04 | 有理rank・区間合同変換の保証には対象がある | [主シェルSH-D03と基準](../specs/structure/main-shell-adopted-specification.md)。求積幾何内積の誤差を連続積分全体の保証と呼ばない |
| N-05 | Kqを同じKで解く自己整合と、独立体積／全面荷重のPDE照合は別 | [R4](../../rd/structure/main-shell/paired-geometry/records/paired-reference-r4-results.md)、[R5荷重](../../rd/structure/main-shell/paired-geometry/records/paired-shell-r5-domain-v2-result.json)。製造解を否定せず、何を検出できるかを分ける |
| N-06 | 残差積和だけ高精度にしても、解と補正の保持・返却丸めで精度を失う | [Q9精度診断](../../rd/structure/main-shell/q9-plate/records/q9-precision-diagnostic-experiment.md)、[精度契約](../specs/foundation/linear-solve-accuracy-contract.md)。格納Aの求解精度はAの組立精度を保証しない |
| N-07 | 合力の求解感度と局所荷重／momentの感度は異なる。条件推定には上界・下側推定の区別が必要 | [精度契約](../specs/foundation/linear-solve-accuracy-contract.md)。小残差・推定rcondから全観測の前進誤差を断定しない |
| N-08 | RMとScully型では、同じεでも端点・分母の式が異なる | [Pteraコア診断](../../rd/aero/unsteady-reference/records/ptera-core-model-diagnostic.md)。理想核比較の成立で有限コアの差を消したことにしない |
| N-09 | 有限格子の辺中点荷重、正則化エネルギー、遠方観測は同じ有限値ではない | [UVLM採用と旧未達](../../rd/aero/adoption/records/uvlm-adoption-completion.md)。共通極限の一致と固定格子の恒等性を区別する |
| N-10 | 共有辺の強度差を先に作る評価と渦輪列を先に足す評価は、相殺の丸めが違う | [残差評価政策](../../rd/aero/unsteady-reference/records/high-frequency-residual-evaluation-policy.md)。同じ原作用素か確認し、旧未達・評価範囲を保って精密比較する |
| N-11 | 幾何依存荷重の転置写像だけでなく、その写像の微分も必要 | [構造統合](../specs/structure/structural-adopted-specification.md)、[全荷重](../../rd/structure/integration/records/structural-adoption-results.md)。追従圧力の非対称接線を対称化しない |
| N-12 | 物質成長・固定年齢窓・同じ離散平衡では、初期微分の境界項と係数が異なる | [J1最終追加規約](../specs/coupled/j1-initialization-output-derivative-addendum.md)。輸送座標を変えて物質速度だけを流用すると偽の時間微分になる |
| N-13 | 受理済み荷重は実際の時間法の循環微分を使う。初期整合式を毎stepの出力再定義に使わない | [J1最終追加規約](../specs/coupled/j1-initialization-output-derivative-addendum.md)。物理釣合い加速度と算法差分を別名で保持する |
| N-14 | 次のBDF2に二時刻の循環で足りても、現在の微分出力の復元には足りない | [保存追加項目](../specs/coupled/j1-initialization-output-derivative-addendum.md)、[再開結果](../../rd/coupled/j1/records/coupled-checkpoint-result.json)。現在の受理微分・荷重・ステンシルを保存し、BDF1で埋めない |
| N-15 | 姿勢／chartを変更すると必要な全履歴も変換対象になる | [SO(3)履歴結果](../../rd/coupled/j1/records/so3-time-history-result.json)。現在状態だけの再中心化を同じ再開としない |
| N-16 | 安定性は同じ有限hの平衡・全履歴時間写像から評価する | [J1採用](../specs/coupled/coupled-adopted-specification.md)、[実PG自由後流](../../rd/coupled/j1/records/pg-uvlm-free-step-result.json)。2方向構造部分空間の結果を全翼mode品質へ広げない |
| N-17 | 周波数順はmode対応ではない。近接群・成長率不確かさ・全探索queueが必要 | [探索対照](../../rd/coupled/j1/records/flutter-search-control-result.json)。既知全スペクトルの制御試験は、実翼でのmode被覆の認定とは別 |
| N-18 | 同じ寸法でも幾何・拘束・精度・親世代が違えばcacheは無効 | [N4](../specs/foundation/numerical-foundation-adoption-specification.md)、[再開対照](../../rd/coupled/j1/records/coupled-checkpoint-result.json)。ポインタや配列長で同一性を代替しない |
| N-19 | 有界patchを繰り返す容量と、接続した大域系の分解費用は異なる | [R4容量の限定](../../rd/structure/main-shell/paired-geometry/records/mass-paired-capacity-result.json)、[数値基盤容量](../../rd/foundation/numerics/records/numerical-foundation-adoption-completion.md)。合成100kを製品性能へ読み替えない |
| N-20 | 実験sourceの入手と、比較に使える不確かさ・幾何・構成の認定は別 | [S10診断](../../rd/acceptance/s10/records/s10-reference-metadata-result.json)、[A5](../requirements/practical-acceptance-plan.md)。sampling rateや他コード間ばらつきで較正精度を作らない |
| N-21 | 単位の内部整合が確認できても較正は証明できない。10g／11g等は異なる物の重量かもしれない | [S10診断](../../rd/acceptance/s10/records/s10-reference-metadata-result.json)、[A5](../requirements/practical-acceptance-plan.md)。誤植と断定せず構成対応を確認。sensor xをyから創作しない |
| N-22 | 計画被覆・hash・限定試作の成功は、製品の実装・物理・性能の合格とは別 | [5の終了](../../rd/acceptance/plan/records/practical-acceptance-plan-completion.md)、[開始判定](../plans/implementation-start-review.md)。今回追加したhashツールも科学的合否を判定しない |

## 再利用先

- 原因診断・実験・採否：`solver-rd`。
- 要求から式・時刻・主要I/F・不変条件・ケースを導く：`solver-design-contracts`。
- 実装者から受け取った変更・証拠を審査する：`solver-implementation-review`。
- 常に守る境界：[SE-01〜14](solver-engineering-rules.md)。UVLASへの具体適用：[レビュー・プロファイル](uvlas-review-profile.md)。

同じ問題が実装後にも現れる保証はなく、ここにない問題がない保証もない。新しい実害・反例を得たときに、適用条件と根拠を付けて知識を更新する。個別の例や利用者の今回の選択を、全プロジェクトの絶対規則へ昇格しない。
