# UVLASの実装レビュー・プロファイル

2026-11-06。汎用Skillをこの製品へ適用するときの入口。ユーザーは製品本体の実装を別の担当者へ依頼し、このアシスタントにレビューを依頼する予定。今回の知識整理では本体を実装しない。将来の作業はその時の依頼範囲に従う。

## 決定の種類

ユーザー確定の製品目的・品質・範囲と、アシスタントが根拠を付けて採用した技術仕様を区別する。後者も無条件に正しいとはしない。矛盾や反例が見つかれば仕様欠陥として起票し、根拠・影響・必要な改訂をレビューする。旧文書の不備を継承するために正しい修正を拒否しない。

| 項目 | 正本／扱い |
|---|---|
| 目的・範囲・利用経路 | [製品要求](../requirements/product-requirements.md)。実務品質、静的・非線形過渡・変形平衡周りflutter。学習目的へ置換しない |
| ユーザー決定の履歴 | [開始時の決定](../requirements/project-start.md)。C++23、両OS、CLI／Python、Nastran限定＋Gmsh、再開、配布、Apache-2.0等 |
| 精度・実験証拠の基準 | [受入精度](../requirements/acceptance-criteria.md)。各出力量の基準と数値誤差の合計予算。新規Ruleで勝手に緩めない |
| 構造採用 | [構造統合](../specs/structure/structural-adopted-specification.md)、[主シェル](../specs/structure/main-shell-adopted-specification.md)、[梁](../specs/structure/beam-adopted-specification.md)。履歴のQ9／MMS等を既定へ戻さない |
| 空力採用 | [UVLM統合](../specs/aero/uvlm-integrated-specification.md)。有限核・観測・W-P1・全荷重・共有／TE所有・誤差認定 |
| 時間・強連成・安定性 | [J1統合](../specs/coupled/coupled-adopted-specification.md)。輸送別初期微分と受理後荷重は[最終追加規約](../specs/coupled/j1-initialization-output-derivative-addendum.md)を優先 |
| 型・数値・失敗 | [N4](../specs/foundation/numerical-foundation-adoption-specification.md)、[IO4](../specs/foundation/io-execution-adoption-specification.md)。旧案の密全域nullspaceやBDF1 accessorを流用しない |
| 受入と実装順 | [A5](../requirements/practical-acceptance-plan.md)、[開始審査IM-00〜08／R-01〜07](../plans/implementation-start-review.md)。対象パッケージの依存と未受入を確認 |

重要な製品規約は既に文書化している。ただし「規約がある」「開始審査がある」ことは、その内容に矛盾がない保証でも全パッケージの無条件承認でもない。開始審査の対象は最初のIM-00。レビュー時点のコードと製品検証を確認して後続へ進める。

## 優先的に探す回帰

1. PGの参照係数と全3D材料・一定シェル質量・全位置トレースを、別の要素／応力仮定へ無断で変えていないか。接続後の許容空間を検査する。
2. 可変法線・厚さ・層・面測度と双対荷重の定義、全幾何／材料接線が揃うか。局所Gram行列の保証を連続全体へ広げないか。
3. RM核と参照コードの別核を同じεだから等価と扱わないか。共有辺の循環差・TE所有・荷重の物理IDが一致するか。
4. J1のr/u・Γ・free Xを同じ試行で評価するか。新生点、物質／年齢輸送、初期付加質量の政策を区別するか。
5. 受理後のD_hΓを保存し、別初期化・BDF1で荷重を再生成しないか。D_hvと物理釣合い加速度を区別するか。
6. 全履歴Pの固定点・微分・固有modeを同じ離散平衡で評価するか。少数modeや周波数順で全安定を断定しないか。
7. hi/lo保持の出口、原系残差、Python view・Workspace・世代、全履歴checkpointと一回確定を壊していないか。
8. S10の不確かさ・sensor位置・構成・marker対応を推測で埋めないか。合成5000／100k試作を実PG・全自由後流の性能合格にしないか。

## 受け取るレビュー材料

対象版／差分、IM-ID、変更した式・契約、検証した製品経路とfixture、実行条件・結果、失敗／未確認、性能内訳、採用仕様からの逸脱を求める。欠けていてもコードから調べられる部分は先に進める。

仕様・実装・検証参照が互いに同じ誤りを共有していないかを確認する。とくに19ケース／37要求の計画被覆、1072ソース／687成果物の監査件数は、製品の科学的な合格数ではない。現時点の本体実装済み・実験認定・製品性能を推定する数字にしない。

## 規約改訂

問題箇所、最小反例、変える数式／時刻／状態／精度、影響する既存結果と再検証を明記する。採否済みの証拠ファイルは保持し、新仕様・訂正・新結果を別版で追加する。現状表示の入口はREADMEとロードマップへ同期する。
