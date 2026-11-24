# PG 局所独立化の開始審査

[局所独立化の設計 Issue](https://github.com/YUKIKEDA/UVLAS/issues/133) の完了条件を、採用した契約と製品成果に照らした。初回審査で止めた累積資源方式の判断と証拠を保存し、その後に採用した有限資源方式で実装を開始できると判定する。

| 終了条件 | 確認した内容 |
| --- | --- |
| 入力実体、行キー、独立列、全再構成、並進・C・アンカー | 同じグラフ実値と活動域を要求する。正面積セルの係数行、必須場の優先、所有する B・T・S・W・Z と元写像 G・I・J を定めた。全背景のアンカーは活動域と別に評価する。[有限資源の準備契約](../../specs/structure/pg-rank-bounded-preparation-contract.md) |
| exact 消去、順序、切断、幅、失敗 | 順序付きの非零小行列式で独立選択を証明する。従属列には全元行の厳密零と、後から採った座標の零を要求する。固定幅整数の上限、再検証時の所有値保持、取消・再入・予算超過の診断を定めた。[返却、幅、再検証](../../specs/structure/pg-rank-bounded-preparation-contract.md) |
| 全 call の資源 | 同一 Plan の候補、因子、証明、返却に必要な領域を確保前に検査する。固定幅と最大4素数・128桁で作業量を囲う。入力生成で借りる外部 Exact 資源は Estimate の範囲に置き、未証明の Proven を受け付けない。[事前計画](../../specs/structure/pg-rank-bounded-preparation-contract.md#事前計画と処理上限) [入力との同時保持](pg-rank-resource-coresident/verification.md) |
| 独立計算と負対照 | Fraction による n=1/2/4 の9系列で、全活動・アンカー外・アンカー帯の全原列を照合した。極小の非零成分を含む6変異を検出し、独立レビューが180行列と467従属写像を追加照合した。[設計検証](pg-local-rank-design/verification.md) |
| 旧 n=4 U 字、48 GiB、指標の区別 | 旧分割の 205,783,717,415 B は事前見積りであり、同じ制限の確保前拒否を再現した。実累積課金、同時領域の見積り、プロセスピークを別に記録した。新しい有限資源方式は別の証拠で確認した。[旧分割と実務容量の境界](../../specs/structure/pg-local-rank-contract.md) [有限資源方式の製品検証](pg-rank-bounded-preparation/verification.md) |
| 独立レビューと開始判定 | 設計プローブの別ツリー再実行、[局所 rank の数学契約](https://github.com/YUKIKEDA/UVLAS/pull/141)、[有限資源方式の製品化](https://github.com/YUKIKEDA/UVLAS/issues/264) を照合した。最終の製品検証は実辞書と全原列を用いる。既存の凍結証拠を保持した。[入力の再検証](pg-rank-bounded-preparation/validation-verification.md) |

判定の対象は局所係数空間を作る工程である。R-02 の実規模と IM-01 の終了審査には、物理幾何・条件付け・全厚拘束・疎チャートを含む独立した製品証拠を使う。
