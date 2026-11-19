# PG分割の辺出現上界：製品検証

2026-11-18。Issue #146、SH-D07/N4-02/N4-08。採用契約pg-partition-incidence-budget-contractの全補正項を、製品資源計画へ移した。製品の幾何・exact述語・返却順・bit上限・外部推定式は変更していない。

固定製品commitは6049a768a458abe935c9d9bad7bf3a09a084d7f2、旧版は89c3923c0d0648910ad700390f7d04d77eec8e40。Windows/MSVCの既存依存prefixを使い、影響する通常・資源の2targetだけをDebug/Releaseで実行した。両構成とも2/2成功。通常8例30,627 assertions、資源6例189 assertions。既存の全n独立面積/moment・全層/全外面/trace、所有、全quota等号/1不足、overflow、確保失敗、取消、再入、別Workspace並列を含む。

新しいn4 U形・2層・名付き部分境界1・refinement1の計画は、入力件数を独立に整数代入したliteralと照合した。clip容量8,712、numeric 228,633,355、topology 184,228,797。最大値を実測から縮めていない。新例でも予測等号で準備が成功し、numeric/topology/owned/externalの各1不足が確保前にResourceExceededとなり、旧内容IDが不変だった。

## 固定コピーの旧新比較

専用runnerは各版をgit archiveした隔離コピーで製品TUをコンパイルし、固定した同じprobeと対象試験objectへリンクした。元worktreeのファイルは変更しない。同じn4入力に対する比較は次のとおり。

| 指標 | 旧版 | 新版 |
|---|---:|---:|
| 所有累積payload予測 B | 205,783,717,415 | 14,655,586,727 |
| 所有累積payload実績 B | 33,873,559 | 30,025,495 |
| Boost外部予測 B | 19,978,141,584 | 19,978,141,584 |
| numeric予測 | 3,214,950,271 | 228,633,355 |
| topology予測 | 464,517,909 | 184,228,797 |
| 実numeric / topology | 178,891 / 37,664 | 178,891 / 37,664 |
| 48 GiBの所有Plan | 拒否・確保0 B | 成功 |

旧版の高cap成功結果と新版の成功結果を、全1,208 recordのcanonical snapshotでexact比較した。頂点、cell面積/活動/成分、全片/辺/接続/由来、三角形、全体積/水平面/側面/trace/接触、元graph全係数とcell、元source全variant、厚さ全値/三角形、層ID/端点、名前と内容IDを走査し、全て一致した。forecast・timing・確保量は数学結果の比較から除く。canonical snapshot SHA-256は `5ed39b823027928161d6dca3a747c55af6fa0f50ed356c83364aca631a70dfc4`。

この48 GiBは所有Planのcapであり、RSS hard capではない。累積payloadと同時Boost予測を区別する。プロセスRSS・全モデル準備・実務性能の認定には使わない。R-02と局所rankの資源方式は後続の既存gateに残る。

## 変異と訂正

固定コピー上で正常controlはcompile/link/test成功。厳密異符号判定を非厳密へ変えて零長辺を生成、clip容量をMへ過少化、候補点走査11VPを除去、原子辺sort項を除去する4変異は、全てcompile/link成功後にtest exit42で検出した。過少見積りに対し余剰quotaを使って成功するだけでは合格にせず、独立計画literalも検査する。

snapshot拡張前の既存試験は、refinement版だけを変更して幾何結果が同じことを調べていた。全所有snapshotへ拡張した初回Debugでは、版/内容IDが異なるためその比較が1件失敗した。該当試験を幾何同値と所有内容の相違の二つへ明示して修正し、通常targetだけを再検証した。最終固定sourceではD/R両targetが成功している。

実行用producerとrunner、隔離のcompile/link/testログ、旧新の全snapshotをこのdirectoryへ保存した。製品全体のsuiteは変更外なので実行していない。Linux実行はこの環境では未実行でCIの対象とし、CIは待機しない。

対象3TUはDebug各1entryの専用compile DBを作り、clang-tidy 21.1.1でexit0。第三者headerの警告194,362件と既存NOLINT150件は抑制された。変更5ファイルはclang-format 21.1.2へ適合した。日本語lintは0件。独立レビュー担当が計画式と全所有snapshotの走査範囲を通読した。
