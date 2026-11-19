# private PCHIP節点locatorの検証

2026-11-19、Issue #194、親 #19。固定sourceは `5ce67fa2b1ed8753aed292a8c7c814c80971b1d7`、cpp SHAは `d01716b024cb202f13b5909702630f16a5aa0d2891607e3a32fcf9891f0f9631`。既存の所有係数から閉定義域内をexact検索し、左右または唯一の内側候補を既存IntervalRequestへ渡すprivate部品である。

| 完了条件 | 確認した内容 |
|---|---|
| 検索・候補 | 区間内は左右とも同一区間、内部節点はj−1/j、端点は存在する内側だけ。不存在側・内部節点の唯一要求・未知side・範囲外をInvalidInputで拒否 |
| 全call資源 | source時刻の最大語幅を不確保走査し、queryと最大probe分のslot/byte/workをchecked予約。heapは準備済み有理arenaの3確保だけで、結果は借用なしのスカラー値 |
| 独立正負・接続 | 凍結360queryと全節点/間隙、2の冪前後の節点数を独立cpp_rational線形検索で許容差0比較。既存IntervalRequestへ左右候補を渡し、内部節点で値1/一階3/2は一致、二階0/4を平均せず保持 |
| 故障・決定性 | 3確保全fault、累積課金、各cap等号/1不足、最終取消、再入/再利用、旧値、所有source消滅後、1/2/4callerの結果/counter一致を確認 |
| 対象検証 | Debug normal5case/6466assert、独立resource filter2case/96assert、変更3TUのWindows clang-tidy、意味変異、独立通読とhash/traceを記録 |

検索の比較上限は $$ \operatorname{bit\_width}(n)+3 $$ （範囲2＋二分探索＋一致1）である。一方、容量の決定には全節点の事前走査O(n)を行う。全callを対数時間とは主張しない。全visitは取消と予算を検査し、失敗は元のcause/detailsを値保持する。結果はsourceの本人性、cache identity、親曲線や連続性の資格を持たない。

独立oracleは元fixtureのexact時刻をBoost cpp_rationalで線形走査し、製品の比較/探索を使わない。低位語・大原点・正の時間尺度を含む既存360queryを再利用した。先行係数suiteは回さず、係数作成は入力setupのみ。境界policy、Hold/Periodic、Automatic-D2/C1、公開curve/jetは既存IM03-0の別単位に残る。

Windows11、MSVC19.51.36260.0、C++23。通常Debugの全5caseがcontrolである。独立担当は同cppの私有複製で二分探索の `order<=0` を `<0` に変え、compile/link成功後、線形検索caseの14assert中5件が失敗した（Catch42、seed194）。例はtwo_points/baseの時刻0が節点から外れる誤分類。分離リンクの別controlは実施しておらず、実ツール観測をindependent-observations.jsonへ転記した。

静的検査はWindows clang-tidy21.1.1。製品/normal2TUと独立担当のresource TUが成功し、数値式・予算・候補も独立通読した。新規stdoutだけを既存redactorで伏字/LF化し、前後hashをpublication.jsonへ保存。旧仕様・fixture・証拠を保持した。

実行対象は `uvlas_test_aero_pchip_locator --rng-seed 194 --abort --reporter compact` と、既存資源exeの `Private PCHIP locator* --rng-seed 194 --reporter compact`。必要なDebug/静的/正負/変異に絞り、旧suite・ローカルRelease・追加bench・CI待機は行っていない。LinuxとReleaseの新規成功や実務規模性能は主張しない。
