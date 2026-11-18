# 公開有理化の確保前計画

PR #154の公開`materialize_rational`について、converterから独立に呼んでも全工程を確保前に計画するよう訂正した。従来は`1e100 m`の値上限が64または332bitでも、字句bufferの5bytesを課金した後に拒否していた。旧版の実検出は保存し、訂正版の成功と区別する。

元の十進数または有理数の両成分、単位のparse、冪のseed・逆数・各積／二乗、最後の積を一つの行程として走査する。既存exact backendのadmitをそのまま呼び、代数work、backend固有のguard limbs、capacity・GCD・Karatsubaを含む外部容量を照合する。分母parse中の分子、単位処理中の元結果、seed処理中のbase、二乗中の冪結果も同時生存量へ加える。

通常はbit上界と外部容量の保守的な見積もりで認定する。上界が要求の境界をまたぐ場合だけ、allocatorを持たないchecked 1024bit整数で同じ行程の値と約分を証明する。これにより`1 m`の値上限1bit、`1e100 m`の333bit、逆Gramの10bitを上界だけで誤拒否しない。固定証明に収まらず、保守上界でも認定できない境界は確保前の`ResourceExceeded`とする。この判定限界は公開APIにも明記した。より大きな値も十分な保守上界内なら処理し、1024bitを有理化API全体の値上限にはしていない。

成功時に生じるcanonical文字容量をC、非零RationalStorageの個数をKとして、累積payloadは次の量を最初の確保前にPlanの残量へ照合する。

$$
P_{\mathrm{cum}}=C+K\,\mathrm{sizeof}(\mathrm{RationalStorage}).
$$

converter側には固定証明の名前付きframeと追加workを計上した。`certificate_stack_bytes`はC++値payloadの保守的な容量であり、process stackやRSSの実測値ではない。外部limb容量、同時scratch、成功確保の累積量を同じ量として扱わない。

固定source `4ff65ef` のWindows Debug/Releaseで、SI正常18cases・39,040assertionsと資源6cases・3,940assertionsが成功した。公開入口の回帰は、値・backend guard・字句・外部容量・累積Planの不足および後段の単位上限不足について、監査したglobal newの回数とPlanの追加課金がともに0であることを確認する。guardの境界は32/64bitのbackend limb幅を既存exact実装から取得し、固定した576bitを両環境へ当てはめない。

その後の取消handleへの限定lint抑制はコメントだけの変更である。実MSVCの`request_stop()`が非constであることを確認し、libstdc++のconst宣言との差を説明した。数値経路は再実行せず、そのTUのコンパイルと静的検査を通した。Linux CI成功や全体ソルバーの数値認定は、この記録では主張しない。
