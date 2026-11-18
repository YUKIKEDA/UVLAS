# PG多角形と全積分分割の製品検証

Issue #125、親 #17。数値sourceは6369a0b1227bf31db246f00f982a98bc7e2eff64。採用した[共通分割契約](../../../specs/structure/pg-integration-partition-contract.md)の製品適合を扱う。物理基底・求積・rankや、IM-01全体の合格ではない。

## 元の領域から全厚さの同じ分割を返す

元support直線を整数canonical化し、元線対のCramer式から有理交点を作る。全half-cell・厚さ三角形・名付き境界と追加区切りを共通細分し、正面積の凸片、cell内成分、共有原子辺、層別体積、上下面の両material trace、物理外周側面を所有する。内部補助辺を物理外周へ露出しない。領域外の頂点が全正方形の交点表に残っても、積分する三角形と体積は実領域内だけである。

公開getterが返すactivity spanは報告用であり、資格を持つのは専用factoryだけが作るPgIntegrationPartition本体である。元graph・領域・厚さ・層・名前・refinementの実内容も照合し、digestだけで互換としない。

## 独立な全cell照合と失敗系列を通した

Windows 11、MSVC 19.51.36260.0、C++23、/fp:strict、sse2。固定fixtureは4領域×n=1/2/4の12例であり、全336cellの面積と一次moment2成分、計1008値を許容差0で照合した。製品のsupport細分に対し、参照は垂直断面とexact Simpson積分を使う。手計算した全領域momentとも一致した。 $$ 10^{-400} $$ の正面積、1cell内の2つの非連結成分、異なるgradientの連続P1厚さ、T-junction、部分名付き辺、refinementの版、全層界面両側、零接触を含む。

| 構成 | normal | resource | CTest |
|---|---|---|---|
| Debug | 7 case、30,533 assertion | 5 case、169 assertion | 2/2成功、3.53秒 |
| Release | 同じ7 case、30,533 assertion | 同じ5 case、169 assertion | 2/2成功、0.48秒 |

幾何不適格・layer/ID/shape違反、全logical quotaとbit/work/external/payloadの同値・1不足、Proven要求の拒否、checked overflow、9選択割当点とBoost外部割当の故障、allocator再入、開始前と初回確保後の取消を検査した。所有割当はResourcePlanの予約を伴い、それ以外の本経路の割当が固定Boostであることをsourceと照合して外部故障を注入した。旧結果、元入力、lease解放も確認した。

並列性はsource identityだけで判定せず、実返却の頂点・面積・成分・全topology/traceをcanonicalに直列化し、直列・独立2workerでbyte一致を確認した。計測Enabled/Disabledでも同じ数学payloadとなる。

隔離controlは固定commitのbyte複製を使い、normal/auditの両経路でcompile/link/testが成功した。activity、H共辺連続性、内部辺の外面化、層upper側、微小正面積の閾値化、refinement版hash、quota admissionの7変異は、全てcompile/link成功後にCatchの失敗で検出した。build不能やloader失敗を検出に数えていない。v1 harnessはsystem検索libをbuild内ファイルと誤認し、compile前に停止した。元producerを保持し、v2で修正した。

## 実測33.9MBと最悪上界205.8GBを区別する

単独のRelease製品prepareを、n=4、8頂点の凹領域、厚さ2三角形、2層、部分境界1、追加区切り1について3外側×3内側の9回測った。i5-13400F、logical CPU16、物理メモリ51,343,966,208B。独立したCLIの--build-infoは同じ固定sourceとuncommitted_changes=falseを返した。同時負荷を統制した実務性能の認定ではない。

wholeは39.733〜40.653ms。最初の標本では、元入力0.211ms、support/交点0.468ms、平面分割30.814ms、trace等2.858msであり、残りは計画・確保等である。主な工程は平面分割だった。実返却は135頂点、31片、57三角形、83辺、114体積、228水平面、94側面となった。

| 量 | byte | 意味 |
|---|---|---|
| actual累積payload | 33,873,559 | 成功したResourcePlan確保の累積。live heapとは異なる |
| owned forecast | 205,783,717,415 | 全入力最悪loop×scalar storage等による保守的累積上界 |
| external forecast | 19,978,141,584 | 固定Boostの同時limb payload推定。allocator/RSSのhard capではない |
| OS lifetime peak（最初の標本） | 32,075,776 | process開始後のRSS最大。call単独peakではない |

予測と実測の差は、契約5節の全入力最悪上界を維持したためである。例えば最終原子化の全頂点照合を数える11VPを、各候補のscalar storageまで累積する。一方、この実行のnumeric opcodeは178,891、保持有理数の最大bitは13だった。実測を使って恣意的に上界を縮めたり、quotaを緩めたりしていない。都度同次化・符号再評価・接触の両端比較、文字列/内容hashのbyte仕事と、moveで順序を確定するhandle配列を実装の計画へ加算した。

このforecastを厳しいbyte budgetへ適用すると、実際の確保量が収まる場合も事前拒否され得る。48GiBの全体実務予算を満たす一般準備の認定ではない。準備費用・疎graph/fill・実務規模は、既存IM01-0のR-02 gateで判定する。#17は継続する。

## 後続へ渡す資格を限定する

このPRは、表現された有限有理領域の分割・活動根拠、所有・失敗・資源保護の適合までを主張する。任意cutの物理求積精度、法線管、活動空間rank、係数条件付け、PreparedPatch、物理M/Kは未認定である。次のPG作業は既存の作業分解にある活動空間rankと参照係数・幾何・求積である。LinuxとsanitizerはCIへ渡し、未実行を成功と書かない。
