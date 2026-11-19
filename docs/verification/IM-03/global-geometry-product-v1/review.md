# 大域幾何・TE所有の製品検証

2026-11-18。Issue #147、親 #19。最終製品sourceは693b963。採用した[開始契約](../../../specs/aero/global-geometry-te-contract.md)と[二回走査・逆aliasの補正](../../../specs/aero/global-geometry-te-review-addendum.md)を実装した。物理写像providerと状態間経路の資格、循環、荷重、後流の実放出・commitは生成しない。

## 完了条件との対応

| Issueの行 | 製品と確認した内容 |
|---|---|
| 所有topology・TE | 元OwnedGeometryを所有共有し、全face role、直接edge alias、元端点と正規基底、D1D2、T、Bound継ぎ目/quad対角拒否、epoch付き唯一キーを実装。raw alias ID順とowner順が逆になる点・辺aliasを通した |
| 全面対・候補と資源 | 閉AABB balanced BVHの二回走査、2B先行予約、有限三角形の非平行/共面/点/線分/面積交差と許容集合を実装。候補はexact重複除去。K=0/3で二回分の価格を独立に照合し、work/bytes等号と1不足を検査 |
| 全alias運動学 | 全CSR列、明示時間項、共通DOFの元語・単位・基底・版を所有。exact作用を先に無確保countして、実際に必要なpacked limb帯だけを確保。三低位帯、最大/最小元語、binary64範囲外の作用、各列/gt/速度assertion不一致、不正CSRを検査 |
| 失敗・所有・再入 | 全確保fault、workspace再入、早期/最終allocator取消、元diagnostics、親不変、移動元拒否、入力/Plan/workspace破棄後の寿命を検査。独立workspaceからの1/2/4同時callで全検査counterが一致 |
| 独立検証・記録 | 下記のFraction参照、固定source変異、対象D/R、静的検査、最小工程計測とhash/traceを記録。rootが設計・数値body・資源計画を独立通読した |

Windows 11、MSVC 19.51.36260.0、C++23、sse2、`/fp:strict`。固定製品source 693b963のDebug/Releaseともnormalは12件4921 assertion、resourceは6件126 assertionが成功した。生成したbuild identityは693b963、dirty=0。Release監査構成の既存 `/DNDEBUG` と `/UNDEBUG` によるD9025はbuild logに出るが、buildはexit0である。

## 参照・許容差・変異

凍結Fraction oracleのSHAを固定し、重心座標の非負多面体から求めた12例を製品fixtureへ変換した。変換後の座標はbinary64で厳密に表せることを確認。別directoryへの再生成はSHA e1d32d17d82ce3042f8121a58129e051ba916234a5af67ab8afc4ecd72cd1d4cで一致した。製品の平面切断/共面投影コードを参照生成に使っていない。

全12例×6×6置換×4配置（元尺度、2^-900、2^900、共通原点2^900）を、許容差0で照合した。共有頂点を越す交差、同端点別EdgeId、重複面、包含、非零の低位gapを含む。固定非平面quad、三面共有辺、180度回転も製品で確認した。零面積は先行するA1 factory自身が拒否するため、その入口条件を再利用した。

固定数値source16d46c6の私有byte複製でcontrolと10変異をcompile/linkした。controlは全てexit0。二回目の予約欠落、逆alias基底省略、T符号、内部TE許可、共有点による全面対skip、開箱、J値比較省略、gt比較省略、速度assertion無視、epoch無視は、全てcompile/link成功後にCatchのassertion失敗42となった。compile失敗を感度の証拠に数えていない。

旧producerは最後の期待Catch42をshellへ残す。出力summaryは11件完了しており、新しい`run-verified.ps1`が名前・整数型・各終了状態・ログ存在を検査してexit0を返す。今回は`-ValidateExisting`で完了済み結果を検査し、10変異を重複実行していない。

## 固定後の限定追補

最終shared controlblock allocationからstopを要求する私有probeでは、16d46c6が成功を返し、3 assertion中1件が失敗した。f501fbfでallocator直後の`run.step(0)`を加え、公開前に取消を検査した。追加resource回帰は固定製品source 693b963で成功し、leaseの回復も確認した。

9491a3eではnoexceptのword比較を`get_if`に変えた。693b963では、Linuxにも必要な符号変換警告で実診断されたBVH極値/coordinateのordinal型をA1の`core::Index`へ合わせた。元の数値式・価格・成功値は変えていない。後者は有効なordinalが非負であることをA1が保証し、64bitのサイズ/値も同じである。10変異は16dの数値bodyの証拠として明示して再利用し、追補は該当D/Rと静的診断で確認した。

対象6TUの静的検査ではnoexceptの指摘を記録し、変更したglobal TUで解消を確認した。追加取消resource TUも確認した。さらに符号変換の出た3TUは必須warningを有効にした限定診断でexit0となった。これはWindows上のClang 21.1.1であり、GCC実行の代用とは主張しない。LinuxはPRの該当CIへ渡し、完了を待たない。

静的検査の最後に、取消回帰内のローカル`Setup`だけをconstへ修正した。呼ぶ関数はいずれも既にconst受取りであり、生成する入力・処理・値は変わらない。実行結果は693b963のD/Rを再利用し、この変更はresource TUの静的検査exit0で確認した。ライブラリの数値sourceを変えていない。

## 工程と費用の境界

9491a3eのclean Release、1 worker、2 facets/5 original points/4 physical owners/2 TE/1 DOFで各APIを16回測定した。N4のStageRecorderとprocess RSSを使った。結果は`measure-release.json`にある。

| 処理 | 1 callの累積payload B | 予約primitive単位 | 消費primitive単位 |
|---|---:|---:|---:|
| topology | 914 | 31,429,164 | 13,411,431 |
| contact | 704 | 4,203,422,364 | 90,406,969 |
| kinematics | 2,016 | 86,621,225 | 9,457,769 |

primitive単位は固定limb容量を使った課金単位であり、CPU命令数ではない。opaque controlblock等は明示したestimateで受け入れ、未知の厳密上界を0扱いしない。必要なproofのないRequireProvenBoundは拒否する。N4のpayload、allocator metadata、process lifetime peak RSSは別の量である。

この小入力から実規模の性能を外挿しない。BVHの事前上界は二乗であり、R-03の実務入力/時間/RSS受入はIM03-0の後続で測る。内部面対の並列化や速度改善は主張せず、今回の並列確認は独立callの決定性と親の不変性である。

## 再現

専用CMake領域で`uvlas_test_aero_global_geometry`と`uvlas_test_aero_global_geometry_resources`をD/Rでbuildして実行する。変異sourceは16d46c6に固定する。対応sourceのcheckoutからMSVC環境内で`run-verified.ps1 -Root <checkout> -Output <new-directory> -SourceCommit 16d46c6...`を実行する。既存出力への上書きは拒否する。

controllerはこの記録の版を使い、Rootに16d46c6の別checkoutを指定する。事前にそのRootの`out/issue147-build`でDebug normalをbuildする。`out/issue147-tidy-db/compile_commands.json`は、同buildのcompile databaseから4製品cppについて`uvlas_aero.dir`かつDebugのentryを一つずつ選んだJSON配列である。この選択はMSVCのcommand文字列と依存prefixをそのまま保持する。今回の測定probeも、そのRelease commandで一時compile/linkした。製品公開targetは増やしていない。

最終修正が新しいaero処理とsealed friendだけなので、無関係な構造・IO・P1核suiteは再実行していない。旧仕様・R&D・既存証拠はbyte保持した。新しい記録はLFで、属性の例外を足していない。
