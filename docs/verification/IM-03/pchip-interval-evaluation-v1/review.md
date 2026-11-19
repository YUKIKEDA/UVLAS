# private exact PCHIP指定区間評価の検証

2026-11-19、Issue #185、親 #19。固定sourceは `ab210fd589dd9bd8ade5760dd06af50f013e5f2b`、数値cpp SHAは `8b37a6dfd7620f3629409693450b18d9e087956b35831c70070d6f935ccb1665`。[採用契約](../../../specs/aero/pchip-time-curve-contract.md) §2の値・物理時刻一階/二階を、指定した閉区間でexactに評価するprivate部品である。

## 完了条件との対応

| Issueの行 | 実装と確認 |
|---|---|
| 採用式と全要求 | 同じHorner graphで値、要求dt/dttを計算。区間端点・内点、要求0/1/2、全成分を独立参照と許容差0で照合。未要求の導関数はInvalidInputで返し、0を偽装しない |
| countと全call資源 | 同じ式driverをcount/数値で使い、stencil/query import、全slot、未約分幅、32bitワード丸め、結果header/record/limbとcopy/workを確保前に予約。byte/work/slot/input-bit/intermediate-bitの等号と1不足を確認 |
| 所有・失敗 | query・局所u・選択区間・要求値を独立所有。source/Workspace/Plan破棄後、move後のborrowを確認。範囲外・未知要求・inactive・非canonical語を拒否し、元cause/detailsと旧結果を保持 |
| 独立数値検証 | 凍結60sourceの全120区間×u=0,1/3,1の360queryを、全要求orderで比較。新oracleは元標本とslopeからHermite基底と解析的な時間微分で生成し、製品係数を使わない |
| 検証と提出証拠 | 全7確保fault、累積課金、最終allocator取消、再入/再利用、1/2/4callerのcanonical語とcounter一致、D/R・静的・変異・最小工程・独立レビューを完了 |

局所uは $$ u=(t-t_i)/(t_{i+1}-t_i) $$ であり、Periodicの正規化物理時刻ではない。区間検索、outside政策、side/C1、issuer/世代と公開ExactPchipJet、有限UVLM投影は既存IM03-0の次段である。本結果はそれらの資格を生成せず、親曲線を保持する公開jetとは別型である。

## 参照と実行

元fixture SHAは `c8b8afb2b25acd691eef0060b8a38afe6dedf56eb9e6bbe8e10cd877c903c1e8`。新fixture SHAは `2d08ac65e4485272a82cdc18aa50f591b8daa5d40dce3b80829847f8c23ca563`。独立担当がline/cubic/constantの閉形式9件、別directoryへのbyte再生成、原本hash不一致と既存出力の拒否を確認した。producerはRuff0.17.4/formatとPyright1.1.408 strictが成功した。

Windows11、MSVC19.51.36260.0、C++23、sse2。Debug/Releaseとも評価normal4case・38585assert、評価resource filter2case・272assert成功。seed185。製品/normal2TUのWindows clang-tidy21.1.1と、独立担当のresource TU解析が成功した。未変更の先行係数suiteは実行せず、既存監査TUへ新規2caseだけ追加してフィルタした。Linux/GCC14/sanitizerはCIへ渡し、成功を先書きしない。

無変異controlは通常Debug成功に含まれるHermite case38450assertである。独立担当が同cppの私有複製で一階微分の1/hを除き、MSVC Debugでcompile/link成功後、38450assert中242件が失敗した（Catch exit42）。例はtwo_points/base、区間0、u=0のdtが4、exact参照2。分離リンクした別controlは作っておらず、その実行結果とは区別する。観測は実行時ツールからindependent-observations.jsonへ転記した。

typed係数は既に不変のcanonical ownerで、queryはarenaへのcopy後にGCD検査する。数値loop中の追加C++ newはなく、所有結果を含む全callの計画確保は7回。C heap全体を動的監査したとは主張しない。既存slot契約のinput capは固定literal2/3/6も含み、Secondの境界3bitを検査した。独立担当が採用式、所有、資源と本報告を通読した。

## 小入力と再現

固定Release libraryへmeasure.cppを/O2・/fp:strictでリンク。Intel Core i5-13400F、1caller、3節点1成分の係数を準備済み、区間0・時刻1/2・Secondを1000回。count合計1,512,600ns、全評価/破棄8,580,200ns。後者は再count・7確保・数値・copyを含み、差を純数値時間とは呼ばない。計画2280B、31slot、252164stepに対し実31slot・1069step、B25/w1。RSSは4,919,296→4,927,488→4,927,488Bで、最後はprocess生涯peak。allocator overhead・工程内瞬間peak・実務規模性能の認定ではない。build_infoはab210fd、uncommitted=false。

既存C3専用Ninja Multi-Config buildを使用。対象normalは `uvlas_test_aero_pchip_interval`、資源は既存 `uvlas_test_aero_pchip_coefficients_resources` に `Private PCHIP interval*` filterを渡す。seed185/compact reporter。新stdoutは既存redactorで機械pathを伏せ、前後hashとLF化をpublication.jsonへ記録した。既存採用仕様・証拠は保持した。
