# PCHIPを支えるprivate有理slotの製品検証

2026-11-19。Issue #166、親 #19。固定製品sourceは `b9d751a726322314a25ab67465b85a88ef4295e1`。[採用契約](../../../specs/aero/pchip-time-curve-contract.md)の整数primitive・有理演算・資源仕様を、可変幅のprivate製品基盤へ移した。公開PCHIP APIは追加していない。係数準備、有理jet、物理入力のissuer・世代・保存復元、有限精度UVLM投影は既存IM03-0の後続である。

## 完了条件との対応

| Issueの行 | 実装と確認した内容 |
|---|---|
| 非確保counterと全call計画 | 固定opcode配列で未約分bit上界・SSA数を数え、carry 1bitと32bitワード丸め後の全payload・workをchecked演算で計画。0cap、等号と1不足、UINT64_MAX overflow、失敗したcountの非変更を試験した |
| 整数・canonical有理 | 符号付きmagnitudeのcarry加減・schoolbook積・restoring商余り・Euclidean GCDと全9opcodeを実装。非canonical入力はGCDを検査して拒否。数値実行と通常の拒否・取消は例外によらずstatus/Resultで伝播する |
| 所有・失敗・再利用 | state/limb/metadataの全payloadを最初の確保前に照合し、成功分を累積課金。非空Workspaceの3確保fault、最終確保からの取消、入れ子確保、診断保持、再入、旧slot、移動、再lease、取消後再利用を確認した |
| 独立算術・費用 | cpp_int/cpp_rationalと手計算を許容差0で照合。各opcodeの実訪問数を採用cost式の独立転記と比較し、全slot数・実payloadも照合した |
| 監査・並列・提出証拠 | 専用processの確保監査、独立Workspaceの1/2/4同時call、下記2変異、対象D/R・4TU静的検査と小入力記録を揃えた。独立レビュー担当がAPI・演算・失敗伝播・日本語報告を通読し、提出を妨げる指摘はなかった |

Workspaceは準備済み領域を所有し、Session・Slot・RationalViewはその借用である。最後のWorkspace ownerはSessionより長く生存する必要がある。epochは同じ生存Workspaceの再leaseで古いslotを拒否するために使い、Workspace破棄後のslot再利用を許可しない。曲線や物理状態の資格を生成する入口はない。

## 参照と変異

正負・零を含む25有理対、INT64_MIN/MAX、1 + 2^-200 + 2^-400、3200/6400/19000bit等の密limbを使った。独立参照側は製品primitiveを呼ばずBoostのcpp_int/cpp_rationalで値を計算する。製品の返却limbも独立にGCD・零の符号・leading zeroを検査するため、参照有理の正規化で不正な返却表現を隠さない。32bit境界、bit/work/slot/byte不足、不正有理、除算零を固定負例にした。

独立レビュー担当が別の一時directoryで最終整数source SHA256 `7ce3807a949b725e003d074ee5fe170a8642856353798eae357e43e2eaa7bbd5` を複製し、既存normal Debug objectとリンクした。制御2caseは1215assertで成功。負の非整数floor補正の符号条件を反転すると1055assert中21件、積のcarry代入を0にすると160assert中4件が失敗した。全3個体はcompile/link成功、Catch exitは0/42/42であり、compile失敗を検出と数えていない。生ログは独立担当の実行記録にあり、ここには結果を転記し、実行時の変異手順に入力source hash照合を加えた再現用producerを保存した。再現用producer自身を再実行した結果とは区別する。

## 構成と監査の範囲

Windows 11、MSVC 19.51.36260.0、C++23、sse2。Debug/Releaseともnormal9case/1352assert、resource5case/68assertが成功した。乱数seedは166。全変更4TUのclang-tidy21.1.1に `-Wsign-compare -Wconversion -Wrange-loop-construct` も渡した。3TUの結果とnormal最終追加試験後の再検査を束ねており、全て成功。`tidy-final.txt` 内のnormal旧指摘はmoved-from仕様試験の除外位置だけであり、最終 `tidy-opcost-final.txt` で解消している。

確保監査はC++ global new/delete群を対象とし、準備済みの正常算術・不正入力・除算零・取消・work拒否で0回だった。C heapそのものを動的に計測したとは主張しない。数値call graphには動的backend、コンテナ拡張、throw/catchを含めず、通常拒否の例外オブジェクト確保も避けた。factoryの所有確保は既存ResourcePlanの値failure境界を使う。

Linux/GCC14・sanitizerはこの記録では未実行で、該当CIへ渡す。新private実装と対象試験だけの変更のため、無関係な全suiteは繰り返さない。既存採用式・R&D・基準結果のbytesは保持した。

## 小入力の工程と資源

[probe](measure.cpp)を固定Release libraryへ `/O2 /MD /std:c++23preview /W4 /WX` でリンク。Intel Core i5-13400F、1caller、各工程10000回、同じ−7/3 + 5/11のfloorを使った。工程計時は参照生成と出力を除き、確保工程には破棄を含む。評価にはlease・入力検査・加算・floor・返却値照合を含む。

[生の計測値](measure.json)ではcount合計1,637,700ns、確保破棄1,630,500ns、評価5,034,000ns。B=8、w=1、4slot、計画payload792B、予約3028stepに対し実229stepだった。工程境界のRSSは4,804,608→4,853,760→4,861,952B。最後の値は同時にprocess生涯peakであり、工程内の瞬間peakやallocator overheadをpayloadから求めた値ではない。比較対象のない小入力記録であり、改善率・製品規模・実務性能を認定しない。

## 編集事故と復元

静的指摘の整形中、PowerShellの単一入れ子配列が文字列へ展開され、未commit新規3ファイルへ文字単位の置換をかけた。既存の退避・置換件数assert規約を守らなかった編集ミスである。hppは独立担当の完全読取り、cppは元編集原文と既存cast訂正、resource試験は元原文と無傷の識別子から復元した。破損写像を再適用して3ファイル全体の一致を確認し、曖昧な字句は元原文と照合した。以後は元byteを退避し各置換が1件であることをassertしてから書き、固定sourceで対象D/Rをやり直した。事故前のobjectや成功数を最終検証の代わりにしていない。

## 再現

`out/issue166-build` はNinja Multi-Config、manifest install OFF、既存issue99 vcpkg prefixを使用した。対象buildは `cmake --build out/issue166-build --config Debug|Release --target uvlas_test_aero_pchip_rational_slots uvlas_test_aero_pchip_rational_slots_resources`。各exeに `--rng-seed 166 --abort --reporter compact` を渡す。変異はMSVC環境で `run-mutations.ps1 -Scratch <新規空directory> -Repo <固定source> -Packages <x64-windows prefix>` を実行する。新規stdout artifactは保存前にLF化し末尾空行を整理した。既存証拠は再整形していない。
