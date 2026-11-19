# private exact PCHIP係数準備の検証

2026-11-19、Issue #176、親 #19。固定製品sourceは `ba3e8159c641030dad6ef4ee3a58fc906b742b55`。[採用契約](../../../specs/aero/pchip-time-curve-contract.md) §2/5/6の全slope枝と無次元u係数を実装した。数学的な時刻・元値・係数を所有するprivate部品であり、政策/C1・events・issuer/世代・保存復元・jet・有限UVLM投影の資格は生成しない。これらは既存IM03-0の後続である。

## 完了条件との対応

| Issueの行 | 確認した内容 |
|---|---|
| 全式・shape | n>=2、m>0、厳密増加時刻、2点・両端・内部の全slope枝、係数のinterval/component/power順を実装。60例の全係数と独立Hermite端値/端傾きを許容差0で照合 |
| 不確保countとadmission | 同じ式driverでcount/数値を実行。分岐はopcode回数・bit・返却分子/分母幅をmax合流し、非Compare回数の和をSの保守上界とする。prefixを二重計上しない67stepの手計算例、nm/R/2Rwのchecked処理、byte/work/input-bit/intermediate-bit/slotの等号と1不足を検査 |
| 所有・失敗 | snapshotはadmission後、GCD検査前にcopy。7確保全fault、累積課金、最終allocator取消、再入、再利用、元値保持、診断、moveと元入力/Workspace/Plan破棄後の寿命を確認。曲線資格への昇格APIなし |
| 独立参照 | 下記の凍結Fraction参照と4×4 exact Hermite連立解。零・符号・非等間隔・clip等号・極近時刻・大原点・正時間尺度・負値尺度・3成分を含む |
| 実装検証 | D/R、変更3TU静的検査、7回以外のC++ newなし、1/2/4callerのcanonical語・work/slot/payload一致、独立変異、工程記録、独立レビューを完了 |

借用は同期call中不変とする。外側Workspaceは排他lease、call内の全scratchと所有結果は最初の確保前に合算する。scratchは終了時に解放し、結果は独立所有する。payloadはC++ object/arrayの計画値で、allocator metadata・stack frame・RSSの上界ではない。

## 参照・構成・変異

凍結v2 resultのSHAは `497287324f6084d162a41aa98bfad9993adb91287fb78e6e5467adef9a2e22c5`。新producerは元12例の係数を製品式によらないHermite4条件の連立解で96件照合し、exactな変換と3成分化で60例を生成した。fixture SHAは `c8b8afb2b25acd691eef0060b8a38afe6dedf56eb9e6bbe8e10cd877c903c1e8`。別directoryの再生成と独立レビュー担当の別作業場所でbytesが一致し、上書き・原本hash不一致を拒否した。producerはRuffと明示includeのstrict Pyrightを通過した。

Windows11/MSVC19.51.36260.0/C++23/sse2で、最終Debug/Releaseともnormal6case・2145assert、resource4case・112assert成功。seed176。変更3TUのclang-tidy21.1.1ではunused-but-set-variable/sign-compare/conversion/range-loop警告も有効にした。計測訂正後の製品TUも再検査して成功した。Linux/GCC14・sanitizerはCIへ渡し、成功を先書きしない。数値loopは通常の拒否・取消もstatus/Resultで返す。C++ new監査をC heap全体の動的監査とは呼ばない。

独立担当が数値source `b16c908`、cpp SHA `843e4d5fa5b2d6f6dd8bf48aa66fa3d7a6718d92f4cd2bd80f55e46de07f83fb` の私有複製をMSVC Debugでcompile/linkした。fixture caseのcontrolは2078assert成功。三次項のa−bをb−aにすると342件、端点の符号不一致条件を反転すると330件のassertion失敗を検出した。compile/linkは全て成功、Catch exitは0/42/42。生出力は独立担当の実行時ツール記録にあり、ここには結果を転記した。`run-mutations.ps1`は同じ固定版と対応objectsを要求する再現用であり、このスクリプト自身を再実行した結果とは区別する。

その後のba3e815は、既に予約していたforecast構造走査10visitをactual workへ保持した計測訂正である。数値式・係数・slot・予約値は不変。差分と最終D/Rで確認し、独立変異とfixture再生成は重複実行していない。独立担当は型/数値式/分岐合流を通読し、方式の変更や利用者判断が必要な穴は見つからなかった。

## 小入力の工程と再現

`measure.cpp`を固定Release libraryへ/O2・/fp:strictでリンク。Intel Core i5-13400F、1caller、3節点1成分、1000回。count合計4,764,700ns、全prepare＋破棄15,639,200ns。後者は再count・7確保・数値・copyを含み、差を純数値時間とは扱わない。計画7928B、128slot、17,216,718stepに対し実123slot・3035step。B75/w3。workは採用primitiveのvisit上界でありCPU命令数ではない。RSSは4,345,856→4,354,048→4,354,048B、最後の値がprocess生涯peak。工程内瞬間peak・実務規模性能は認定しない。build_infoはba3e815、uncommitted=falseである。

buildは既存C3専用 `out/issue166-build` のNinja Multi-Config、manifest install OFF、既存issue99 vcpkg prefixを利用。対象は `uvlas_test_aero_pchip_coefficients` と `_resources`、各exeへ `--rng-seed 176 --abort --reporter compact`。旧slotの既存演算bodyは変更していないため、成功済みの先行検証を再利用し無関係なsuiteは回していない。新規stdoutは既存redactorで公開用copyにし、前後hashとLF化をpublication.jsonへ記録した。既存採用仕様・証拠は書き換えていない。
