# PR #161の3指摘への訂正記録 v2

対象のレビューheadはaffe9c6。修正sourceは01bee09である。v1 producer・結果・報告・manifestは旧head時点の証拠としてbyte保持し、現行契約と新v2 producerをこの記録で再発行した。旧v1 manifestのsource hashを現在の契約へ付け替えていない。

| 指摘 | 修正と固定した期待値 |
|---|---|
| 4623229291 Periodic | period省略はInvalidInput。全ての正常例にperiodを明示した。節点(0,1,2,3)、値(0,1,-1,0)、period=3でt=4の値は1。省略とperiod=2は拒否 |
| 4623229298 D2 Automatic | 内部節点で左右exact一致時だけ許可する。節点(0,1,2)・値(0,1,2)の中央は左右D2=0でAutomatic成功。値(0,1,4)では左右D2=0/4でInvalidInput。Hold/Periodicの表も同じ条件へ統一 |
| 4623229303 入力世代 | 呼出側の準備hostがissuer識別値と入力世代を明示供給し、不変準備入力・curve/jet・保存値へ保持する。評価I/Fに独立expected_preparation_inputを明記。cache/restoreはstampと元入力実内容を別々に照合し、旧保存値へ現在stampを付け直さない |

既存state-contractとcore::Model::createは、内容fingerprintが同じでも準備時の世代を別に照合する。PCHIPもこれに従い、実在しないissuer serviceを内部に仮定しない。issuer名前空間と世代非再利用を管理するのは準備hostである。同じ不変ownerを共有する経路と、別所有copyをexact比較する経路を区別した。後者の比較結果を再利用するなら両owner/stampを保持し、そのowner対だけへ適用する。各queryで全配列走査することを一律に要求しない。

v2の設計witnessは、元の165 exact等値・16,129小整数GCD対に加え、Periodic省略を含む11拒否と10意味変異を実行した。D2固定正負例も実行し、実返却値をresult.jsonへ記録した。これはC++のcompile/link変異や実allocator試験ではない。

入力の同内容copy・同stampなら保存bytesがそのまま復元できた。同内容で世代だけ3→4、issuerだけ17→18、同stampで値だけを変更する各入力を、cache/restoreの両方でStaleGenerationと確認した（6拒否）。係数を壊すrestoreはInvalidInputの別1件である。issuer/generationは固定した準備host入力であり、hashや時刻から作っていない。この限定Python照合模型は、製品hostの発行/再開管理や実checkpoint容器を実装したものではない。

CPython3.12.8/Windows11。新producerはRuff check/format、明示includeのstrict Pyrightが成功。別出力で再生成し、SHA256 497287324f6084d162a41aa98bfad9993adb91287fb78e6e5467adef9a2e22c5のbyte一致を確認した。親担当も別directoryで独立再実行し同じ結果/hashを確認した。

12曲線のscalar_fixturesはv1と全値一致し、二次照合producerも不変である。このため変更のないSciPy11例121値の対照はv1の実行結果を明示再利用し、重複実行していない。再現時は`python witness.py --output <new-path>`を使う。出力上書きは拒否する。日本語とhash/traceの対象検査を行い、C++/CMake/依存/旧R&Dは変更していないため製品suiteを回さない。新CIの完了も待たない。
