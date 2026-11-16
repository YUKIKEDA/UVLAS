# 厳格な入力ID・index変換の実装契約

2026-11-17。Issue #85、IM06-1。[入力文書契約](input-document-contract.md)の整数・ID・失敗のうち、parserを必要としない変換だけを具体化する。[N4](numerical-foundation-adoption-specification.md)の強い型と値としての失敗を継承し、数値・物理範囲・JSONの意味を変更しない。

## 同期借用から所有する値へ変換する

製品入口は [input_scalar.hpp](../../../include/uvlas/io/input_scalar.hpp)。`read_physical_id<Id>`は復号済みの長さ付きASCII文字列、必須のcaller-owned optional診断slot、任意の既知FailureDetailsを受け、`core::Result<Id>`を返す。IdはDofId/PanelId/WakeVertexId/WakeFaceId/LayerId/ObservationId/ScopeId/ComponentIdだけに制約する。SourceID・世代・試行・確定番号・Indexへ暗黙変換しない。

入力は呼出し中だけ借用し、結果と診断は入力を参照しない。0は有効なIDで、欠落ではない。Model membershipや同数IDの一意性は別の検査であり、この変換はそれらを認定しない。

`read_index`は、int64/uint64のどちらかだけを明示構築できるNativeIntegerから、0〜INT64_MAXのcore::Indexを返す。bool/float/string/より狭い整数からの暗黙昇格を避けるため、constructorは型推論したexactなint64/uint64だけを許す。配列length内の確認と、寸法積/byte数の検査は呼出側の既存checked処理へ接続する。

SourceSnapshot/hash・JSON構文/tree・runtime node kindはここで生成も認定もしない。将来の文書adapterが元kindを確認し、同じraw sourceとSourceOriginを所有して、変換した値に対応づける。native整数0から元tokenが`0`だったと推定せず、`-0`等の字面は元文書が保持する。

## 固定文法、順序、失敗

IDは非空のASCII digit列だけを受け、先頭zeroを許す。符号、前後空白、小数点、指数、nonASCII、embedded NULを拒否する。 $$ U=2^{64}-1 $$ とし、左から一桁dを読むごとに $$ v\le(U-d)/10 $$ を確かめてから $$ v\leftarrow10v+d $$ とする。double、locale digit判定、符号を許す汎用文字列変換へ置換しない。

検査順はdiagnostic reset→empty→各byteのASCII判定→その桁でのoverflow判定→累積。最初の不正byteまたはoverflowで止める。MAX+1にXが続けばIdOverflow/19、MAXにXが続けばInvalidIdCharacter/20である。

indexはsigned/unsignedの区別を保持する。signed負値は符号反転/unsigned化する前に拒否し、INT64_MINのabsを作らない。unsignedはINT64_MAX以下と確認してからcastする。

| 条件 | code/補助reason/位置 |
|---|---|
| 空ID | InvalidInput / EmptyId / 位置なし |
| ASCII digit以外 | InvalidInput / InvalidIdCharacter / 最初のdecoded-byte ordinal |
| uint64範囲超過 | InvalidInput / IdOverflow / 最初に超える桁のordinal |
| native signed負値 | InvalidInput / NegativeIndex / 位置なし |
| native unsignedのIndex範囲超過 | InvalidInput / IndexOverflow / 位置なし |

operationはPrepareModel。core::FailureDetailsの既知欄は値として全て保持し、変換対象と別のdomainのproblem IDを作らない。補助診断の位置は復号ID文字列内で、元file/JSONのbyte offsetではない。補助診断は確保不要のenumとoptional位置を所有する。

正常な入口で診断slotをresetし、成功は空、失敗は現在の原因で置換する。過去のFailureと診断はcallerが所有copyして残せる。出力slotのalias/同時共有を許さず、独立slotでは共有するconst入力を並行に読める。入力、既知context、既存結果は変更しない。

本体は整数・byte・固定値だけを使い、追加確保をせず例外を外へ出さない。NativeIntegerのprimitive variantはvalueless状態にならないが、標準visitorに宣言された例外経路も既存guard_entryで封じる。予期しない標準library例外は既存のInvalidInput/bad_alloc時ResourceExceededを維持し、位置なしのUnexpectedException補助診断を付ける。通常の入力拒否をこのfallbackへ落とさない。core::Result、noexceptと既存Failureの意味を維持する。引数shape/表現は型で固定され、未知の数値kindをroundしてからこの入口へ通すadapterを作らない。

## 独立期待値と検証範囲

試験にはIO4-01/N4-02/N4-08、所有にN4-03、確保にN4-09の該当タグを付ける。比較は整数・byte・enumの厳密一致で、許容差は0。期待値を製品変換やJSON parserから生成しない。

固定系列は各physical domainの0/1/UINT64_MAX、2^53の前後、leading zeroと長いzero列、MAX+1/超過後の不正文字、全不正文字群、NULの明示長、signed INT64_MIN/−1/0/MAX、unsigned MAX/その次/UINT64_MAXである。型制約の正常対照と全IDdomain間の暗黙混用拒否、bool/float/stringからNativeIntegerへの構築拒否もコンパイル時に確認する。

診断slotの再利用/成功reset、FailureDetails全欄、入力破棄後の所有診断と旧結果、独立並列呼出しを別に確認する。Cpp allocation auditは対象区間だけを測り、Catch・fixture・thread生成を外へ置く。監査が明示的な確保を検出する対照も用意する。実装にC heap/Eigen/外部確保がないことを読み、宣言した監査範囲と照合する。

隔離変異は、double往復、符号許容、overflowまたはindex境界検査の欠落、leading-zero拒否等から独立要件を壊すものを選ぶ。control成功とcompile/link成功を先に確認し、実行assertの失敗だけを検出に数える。実行版・入力hash・元archive経路へ迂回していないことを記録する。

この部品の適合は、入力文書kernel、有効case、CLI/Python、mesh、prepare/run/resume、配布やIM-06全体の完了を意味しない。それらは親Issue #22が持つ。既存の数値solver・物理kernel・R&D参照の再実行は、この変換追加の検証へ含めない。
