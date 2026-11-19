# Reject/Holdのexact domain分類

固定製品source `fa5213e3a939e3fa0681eae2aa1d02d49f204aad` は、採用済みAE-07のReject/Holdについて、exact時刻の位置と内側・端点定数候補を返すprivate部品である。Hold端点では両候補を残す。Periodic、Automaticの導関数一致、C1、最終side選択、世代や公開曲線の資格はIM03-0の後続であり、この結果をそれらの資格には使わない。

## 完了条件と検証

| 条件 | 対応と結果 |
|---|---|
| exact分類・拒否 | 既存60sourceの両端、全節点、両端の前後2^-400を独立cpp_rational比較で判定。Reject外側はInvalidInput、未知policyはUnsupported、非canonical入力とinactive sourceも拒否。許容差0 |
| 既存入口への接続 | 内側候補を既存locatorからIntervalRequestへ渡し、評価のexact時刻と端点値を確認。定数候補は元CoefficientBlockの端点値と全成分一致 |
| 資源・所有・失敗 | 非確保countの後に3確保。全確保故障、累積課金、byte/work/slot/bit境界、最後の確保取消、再入拒否と再利用、元結果の保持を専用2ケースで確認 |
| 対象Debug・静的・変異 | normal 2ケース18,783 assertions、資源2ケース84 assertions、Windows clang-tidy21.1.1の製品/normal/資源3TUが成功。隔離端点誤分類はcompile/link後にCatch42、15assertions中1失敗 |

独立レビュー担当が資源試験・静的検査と変異を実行した。変異では最初の端点をBeforeへ誤分類し、正常時に受理されるRejectを拒否した。無変異controlは通常Debug試験であり、別の隔離controlを実行したとは主張しない。変異の元にした製品sourceと通常試験の製品sourceのSHA-256は一致する。転記した観測の範囲はindependent-observations.jsonに記録した。

## 計画と境界

queryと両端の3入力、比較2回だけを数える。事前の形状検査と計画構築は5訪問、実行のimport・比較・完了は6訪問で、既存slot演算費用に11訪問を加えて全callを予約する。全節点は再走査しない。所有済み係数の厳密増加を前提に使い、元入力を変えない。返却は借用を持たないscalarであり、別sourceのidentityやcacheとして使えない。

旧採用仕様・fixture・証拠は変更していない。対象外の旧suite、追加計測、ローカルReleaseとLinuxは実行していない。実務規模の性能、全曲線やUVLM接続の完了は主張しない。今回の変更は新しい局所部品とAeroの局所CMake登録だけであり、対象Debugと静的検査に限定した。
