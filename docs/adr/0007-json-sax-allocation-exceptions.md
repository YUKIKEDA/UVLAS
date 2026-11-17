# 0007: JSON SAXの割当失敗を例外境界へ返す

- Status: Accepted
- Date: 2026-11-17
- Issue: #98

## Context

[入力文書契約](../specs/foundation/input-document-contract.md)は、構文解析の割当失敗を`ResourceExceeded`へ変換し、プロセスを終了させない。Issue #93の開始審査で、固定nlohmann/json 3.12.0をMSVC Debugで使う場合に、この境界へ例外が届かない経路を確認した。

MSVC 14.51のDebug用`vector`と`basic_string`は、空のdefault constructorやmove constructorの内部でもiterator proxyを確保する。それらの一部は`noexcept`である。独立した最小プログラムでは、正常時exit 0、default構築の`operator new`へ`bad_alloc`を注入した場合は外側catchへ届かずterminate handlerのexit 77、同じ空値の明示count構築ではcatchのexit 42となった。

parserの`states`だけを直しても不十分である。lexerのtoken領域と外側constructor、例外文字列の構築・返却も通る。名前付き文字列の返却は、NRVOが行われなければMSVCのmove constructorへ入る。NRVOを無効にした対照で、この未修正経路のterminateも確認した。

## Decision

上流版3.12.0を保ち、vcpkg overlay portでSAXが通る3headerの9箇所だけを補修する。

| 対象 | 補修 |
|---|---|
| `detail/input/parser.hpp` | 空のstatesをthrow可能なcount constructorで作り、例外説明の返却を直接prvalueのcopy構築にする |
| `detail/input/lexer.hpp` | token用vector/stringとerror文字列を同じ空値の明示constructorで作る。外側の過剰な`noexcept`を外し、error文字列をNRVOに依存せず返す |
| `detail/string_concat.hpp` | 結果を同じ空値の明示constructorで作り、直接prvalueのcopy構築で返す |

標準ライブラリやABIを変更せず、Debugのiterator検査を全体で無効にしない。独自JSON parser、手置きheader、依存の共有prefixへの上書きも採用しない。JSONの文法・イベント・値・失敗codeは変えず、割当失敗だけを既存guardまで伝播可能にする。

baselineのnlohmann-json 3.12.0 port #2の取得先・SHA512、既存のchar8_t/optional修正、設定、MIT通知を維持する。UVLASの補修を加えたoverlayのport revisionで識別する。導入は専用prefixへ行い、実際のportから消費した回帰試験をIssue #99で閉じる。入力文書本体の開始gateは、その導入・検証と外部資源forecastの審査である。

## Verification

私有の依存copyを使った限定試作では、正常の全kind、不正surrogate、数値overflow、深いstack、長tokenと末尾不正の5literalについて、原依存と補修候補の成否・callback数・通常のerror ID（101/406）が一致した。全イベント列の値を比較した試験ではない。最終候補はNRVO無効のMSVC Debugで、controlが通った後の全200個のC++new地点への注入が外側catchへ届き、terminateは0件だった。これは選定したSAX経路の証拠であり、全nlohmann機能の保証や正式port導入済みの主張ではない。

Issue #99では正式port経由で同じ性質を回帰試験にし、影響する固定版・build-info・計測JSONのconsumerも確認する。LinuxはWindowsの結果から成功と推定しない。通常値の一致は許容差0で確認し、数値方式・物理の受入精度は変更しない。

MSVC CRTが例外messageを`malloc`で複製する経路はC++newの試験と別である。文書入口は固定sourceに基づいてmessage保存失敗のfallbackを識別し、確保不要の診断へ翻訳する。その分類対照と、実CRT mallocへ故障を注入した実証を混同しない。

## Consequences

返却文字列のcopyが増えるため、parser資源forecastは原文字列とcopyの同時領域を含める。実行時間・確保量が原依存と同じとは主張しない。ここは入力準備の処理であり、物理反復kernelの無確保を緩める判断ではない。

将来の依存更新では、補修の各経路とMSVC/STLの例外宣言を再確認する。上流で解消した場合も、通常値と故障経路の回帰を確認してから補修を外す。過去の依存source hashやR&D・検証記録は、その版のまま保持する。
