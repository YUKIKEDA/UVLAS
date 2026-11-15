# 0006: pybind11はheader-only port、CPythonは外部3.12を使う

- Status: Accepted
- Date: 2026-11-15
- Issue: #14

## Context

IO4-02はCPython 3.12、NumPy 2.x、pybind11 3.1.0を採用している。固定baselineのpybind11 portは `PYBIND11_TEST=OFF` と `PYBIND11_NOPYTHON=ON` でPythonを探索しないheader-only exportを作るが、manifestにpython3への無条件依存がある。この依存をそのまま導入すると、今回のbaselineでは別のCPython 3.14とそのruntime依存を構築する。製品が使う外部3.12とは別で、採用されたbindingの構築に必要ない。

## Decision

`cmake/vcpkg-ports/pybind11` にbaselineのportを写し、manifestから未使用のpython3依存だけを取り除く。pybind11の版3.1.0、GitHubの取得先・SHA512、portfile、usage、上流LICENSEの導入処理は維持する。製品を手動配置のヘッダやpip版pybind11へ切り替えない。

consumerのCMakeは外部 `Python3 3.12 EXACT` のInterpreter、Development.Module、NumPyをまとめて発見し、pybind11のFindPython経路で発見したinterpreterのModule targetに接続する。NumPyのmajorは2に制限する。Python featureには専用のvcpkg導入prefixを使い、通常のtestsやBenchmark構成を変更しない。

これは新しいruntimeの採用ではなく、採用済みheader-only依存と外部interpreterの役割をbuild graphへ反映する調整である。通常構成のPython機能はOFFのまま。wheelとruntime配布はIM-06で扱う。

## Verification

独立レビューで上流のNOPYTHON分岐とINTERFACE exportを読み、外部3.12でのconsumer compile/link/importも別途確認した。製品はoverlay経由の専用prefixを構築し、外部CPython 3.12.8とNumPy 2.3.5を使った実moduleのcompile/link/import、所有と失敗のCTestで確認する。Linuxは同じoverlayと外部3.12を使う専用CI jobで確認し、Windowsの結果から成功と推定しない。

## Consequences

上流portの更新時は、この差分とNOPYTHONの意味を再確認する。vcpkgが外部Pythonを選ぶことを期待せず、CMakeへ使用するinterpreterを明示する。元portのpython3依存を戻す場合は、採用runtimeと配布範囲を改めて判断する。
