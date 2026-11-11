# 0002: WindowsのPresetはNinjaとし、MSVCは環境から取る

- Status: Accepted
- Date: 2026-11-11
- Issue: #1

## Context

2026-10-29のビルド検証は、Visual Studio 2022 Build Tools（MSVC 14.44）と、Presetのジェネレータ「Visual Studio 17 2022」で行った。2026-11-11の時点で、開発機にあるのはVisual Studio 2026（MSVC 14.51）だけで、導入済みのCMake（3.30.4、4.1.2）はVisual Studio 2026のジェネレータを持たない。既存のPresetでは構成できなかった。

## Decision

WindowsのPresetのジェネレータを「Ninja Multi-Config」にし、コンパイラは環境の `cl` を使う。環境は [tools/dev/msvc.ps1](../../tools/dev/msvc.ps1) が最新のVisual StudioのC++ツールから読み込む。justfileのレシピはこれを通す。

## Consequences

- Visual Studioの版が変わっても、Presetを変えずに済む。LinuxとWindowsが同じNinjaに揃う。
- Windowsのコマンドは、MSVCの環境の中で実行する必要がある。素のシェルから `cmake --preset windows-msvc` を実行すると、コンパイラが見つからない。
- 使うMSVCの版は、導入されている最新のものになる。`cmake/CompilerPolicy.cmake` が19.44未満を拒否する。CIでは版を記録する。
- MSVC 14.44で得た2026-10-29の検証結果を、14.51での結果として扱わない。
