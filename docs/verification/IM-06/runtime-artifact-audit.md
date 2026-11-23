# CLIとPython拡張の初回ランタイム監査

ソース `df11771ce33902c6944cee1671defeb54f844934` からWindows x64とLinux x86_64のReleaseを構成した。CLIとCPython拡張をそれぞれビルド先から一時ディレクトリへ複製し、SHA-256を照合した。開発ツールの検索パスとPythonのユーザー環境を取り除き、CLIの来歴、拡張の読込み、既存の所有テストを確認した。成果物のハッシュ、各依存、通知、観測したABIは[Windows](runtime-artifact-audit/windows.json)と[Linux](runtime-artifact-audit/linux.json)に記録した。

| 項目 | Windows | Linux |
|---|---|---|
| 実行環境 | Windows 11 build 26300、AMD64 | Ubuntu 24.04.5 x86_64、WSL2カーネル |
| ビルド | MSVC 19.51.36260.0、CMake 3.30.4、Ninja Multi-Config 1.13.2 | GNU 14.2.0、CMake 3.31.6、Ninja 1.11.1 |
| vcpkg triplet | `x64-windows` | `x64-linux` |
| 実行時Python | CPython 3.12.8、NumPy 2.3.5 | CPython 3.12.3、NumPy 2.3.5 |
| 拡張ABI | `cp312-win_amd64` | `cpython-312-x86_64-linux-gnu` |
| 隔離起動 | CLI来歴0、所有テスト10件、欠落拡張1、未知オプション2 | CLI来歴0、所有テスト10件、欠落拡張1、未知オプション2 |

## 実行時の依存と通知

- Windowsの両ファイルはAMD64のPE32+。CLIはMSVCP140、VCRUNTIME140、VCRUNTIME140_1、UCRTのAPI set、KERNEL32を参照する。拡張は同系統のランタイムに加えて `python312.dll` を参照する。CLIはSystem32だけをPATHに入れて起動した。Python内で実際に読み込まれたモジュールを照合すると、MSVCPはSystem32、VCRUNTIMEと `python312.dll` はCPython本体の配下にあった。同じ隔離環境の補助ローダー検査ではUCRTのAPI setがSystem32の `ucrtbase.dll` へ解決した。これはPython内の観測であり、短命なCLIプロセスのDLL一覧は採取していない。
- Linuxの両ファイルはx86_64のELF64。両ファイルのRPATHとRUNPATHは空で、Nobleのローダーは全依存を解決した。要求シンボルの最大値は、CLIが GLIBCXX_3.4.32、拡張が GLIBCXX_3.4.31 と CXXABI_1.3.15。glibcについては両方とも GLIBC_2.34 だった。同じCLIをUbuntu 22.04.5へ複製すると、libstdc++の GLIBCXX_3.4.31 と GLIBCXX_3.4.32 が見つからず終了コード1になった。
- 両OSで36件ずつ、vcpkg通知を対応するインストール済みポートの原本と照合した。37ポートのうち、vcpkg内部のアンインストール用ポートは通知ファイルを生成しない。対象は直接のEigen、Spectra、Boost、nlohmann-json、pybind11を含み、ヘッダーのみの依存も確認した。CPythonの通知、NumPy wheelの通知と同梱ライブラリ、Linuxのシステムパッケージの通知は機械可読記録に所在とハッシュを載せた。Windowsのランタイム配布条件はベンダーの資料を配布設計で照合する。 [Visual Studio 2026の再配布情報](https://learn.microsoft.com/en-us/visualstudio/releases/2026/redistribution)

## 再実行

両OSで同じソースとvcpkgの基準 `2ec0f1747d8186a5c9babe276e2ad57cd2b12d2a` を使う。CMakeに `-DBUILD_TESTING=OFF -DUVLAS_BUILD_PYTHON=ON -DUVLAS_WARNINGS_AS_ERRORS=ON -DVCPKG_MANIFEST_FEATURES=python` と、ビルド用CPython 3.12の絶対パス、vcpkg toolchainを渡す。Linuxでは `CC=gcc-14 CXX=g++-14` と `-DCMAKE_BUILD_TYPE=Release` を指定した。Windowsは `tools/dev/msvc.ps1` を通してNinja Multi-Configで構成し、ビルドに `--config Release` を渡した。ビルド対象は `uvlas_cli uvlas_python` の2件。

配布条件と独立した新規venvに `tests/python/requirements.txt` のNumPyを入れる。監査ツールを、そのvenvのPythonまたは同じ版のPythonで実行する。Windowsでは、PEを検査するため VSLANG=1033 を設定し、監査ツール自体を `tools/dev/msvc.ps1` 経由で実行する。起動対象のプロセスは監査ツールが作る環境を使う。

```text
python tools/qa/runtime_artifact_audit.py --cli <uvlas> --extension <_uvlas-module>
  --python <absolute-fresh-venv-python> --ownership-test tests/python/test_ownership.py
  --notices <build>/third-party-notices --vcpkg-share <vcpkg-installed>/<triplet>/share
  --commit df11771ce33902c6944cee1671defeb54f844934
```

LinuxはUbuntu Base 24.04.5の公式amd64アーカイブをSHA-256で照合し、WSL上に別のユーザー空間を作った。GCC 14とCPython 3.12はUbuntuの署名付きパッケージから入れた。ホストの標準GCCは11.4.0であり、この検証のビルドに使っていない。 [Ubuntu Base 24.04](https://cdimage.ubuntu.com/ubuntu-base/releases/24.04/release/)

## 失敗した試行と今回の範囲

Windowsでは旧キャッシュを使った試行の来歴にJSON補修0が出た。現行のマニフェストを使う別のビルドへ切り替え、補修1を照合した。旧バイナリはPythonの最適化モードでも監査が終了コード1で拒否した。通知のコピーから推移依存のboost-coreを外した対照も終了コード1だった。LinuxではHTTPの取得が滞ったためHTTPSに切り替えた。ミラーから補ったパッケージは、署名済みのUbuntu索引にあるSHA-256と照合した。初回はGNU 14.2.0で両ファイルをリンクできたが、WSLの再起動で一時領域が消えた。記録には、保持できる領域でビルドから隔離起動まで再実行した値を使った。

今回確かめたのは、記録した2環境で現行CLIと所有型のPython拡張を起動できること。ソルバーの実行・再開、配布物の導入、manylinuxタグや最小OS、署名と配布先はIM06-0とREL-1に残る。
