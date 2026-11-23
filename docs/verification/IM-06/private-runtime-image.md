# 非公開のランタイム用ディレクトリを作成して実行する

IO4-02の先行検証として、ビルド済みのCLIとCPython拡張を、ライセンスと来歴を伴う一つのディレクトリへ配置した。WindowsとLinuxで、同じビルドから二度作った42ファイルがバイト単位で一致した。別の場所へ移してもCLIを起動でき、新規venvの通常のsite-packagesから拡張を読み込み、既存の所有テスト10件が通過した。

## 構成と使い方

Releaseを構成する際に `UVLAS_BUILD_PYTHON=ON`、`VCPKG_MANIFEST_FEATURES=python`、`UVLAS_ENABLE_PRIVATE_RUNTIME_IMAGE=ON` を指定し、`uvlas_cli` と `uvlas_python` をビルドする。WindowsのNinja Multi-Configでは `--config Release`、LinuxのNinjaでは構成時に `-DCMAKE_BUILD_TYPE=Release` を使う。任意のインストールコンポーネント `uvlas-private-runtime` は、ツールが明示して使う。通常の `cmake --install` ではこのコンポーネントのファイルを作らない。

CPython 3.12で次を実行する。ビルドしたソースの40桁のコミットと、実行環境用のNumPy 2.3.5のwheelを指定する。wheelは試験のvenvにローカルファイルから入れる。WindowsでCMakeがPATHにない場合は `--cmake` に実行ファイルを指定する。出力先の親ディレクトリを先に用意し、出力先自体は未作成にする。

```sh
python tools/qa/private_runtime_image.py stage \
  --build "$build" --vcpkg-share "$installed/$triplet/share" \
  --image "$image" --commit "$commit" --numpy-wheel "$wheel" \
  --ownership-test tests/python/test_ownership.py

python tools/qa/private_runtime_image.py verify \
  --image "$relocated_image" --commit "$commit" --numpy-wheel "$wheel" \
  --ownership-test tests/python/test_ownership.py
```

配置先には `bin/` のCLI、`python/` のネイティブ拡張、`share/uvlas/` のライセンスと通知、`build-info.json` と `inventory.json` が入る。インベントリは自身を除いた各ファイルの相対パスとSHA-256を保持する。検証時はコピー前後のインベントリを照合し、開発ツールのPATHを外してCLIの来歴を再確認する。拡張は新しい検証用venvのsite-packagesへ置く。通常の `import _uvlas` がそのファイルを選んだことを確認する。

## 失敗時にも既存の配置を保持した

- CLIの改変、推移依存の通知の欠落、別コミットの来歴、無効な拡張、コンポーネントの無効化を拒否した。無効な拡張はインベントリと照合できる状態でも、実際のimportで失敗した。
- 既存の出力先に置いたデータは保持された。並行した二件の配置では片方だけが成功し、残った42ファイルは単独配置の結果と同じだった。最終配置の直前に作られた空ディレクトリも保持した。Linuxでは `renameat2` の `RENAME_NOREPLACE`、Windowsでは既存先を置き換えないディレクトリrenameを使う。
- インストール用の一時ディレクトリと、検証で使うvenvを分けた。検証の一時データを片付けてから、完成ディレクトリを出力先へ一度で移す。失敗した試行は出力先や配置用の一時ディレクトリを残さなかった。
- Windowsの最初の並行試験ではOSが既存先を正しく拒否したが、試験側が日本語の診断を英語と照合して失敗した。競合時にツールが一定の診断を返すように直し、両OSで対照を再実行した。

## 検証環境

ソース `099ed73d9d14cf64aff35d2c7b544e82f6fa4263` を使った。両OSのCLIがReleaseと `uncommitted_changes=false` を報告した。vcpkgは先行の[ランタイム監査](runtime-artifact-audit.md)と同じ固定基準を使用した。

| 環境 | ビルド | 実行時Python | 二度の構成 | 移設後の所有テスト |
|---|---|---|---:|---:|
| Windows 11 x64 build 26300 | MSVC 19.51.36260.0、Ninja Multi-Config 1.13.2、CMake 3.30.4 | CPython 3.12.8、NumPy 2.3.5 | 42ファイル一致 | 10件通過 |
| Ubuntu 24.04.5 x86_64、WSL2 | GNU 14.2.0、Ninja 1.11.1、CMake 3.31.6 | CPython 3.12.3、NumPy 2.3.5 | 42ファイル一致 | 10件通過 |

[実行スクリプト](private-runtime-image-probe.py)で[Windows](private-runtime-image/windows.json)と[Linux](private-runtime-image/linux.json)の結果を採取した。今回のディレクトリは非公開の検証用で、C系ランタイムやPython本体は実行環境が供給する。配布物の形式、対応下限、署名、公開先、実際の解析と再開の受入はIM06-0とREL-1に残る。
