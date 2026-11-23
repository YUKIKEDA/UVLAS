# Release 成果物の来歴と CPython 通知の再検証（R-06）

Windows と Linux の実ビルドで、CLI と Python 拡張を生成元の記録に照合し、別ディレクトリでの起動と Python 所有権試験 10 件を各環境で確認した。
同名で内容の違う拡張は移設前に拒否し、CPython の通知が欠けた場合は探したパスを報告して終了した。

## 同じ Release ビルドから照合記録を生成する

Python を有効にしたネイティブビルドで `cmake --build <build-dir> --config Release --target uvlas_runtime_artifact_manifest` を実行する。
CLI の `--build-info` と、同じビルドが生成した CLI・拡張の SHA-256 を `<build-dir>/runtime-artifacts-Release.json` に記録する。
[前回の監査手順](runtime-artifact-audit.md)の引数に `--producer-manifest <build-dir>/runtime-artifacts-Release.json` を追加する。
監査ツールは指定されたコミット、成果物名、SHA-256 を移設前に照合し、移設後のバイト列と CLI のビルド情報も確認する。

## 取り違えと通知の欠落を両環境で再現した

| 環境 | 拡張の差し替え | CPython 通知の欠落 |
|---|---|---|
| Windows 11 AMD64、MSVC 19.51.36260.0、CPython 3.12.8 | 別の実ビルドが生成した同名の `.pyd` を指定し、終了値 1。 | 通知だけを除いた CPython の作業用コピーから仮想環境を作成。ベースの `LICENSE.txt` のパスを報告し、終了値 1。 |
| Ubuntu 24.04.5 x86_64、GCC 14.2、CPython 3.12.3 | `objcopy` で印を付けた同名の `.so` を作成。単独でインポートできることを先に確かめ、監査は終了値 1。 | 別のマウント名前空間で通知のディレクトリを空に見せた。`/usr/share/doc/python3.12/copyright` を報告し、終了値 1。 |

正常系の JSON に、生成元の照合記録の SHA-256 と、実際に読み取った CPython 通知のパス・SHA-256 を保存した。
検証したコミットと証拠ファイルは[検証記録](im06-runtime-artifact-producer-review.json)に記載した。
この記録は R-06 の初回実ロード監査の再検証に限る。配布用ランタイムの選定とリリースの認定は [REL-1](../../roadmap.md) で扱う。
