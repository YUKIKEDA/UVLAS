# CLI と Python wheel を空の仮想環境で試す

コミット `63df6cc00d7180725370201dca3454ad35073cfd` の Windows / Linux Release を、ビルドディレクトリから別の配置へ移した。
ネイティブ拡張をローカル専用の wheel に格納し、空の仮想環境に `pip --no-index` で導入した。
両OSで `pip check`、通常ケースの `validate`、Python の所有権試験10件が通った。
外部 Nastran デッキの SHA-256 照合と、その後の必須 sidecar の未対応診断も確認した。
生成元の記録、wheel、インストール後の実行物とライブラリの SHA-256 は[検証結果](private-wheel-rehearsal-results.json)にまとめた。

## 導入先と依存ライブラリ

| 実行環境 | wheel タグと実行時の依存 |
| --- | --- |
| Ubuntu 24.04.5 Base の再展開から作った Linux x86-64 chroot。CPython 3.12.3 と必要な共有ライブラリを既存の同系統環境から配置し、開発ツールやソースを入れずに実行した。カーネルは WSL 2 と共有する。 | `cp312-cp312-linux_x86_64`。CLI と拡張が要求する最も新しいバージョン付きシンボルは `GLIBC_2.38` と `GLIBCXX_3.4.32`。実行時は glibc 2.39 と `libstdc++.so.6.0.33` を使用した。 |
| Windows 11 Home 10.0.26300、CPython 3.12.8 x64。新しい仮想環境を作り、`PATH` を `System32` に限定した。ホストには開発ツールがあり、独立した Windows 環境での試験は残る。 | `cp312-cp312-win_amd64`。CLI と拡張は MSVC CRT と UCRT を要求する。拡張は CPython から `python312.dll` と `vcruntime140*.dll`、Windows から `msvcp140.dll` と `ucrtbase.dll` を読み込んだ。 |

どちらも NumPy 2.3.5 の既存のバイナリ wheel を使用した。
UVLAS の wheel に収めた38個のライセンス・告知ファイルは、それぞれ同じビルドから配置した CLI 側のファイルと一致した。
CPython の告知ファイルは、各実行環境で読み取って SHA-256 を記録した。
Linux のタグは検証した環境向けで、他の Linux ディストリビューションや manylinux との互換性は未検証。

## 外部ファイルの照合とケースの所有

公開スキーマ例の `coupled_transient` に、三次元 GRID 1点の Nastran デッキと sidecar の参照を加えた。
両ファイルの SHA-256 が合う場合、CLI は `/model/geometry/1/sidecar` を示して終了コード2を返し、Python は `Unsupported` を返した。
デッキの末尾にコメントを足して宣言済みの SHA-256 を変えずに再検査すると、CLI は `/model/geometry/1/deck` の SHA-256 不一致で終了し、Python は `InvalidInput` を返した。
外部参照付きケースの検証成功まで進むには、sidecar の役割の実装が必要となる。

外部参照のない `structural_static` では、CLI の終了コード0と Python のケース取得を確認した。
Python では結果オブジェクトを破棄しても、取得したケースから元の SHA-256 とバイト数を読めた。
Linux ではケースの所有中に元ファイルを削除しても値が一致した。
Windows では所有中の削除が共有違反（`WinError 32`）で拒まれ、ケースを解放してから削除できた。
既存の [Python 所有権試験](../../../tests/python/test_ownership.py)も、インストール先の拡張で実行した。

## 再実行

1. 検証対象のコミットを clean な作業ツリーに置く。既存の [非公開 Release 配置](private-runtime-image-producer-review.md)に従い、`uvlas_runtime_artifact_manifest` をビルドする。`cmake --install <build> --config Release --component uvlas-private-runtime --prefix <stage>` で配置する。
2. ビルド側のネイティブ CPython 3.12 で `python -I -B private-wheel-rehearsal-probe.py build <stage> <producer-json> <numpy-whl> <new-transport>` を実行する。生成元の記録と両実行物の SHA-256 を照合してから、CLI の配置と wheel を作る。wheel は [Wheel 仕様](https://packaging.python.org/en/latest/specifications/binary-distribution-format/)の `METADATA`、`WHEEL`、`RECORD` を持つ。
3. 移設先に `<new-transport>` と[試験スクリプト](private-wheel-rehearsal-probe.py)を移す。[`case-schema-examples.json`](../../specs/foundation/case-schema-examples.json) と `test_ownership.py` を `<fixtures>` に置く。移設先の CPython 3.12 で `python -I -B private-wheel-rehearsal-probe.py probe <transport> <fixtures> <new-workspace>` を実行する。空の仮想環境を作り、wheel からの導入、CLI、外部参照、所有権を確認する。

試したのは R-06 の非公開導入手順である。求解・連成実行や、対応環境全体での配布受入には別の検証が必要となる。
