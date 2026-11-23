# 非公開ランタイムに Release 成果物の来歴を添える

非公開のランタイム用ディレクトリを作る前に、CLI と Python 拡張を同じ Release ビルドの記録に照合する。
この記録も配置先に保存するため、別の場所へ移した後でも元のビルドディレクトリなしで照合できる。

## 構成・配置・移設

[最初の検証で使った構成と実行手順](private-runtime-image.md#構成と使い方)に従い、CLI と拡張を個別にビルドする工程で、次のターゲットを実行する。
Windows は Ninja Multi-Config の Release、Linux は Ninja の `CMAKE_BUILD_TYPE=Release` を使う。

```sh
cmake --build "$build" --config Release --target uvlas_runtime_artifact_manifest
```

`stage` がビルド先の `runtime-artifacts-Release.json` を読み、記録と一致した二つの成果物を配置する。
配置したファイルも照合し、ビルドの記録をインベントリに加える。
`verify` はこの記録と、移設先で起動した CLI のビルド情報を照合する。

## 両環境で再配置と差し替えを確認した

Windows 11 AMD64 と Ubuntu 24.04.5 x86_64 の実ビルドで二度配置し、インベントリを含む43ファイルが一致した。
別ディレクトリへコピーして `verify` を実行し、両環境で CLI の起動と通常の site-packages からのインポートを確認した。Python 所有権試験10件も両環境で通過した。

同名で内容の違う拡張に置き換えると、`stage` は CMake を呼び出す前に終了値1を返し、配置先も一時ディレクトリも残さなかった。
配置後に拡張を取り替えてインベントリのハッシュを更新した場合も、`verify` は生成元の記録との差を検出した。
対照には、Windows では別の実ビルドの `.pyd`、Linux ではインポートできることを確認済みの、印を追加した `.so` を使った。

実行結果と確認したソースは[検証記録](im06-private-runtime-image-producer-review.json)にまとまっている。CPython 本体は実行環境が供給する。
