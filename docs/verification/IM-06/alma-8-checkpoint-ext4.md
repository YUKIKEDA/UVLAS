# AlmaLinux 8.10 の ext4 でチェックポイントを確定して読み戻す

[IO4-06](../../specs/foundation/io-execution-adoption-specification.md#io4-06-checkpoint) の確定保存とチャンク読戻しを、AlmaLinux 8.10 の署名済みカーネル `4.18.0-553.189.1.el8_10.x86_64` と ext4 の組で確かめた。
前回と同じカーネルをゲストで起動し、別の 128 MiB 仮想ブロック装置 `/dev/vda` を ext4、読み書き可で `/save` にマウントした。
[カーネルの署名・来歴](alma-8-kernel-runtime.md) [保存先と試験のログ](alma-8-checkpoint-ext4/boot.log)

## 途中書込みからの読戻し

ネイティブの `save_checkpoint_directory` を使って番号12を保存し、次の保存では書込みプロセスのファイルサイズ上限を12バイトにした。
プロセスが `SIGXFSZ` で終了すると、一時ディレクトリには1つ目のチャンク8バイトと2つ目の途中までの12バイトがあり、`manifest.json` と確定番号13は存在しなかった。
`read_checkpoint_directory` で同じ一時ディレクトリを指定すると `CheckpointCorrupt` を返し、`read_latest_checkpoint` は確定番号12と値42を返した。

同じ親と番号13で保存し直すと、最新の確定番号は13になり、チャンクを型付きで読み戻すと、値は 42、100、200 だった。
番号13への再保存は拒否された。Linux の確定処理には `renameat2(RENAME_NOREPLACE)` を使う。
ゲストが正常にアンマウントした後、ホストから読み取り専用で `e2fsck -fn` を実行して成功した。
さらに `debugfs` で両確定番号と一時ディレクトリの実バイトを読み、Python の `struct.pack` と `hashlib.sha256` で型、形状、値、マニフェストの指紋を照合した。
[保存・読戻しの試験コード](alma-8-checkpoint-ext4/checkpoint-recovery-probe.cpp) [実バイトの照合コード](alma-8-checkpoint-ext4/check-disk.py) [アンマウント後の照合](alma-8-checkpoint-ext4/disk.log)

## ビルドと再実行

GCC 14.2.1 Release の固定候補 `a13bce0b1622b420a4b0a0083a0b8f3d9d70de66` の製品ライブラリーに、既存の `io.checkpoint_directory` と `io.checkpoint_publication` の二つの試験ソース、および上記の連続試験をリンクした。
比較した `fb945dfd08582cc2e1e30f538b2a3f3476e827ad` のチェックポイント実装、Linux の参照ルート、既存試験は候補と同一だった。
ディレクトリ試験3件・120アサーションと保存試験6件・92アサーションが通った。
ディレクトリの作業先を `/save`、保存試験の `TMPDIR` を `/save/tmp` とし、両方を同じ ext4 内で実行した。
[コンパイル・リンクの記録](alma-8-checkpoint-ext4/build.log) [実行物と生ログの指紋](im06-alma-8-checkpoint-ext4.json)

前回の QEMU 起動手順で、小さい initramfs に三つの試験実行ファイルと同じ AlmaLinux の実行時ライブラリーを置く。
ext4 と virtio ブロック装置のモジュールは前回の署名済み `kernel-core` RPM から取り出し、解凍してゲストに置く。
e2fsprogs 1.46.5 で新しい 128 MiB のファイルを `mkfs.ext4` で初期化し、QEMU に `-drive file="$image",if=virtio,format=raw,cache=none` で渡す。
ゲストでは `/dev/vda` のマウントと ext4 の種別を確認してから試験を実行し、完了後にアンマウントする。
[ゲスト起動時の実行手順](alma-8-checkpoint-ext4/init.sh)

検査したのは、この ext4 上でのプロセス中断、確定ディレクトリの公開、指定番号と最新番号の型付き読戻しである。
J1 の履歴を含むソルバー状態の復元、連続計算との bitwise 比較、電源断や別の保存先、正式な対応環境の選定はそれぞれ別の受入項目となる。

## 再実行時の独立照合

再実行には `python3 check-disk-v2.py <ext4-image>` を使う。第1版の検証器とそのログは初回計測の記録として保存した。
第1版は、両確定版のマニフェストのチャンク一覧を空にすると、通常実行では失敗し、`python -O` では成功を出力した。
第2版は全条件を明示的に検査し、確定版とステージングの照合を終えてから結果を出力する。

通常、`python -O`、`PYTHONOPTIMIZE=1` の3通りで、空のチャンク一覧とステージングの第2チャンク欠損をそれぞれ与えた。全6通りが終了コード1で、成功行を出力しなかった。
同じ3通りで元の ext4 イメージの複製を読み、全4行が初回のログとバイト一致した。ログの SHA-256 は `82f8f96bd6a649872d7f82e7b8b24a5687eb7fa23e2f8c111099bf0fe682e74d`、第2版の SHA-256 は `ef29367312f52dd896ad50956f4d3de477e15dc4b0feb8d6adf2d60ca27783bc`。
[第2版の検証器](alma-8-checkpoint-ext4/check-disk-v2.py) [初回のログ](alma-8-checkpoint-ext4/disk.log)
