# Case の取得元からチェックポイントを読む

[FD-08 と IO4-06 の採用契約](../../specs/foundation/case-checkpoint-compatibility-contract.md)に従い、型付き Case と、取得時に開いた読取り根を `ReferenceInput` で保持する。
ファイル入力の起点は Case の親ディレクトリ、辞書入力の起点は明示された `reference_root` である。
字句の正規化と番号の検査を終えた後、指定された一つの確定ディレクトリを取得する。
同じ最終ハンドルからマニフェスト原本とチャンクを読み、原本の SHA-256 を照合してから JSON を解析する。
成功結果はマニフェスト原本、解析したマニフェスト、検証済みのチャンクを所有する。
ホスト指定のディレクトリ読取りは従来の入口から同じ容器検査を使う。

## 独立したマニフェストで照合する

[試験の固定バイト列](../../../tests/io/checkpoint_test_fixture.hpp)は製品の保存処理を呼ばずにファイルへ置いた。
符号なし64ビット値42と binary64 の1.0をそれぞれ8バイトで保存する。
マニフェストには語型、形状、バイト数、チャンクの SHA-256、七つの期待指紋を指定する。
下のコマンドは標準ライブラリーでチャンクとマニフェスト原本の SHA-256 を検算する。

```bash
python3 -B - <<'PY'
from hashlib import sha256
from pathlib import Path
import struct

source = Path("tests/io/checkpoint_test_fixture.hpp").read_text(encoding="utf-8")
manifest = source.split('constexpr std::string_view manifest = R"(', 1)[1].split(')";', 1)[0]
for field in (struct.pack("<Q", 42), struct.pack("<d", 1.0),
              manifest.encode(), manifest.encode() + b" "):
    digest = sha256(field).hexdigest()
    if digest not in source:
        raise ValueError(digest)
    print(digest)
PY
```

[製品試験](../../../tests/io/checkpoint_directory_test.cpp)は、ファイル入力とメモリー入力の双方で、作業ディレクトリを変えた後に同じ番号を読めることを確認する。
Linux では取得済みの根を改名して元のパスに別のディレクトリを置き、元の実体から読めることも確認する。
[Linux のネイティブ取得試験](../../../tests/io/linux_reference_root_test.cpp)は、選択済みの確定ディレクトリそのものを改名し、同名の別ディレクトリを置いた場合も、保持したハンドルから元のマニフェストを読む。
Windows の試験は、途中と終端にジャンクションを置いて拒否を確認し、[ネイティブ取得の試験](../../../tests/io/windows_reference_root_test.cpp)は終端ハンドルを保持した間の祖先と終端の削除共有を検査する。
Linux は同じ場所にシンボリックリンクを置く。

絶対パス、ドライブ相対、UNC、ストリーム指定、根外への移動、末尾区切り、非正準の番号は宣言の不備として拒否する。
途中または最後がファイルである場合とリンクを拒否する。
マニフェストへ空白1バイトを加え、宣言 SHA-256 の不一致が JSON の項目数制限とディレクトリ列挙制限より先に報告されることを確かめる。
指定番号からチャンクを消した場合、同じ親の別番号を選ばず、以前の成功結果から元のバイトと値42を引き続き取得できる。
パスの入力・正規化後のバイト数、成分数、ネイティブ作業域、列挙の上限と取消も検査する。

## 再実行

Linux Release で関係する試験をビルドし、実行する。

```bash
cmake --build --preset linux-release --target uvlas_test_io_checkpoint_directory uvlas_test_io_checkpoint_manifest uvlas_test_io_checkpoint_chunk uvlas_test_io_checkpoint_publication uvlas_test_io_linux_reference_root uvlas_test_io_reference_input uvlas_test_io_reference_path
ctest --preset linux-release --output-on-failure -R '^io[.](checkpoint_(directory|manifest|chunk|publication)|linux_reference_root|reference_(input|path))$'
```

Windows Release では `pwsh tools/dev/msvc.ps1` を通して同等の対象をビルドし、`windows-release` で `linux_reference_root` の代わりに `windows_reference_root` を実行する。

## Ubuntu Release の結果

Ubuntu 24.04.5 LTS、カーネル `7.0.0-44-generic`、GCC 14.2.0 Release で、コードと試験の版 `462f27dd33dcc76b776297acee570e7e156c0adc` をビルドした。
対象の7 CTest が通過した。Catch2 の44ケース、5,055アサーションを確認した。
`io.checkpoint_directory` は5ケース、358アサーション。そのうち `[FD-08]` の2ケース、238アサーションで Case の相対取得を確認した。
`strace` で、Case の根から `runs/checkpoint-12` を取得した記述子が、`manifest.json` と両チャンクを開く起点になっていることを確認した。
マニフェスト SHA-256 の照合だけを除いたビルドでは、ファイル・メモリーの両経路で JSON または列挙の上限が先に報告され、回帰試験の8アサーションが失敗した。

Windows の試験コードは同じ Case の対照と、保持中の削除共有を確認するネイティブ対照を含む。今回の実行は Linux で行った。

## AlmaLinux 8.10 カーネルでの照合

[先行の実行環境と同じ署名済みカーネル](alma-8-kernel-runtime.md)を QEMU のネットワークなしゲストで起動した。
独立のシステムコール試験は、`openat2` に `O_PATH | O_DIRECTORY | O_CLOEXEC` と `RESOLVE_BENEATH | RESOLVE_NO_SYMLINKS` を指定した。
根を保持した後のディレクトリ置換、途中と終端のリンク・通常ファイル、根外と欠落の拒否を確認した。
独立試験の終了値は0だった。
製品試験は、Ubuntu GCC 14.2 Release の実行ファイルと同じユーザー空間ライブラリーをゲストの initramfs に載せ、同じカーネルで実施した。
`4.18.0-553.189.1.el8_10.x86_64` で `[FD-08]` の2ケース、238アサーションが通過し、製品試験の終了値は0だった。
この照合は、対応するネイティブ取得がこのカーネルで動くことを示す。
AlmaLinux の glibc 2.28 での配布候補の認定は、[配布物の実行時検証](linux-release-runtime-qualification.md)に従う。

| 受け取ったファイル | SHA-256 |
| --- | --- |
| 署名済み RPM から抽出したカーネル | `6c68c16ecb83e4884b10e3cb5b9f9814cdff9127264a1ec32e7edfc910565a69` |
| 独立システムコール試験の C ソース | `eb6f3d211055a96cc890bfd3479177de71d84ab0e267c4d4cda1c9e22b5023b9` |
| Ubuntu とゲストで使った製品試験の実行ファイル | `6e42ba314b83f2f0ee2b896f7e9e8ab55e3603d63a50815d6c881f0cbafe8d64` |
| 製品試験を含む initramfs | `934ac921057b0cadce7f12383e22b0bfa7a1784d51b5b18bdd42ef037814bc2e` |
| 製品試験のゲストシリアル出力原本 | `5a6eee6b4eb59c9a8b9298825eeef51fa5113234544ce02aaa2fabca268c703a` |

ゲストでは Ubuntu の試験ファイルを `/product-test` に置く。元の ELF インタープリターを同じパスに置き、`ldd` が示す依存ライブラリーを対応する soname で `/runtime` に置く。
静的 BusyBox の `/init` で procfs と devtmpfs をマウントし、書込み可能な `/work` へ移動して `LD_LIBRARY_PATH=/runtime /product-test '[FD-08]' --colour-mode none` を実行する。
その終了値を出力してからゲストを停止する。
QEMU 8.2.2 (Ubuntu Noble) で TCG を使用し、署名済みカーネルと上記の initramfs を次の条件で起動した。

```bash
qemu-system-x86_64 -accel tcg -machine pc -cpu max -m 768 -smp 1 \
  -nodefaults -display none -serial stdio -monitor none -nic none -no-reboot \
  -kernel "$vmlinuz" -initrd "$initramfs" -append 'console=ttyS0 panic=1 rdinit=/init quiet'
```

準備済み型からの七指紋の導出、復元する環境の認定、J1 と探索の実状態を組み直す処理は、対応する型と復元入口が揃った段階で検証する。
