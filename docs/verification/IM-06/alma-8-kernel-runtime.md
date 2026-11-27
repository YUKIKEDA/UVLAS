# AlmaLinux 8.10 カーネルでの外部参照取得

AlmaLinux 8.10 の `4.18.0-553.189.1.el8_10.x86_64` をネットワークなしの QEMU ゲストで起動した。
実カーネルはルートと相対パスの照合に用いる `openat2` を受理した。
glibc 2.28 用の候補 CLI と CPython 3.12 wheel は、Nastran BDF、Gmsh MSH、両形式の対応表を参照する同一 Case を受理した。
[公式 erratum](https://errata.almalinux.org/8/ALSA-2026-102218.html) [指紋と実測値](im06-alma-8-kernel-runtime.json)

## 同じ配布候補を照合する

GNU Make で実際に使用した四ヘッダーを認定候補と照合した固定ビルド `a13bce0b1622b420a4b0a0083a0b8f3d9d70de66` の移設物を使った。
旧カーネルのゲストと WSL2 の対照で、CLI・拡張・libc・libstdc++・Case と四つの参照の SHA-256 が一致した。
両環境で CLI の終了は0、Python のファイル入口と辞書入口は `ok=true` だった。
[四ヘッダー、CLI、拡張の指紋](im06-linux-parser-input-qualification.json) [旧カーネルのシリアル出力](alma-8-kernel-runtime/boot.log) [WSL2 の対照](alma-8-kernel-runtime/control.log)

ゲストで構文不正の JSON を読むと、CLI は終了2、Python のファイル入口は `FailureCode.InvalidInput` を返した。
この二つのプロセスのロード記録で、CLI と拡張の `std::runtime_error(char const*)@GLIBCXX_3.4.21` が `/lib64/libstdc++.so.6` に結び付いた。
リンク先は `/usr/lib64/libstdc++.so.6.0.25` で、先行の確保モデル監査で照合したライブラリーと同じ SHA-256 だった。

製品検査とは別に、ルートの `RESOLVE_NO_SYMLINKS` と子の `RESOLVE_NO_SYMLINKS | RESOLVE_BENEATH` を実 syscall に渡す二つの対照を両環境で実行した。
ルートのファイル記述子を取得し、その記述子から `wing.bdf` を開いた。
WSL2 の製品トレースではルート、Case、`wing.msh`、`gmsh.json`、`wing.bdf`、`nastran.json` の順に同じ制約を適用した。
[通常実行のトレース](alma-8-kernel-runtime/trace-cli.log)

同じ CLI を WSL2 で再実行し、ルートの `openat2` に `ENOSYS` を返すと、終了2と `manifest acquisition` を得た。
最初の参照 `wing.msh` に `ENOSYS` または `EINVAL` を返すと、終了2と `reference acquisition at /model/geometry/1/mesh` を得た。
Python のファイル入口と辞書入口は `wing.msh` の `ENOSYS` をどちらも `FailureCode.Unsupported` に変換した。
[失敗注入のトレース](alma-8-kernel-runtime/injected-openat2.log)

## 再実行

[先行の移設試験](linux-release-runtime-qualification.md#再実行)で用意した CLI と導入済み wheel を使う。
AlmaLinux 8.10 の `kernel-core-4.18.0-553.189.1.el8_10.x86_64.rpm` を公式 BaseOS から取得する。
固定した PyPA イメージに登録されている AlmaLinux 鍵で `rpmkeys --checksig --verbose` を通し、`rpm2cpio` で `lib/modules/4.18.0-553.189.1.el8_10.x86_64/vmlinuz` を取り出す。
[配布元](https://repo.almalinux.org/almalinux/8/BaseOS/x86_64/os/Packages/kernel-core-4.18.0-553.189.1.el8_10.x86_64.rpm) [署名とカーネルの指紋](alma-8-kernel-runtime/kernel-signature.log)

一時ディレクトリにイメージの `usr/lib64`、CPython 3.12.16、wheel を導入済みの仮想環境、移設版 CLI、Case と参照を置く。
`lib64` は `usr/lib64` へのリンクとし、仮想環境とインタープリターは元と同じ絶対パスに置く。
[ゲストの init](alma-8-kernel-runtime/init.sh) は静的 BusyBox で `proc` と `devtmpfs` をマウントし、`uvlas --build-info`、`uvlas validate`、`_uvlas.load_case` と `_uvlas.case_from_dict` を順に実行する。
[ゲストでの実行手順](alma-8-kernel-runtime/run-probe.sh) [Python と syscall の対照](alma-8-kernel-runtime/probe.py)
同じディレクトリを `cpio -H newc` と gzip で固めたものをゲストの initramfs に使い、元ディレクトリを WSL2 の `unshare -m -n chroot` でも実行する。
Ubuntu 22.04 の QEMU 6.2.0 と依存は、署名済みの apt 索引が示す SHA-256 を照合して隔離領域へ展開した。

```sh
qemu-system-x86_64 -accel tcg -machine pc -cpu max -m 2048 -smp 2 \
  -nodefaults -display none -serial stdio -monitor none -nic none -no-reboot \
  -kernel "$vmlinuz" -initrd "$initramfs" \
  -append 'console=ttyS0 panic=1 rdinit=/init quiet'
```

この結果は、上記の署名済みカーネルと固定実行物の組で `openat2` を使えることを示す。
同じ 4.18 系でも、別の配布元・更新版・ファイルシステムの能力はそれぞれ確認が要る。
ゲストのファイルは initramfs に置いた。
一般の保存先、チェックポイントの確定と再開、正式な対応環境の選定は各ゲートに残る。
