# Linux 配布候補の実行時条件

AlmaLinux 8.10 の glibc 2.28 上で、CPython 3.12 用の候補 wheel と CLI が同じ外部参照付き Case を受理した。この候補を再現するときは、ABI に加えて、入力文書の確保モデルに使うヘッダーと共有ライブラリー、外部参照を開くカーネル機能を照合する。正式な対応環境の選定は REL-1 に残る。 [移設時の Case と入力指紋](im06-private-wheel-external-references.json) [REL-1](https://github.com/YUKIKEDA/UVLAS/issues/26)

## glibc 2.28 の候補

固定元 `3a140d7bfbffef797a3fb7b29187dc458e50f53a` を、PyPA `manylinux_2_28_x86_64` の manifest `sha256:55a36963bc0504cdf004011001c1072306402f0b8b83be5a5eba7068f93eec53` で構築した。旧 Ubuntu 版が glibc 2.35 で要求シンボルの不足により起動できなかったのに対し、このビルドは CLI の起動と拡張の import まで通った。実 Case は未認定の libstdc++ ヘッダー版 `20250110` のため、入力文書の準備で拒否された。 [PyPA の配布画像](https://github.com/pypa/manylinux)

隔離した候補 `a13bce0b1622b420a4b0a0083a0b8f3d9d70de66` で、配布元ソースと確保の照合結果に基づいて入力文書の認定表を追補した。確保試験4件・387アサーションと入力文書試験10件・636アサーションに合格した。`auditwheel 6.8.2` が `manylinux_2_27_x86_64` と `manylinux_2_28_x86_64` の両タグを付けた wheel を空の仮想環境に導入した。CLI、Python ファイル入口、辞書入口が BDF・MSH・両対応表を受理した。CLI とファイル入口は、BDF の改変と正しい SHA-256 を持つ別 Geometry の Gmsh 対応表を拒否し、復元後に再び受理した。所有権試験は10件通った。 [隔離候補の差分](linux-release-runtime-qualification/candidate.diff.txt) [候補の実行結果](linux-release-runtime-qualification/candidate-result.json)

CLI と拡張が要求するシンボルの最大版は `GLIBC_2.26` と `GLIBCXX_3.4.21` だった。実行は glibc 2.28 と WSL2 の組で確認した。 [wheel の ABI 検査](linux-release-runtime-qualification/auditwheel.log) [CLI と拡張の要求シンボル](linux-release-runtime-qualification/qualified-abi.log)

## ヘッダーと例外コンストラクターの来歴

コンパイル時の `__GLIBCXX__` が同じでも、使用されたヘッダーと、実行時に `runtime_error(char const*)@GLIBCXX_3.4.21` を提供したライブラリーは変わり得る。同じ診断バイナリーを AlmaLinux 8.10、同じ環境で AlmaLinux 8.9 の署名済み libstdc++ を指定した状態、Ubuntu 22.04.5 で実行した。実際の例外コンストラクターは指紋が異なる3本にバインドしたが、候補が報告したソース指紋は全て AlmaLinux 8.10 のものだった。8.9 と 8.10 の対照では配布候補の CLI と導入済み拡張にも同じ違いを確認した。 [実バインドと指紋](linux-release-runtime-qualification/provider-proof.log) [CLI と拡張の対照](linux-release-runtime-qualification/alternate-provider.log)

同じ日付の GCC 14.2.1 ソース RPM に含まれる配布パッチ適用前の `vector.tcc` を使う対照でも、候補は実際に使ったファイルと異なる指紋を報告した。配布ビルドの台帳には、対象翻訳単位が実際に使う5ヘッダーを `-H` または依存ファイルから特定し、その SHA-256 を記録する。配布対象の環境では、構文エラーを通した実プロセスの例外バインド先とその指紋を取得し、対応ソースと確保回帰を照合する。今回の対照で判明したのは来歴表示の誤りであり、確保見積りの不足は観測していない。 [ヘッダー対照](linux-release-runtime-qualification/provider-proof.log)

## Linux の外部参照取得

Ubuntu Base 24.04.5 に移設した固定元 `3a140d7b` の CLI と wheel を、WSL2 カーネル `6.18.58.1-microsoft-standard-WSL2` 上のネットワークを遮断した環境で実行した。通常は `openat2` を6回使い、Case と4参照ファイルを受理した。最初の呼出しに `ENOSYS` を返すと CLI は Case 取得を拒否して終了コード2になった。3回目に `ENOSYS` または `EINVAL` を返すと、ルートと Case を開いた後、CLI は `/model/geometry/1/mesh` で終了コード2になった。同じ3回目の `ENOSYS` を受けた Python の `load_case` は `FailureCode.Unsupported` を返した。 [実行ログ](linux-release-runtime-qualification/openat2.log)

配布対象では、`openat2` と使用するパス解決フラグが実際に利用できることを確認する。候補 wheel のタグや glibc の検査ではこの条件を判定できない。この試行は機能がない場合の外部参照取得を確認した。古いカーネルでの実行、保存先ファイルシステムとチェックポイント公開・再開の認定は残る。 [openat2 の仕様](https://man7.org/linux/man-pages/man2/openat2.2.html)

## 再実行

[先行の移設手順](private-wheel-external-references.md#再実行)で作った CLI、wheel、Case を使う。Linux で `strace -f -qq -e trace=openat2` を通常の隔離先起動コマンドの前に付ける。生のトレースで Case が2回目、最初の外部参照が3回目であることを確かめ、`-e inject=openat2:error=ENOSYS:when=3` と `-e inject=openat2:error=EINVAL:when=3` を個別に追加する。Case 取得の対照は `when=1` を使う。Python は導入済みの `_uvlas` から `load_case` を呼び、返却値の `failure.code` を読む。

AlmaLinux の候補は固定元の別 checkout に上記の差分を適用してコミットする。移設 probe の `COMMIT` に再現時のコミットを設定し、同じ画像の GCC Toolset 14、CPython 3.12.16 と `-march=x86-64 -mtune=generic` を使って構築する。最初の wheel を `auditwheel repair --plat manylinux_2_28_x86_64` で処理する。導入時の probe は生成された複数タグ付きの wheel 名と `lib/python3.12/LICENSE.txt` を使い、開発ツールを `PATH` から外して実行する。

実行時来歴は、不正な JSON を読ませたプロセスを `LD_DEBUG=bindings` で観測し、例外コンストラクターを受け持ったパスの指紋を照合する。 [AlmaLinux 8.9 の対照 RPM](https://vault.almalinux.org/8.9/BaseOS/x86_64/os/Packages/libstdc++-8.5.0-20.el8.alma.x86_64.rpm)
