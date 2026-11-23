# チェックポイントの確定保存を実機で検証する

[IO4-06](../../specs/foundation/io-execution-adoption-specification.md#io4-06-checkpoint) の保存手順を、Windows と Linux の実ファイルで確かめた。同じ番号に二つの書込みが競合しても一方だけが確定し、以前の確定先を置き換えなかった。

## 確定名を切り替える前に書込み結果を照合する

保存先の親ディレクトリを開き、そのハンドルを起点に、確定先と同じファイルシステムへ新規の一時ディレクトリを作る。chunk を flush して閉じた後、同じ実体を開き直し、長さと SHA-256 を照合する。全 chunk の照合後に `manifest.json` を書き、そのファイルも flush・close・再照合を通す。

最後に番号から作った `checkpoint-<十進数>` へ一度だけ rename する。Linux は `renameat2(RENAME_NOREPLACE)`、Windows は親ハンドルを指定する `NtSetInformationFile(FileRenameInformation)` を使う。確定名と一時名 `.uvlas-checkpoint-*` を別にした。確定前に失敗した場合は、この呼出しが作った一時ファイルだけを片付ける。プロセスが停止した場合は一時ディレクトリが残り得る。

## 競合と書込み中断を実ストレージで再現した

| 環境 | 観測 |
|---|---|
| Windows 11 x64、MSVC 19.51 Debug、NTFS | chunk と manifest のバイト列、確定済み・空の出力先の保持、競合時の一件成功、SHA-256 不整合時の片付け、取消・容量不足・不正名の拒否を確認。 |
| Ubuntu 24.04.5 x86_64、GCC 14.2 Debug、ASan・UBSan | Windows と同じ対照を実行した。別プロセスのファイル上限を12バイトに制限した。8バイトの最初の chunk が完了し、16バイトの二つ目で13バイト目を書こうとすると `SIGXFSZ` で停止することを確認。完全な最初の chunk と12バイトの二つ目が残っても確定名と manifest はなく、既存の確定内容を保った。同じ番号への再試行は成功した。 |

両環境で `io.checkpoint_chunk` と `io.checkpoint_publication` が合格し、保存試験の12回繰返しも通った。Windows の最初の対照では Win32 の rename 入口が相対名を拒否したため、親ハンドルを指定できるネイティブ入口に替えて同じ試験を再実行した。

この検証は操作者が管理する保存先でのファイル書込みとプロセス中断を対象とする。電源断後の保持、最新確定番号の探索、物理状態の再開は別の検証で扱う。[Linux の rename](https://man7.org/linux/man-pages/man2/renameat2.2.html) [Windows の rename 情報](https://learn.microsoft.com/windows-hardware/drivers/ddi/ntifs/ns-ntifs-_file_rename_information)
