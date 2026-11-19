# Windowsでの相対openと親の共有条件

Windows 11 x64（build 26300、ローカルNTFS）、MSVC 19.51。採否判定に使う [最小対照のソース](reference-snapshot-windows-probe.cpp) を2通りにビルドした。通常の親ディレクトリとジャンクションをそれぞれ同じ根から相対openし、返った属性で区別できることを確認した。通常ファイルとそのハードリンクは同じvolumeとfile IDを返した。

| 親ディレクトリの権限 | 共有条件 | ハンドル生存中のrename | 解放後 |
|---|---|---|---|
| 属性参照と同期 | 読取り・書込みの共有 | 成功（負対照、probe exit 6） | 対象外 |
| 上記と `FILE_LIST_DIRECTORY` | 読取り・書込みの共有 | 拒否（Windows error 32） | 成功（probe exit 0） |

再現には同じローカルNTFS上の新規ディレクトリを二つ使う。各根に `parent/original.txt` とそのハードリンク `parent/alias.txt` を作り、別ディレクトリを指す `junction` を置く。開発版MSVCでソースを `/W4 /WX /std:c++23preview` と `ntdll.lib` でビルドする。片方だけ `/DDO_LIST` を加える。それぞれのバイナリへ、独立した根、その `parent`、存在していない同じ根内の `parent-renamed` の絶対パスを渡す。対照は元の一時ディレクトリの名前を変更する。

この対照で確認したのは、親ハンドルの共有条件、ジャンクションの識別と実体IDである。別のWindows版、ネットワークドライブ、Linux、製品の読取り・変更検出・全資源上限は、開始契約に記載した後続の実機試験で扱う。
