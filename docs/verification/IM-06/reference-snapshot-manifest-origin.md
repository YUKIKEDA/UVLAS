# manifestの取得元による実体照合

4例とも同じ根、呼出側ID `1`、manifestの元bytes（16進 `616263`）、長さ `3`、SHA-256 `ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad` を使う。外部参照の宣言は1件で、その宣言digestも同じ値とする。実体照合だけを確認するため、人工の元bytesと宣言を用い、ファイル形式固有の解釈は行わない。

| 取得元 | 参照パス | 実体 | 期待結果 |
|---|---|---|---|
| `File` | `manifest` | 元ハンドルと同じ | `Cycle` |
| `File` | `alias` | 元ファイルのハードリンク | `Cycle` |
| `File` | `copy` | 同じbytesを別々に書いたファイル | 成功し、参照先の3 bytesを所有 |
| `Memory` | `alias` | 根にある通常ファイル。取得元にファイルIDはない | 成功し、参照先の3 bytesを所有 |

新規の一時ディレクトリを根として使う。`manifest` を作り、`alias` をそのハードリンクとして作る。`copy` には別途3 bytesを書き込む。開いた元ハンドルからbytesを読み、候補も一つずつ開く。`File` の期待結果は元ハンドルと候補ハンドルのIDの比較から得る。`Memory` の入力には同じbytesを渡す。成功例は元bytesとdigestも照合する。`File` のハードリンクを許してしまう対照と、IDの代わりにdigestを比較して `copy` を拒否する対照が、この表と食い違うことを確かめる。

レビュー修正時に Windows 11 のローカルNTFS上で4例を照合した。Pythonの `os.fstat` で開いたハンドルの `(st_dev, st_ino)` を比較し、`hashlib.sha256` で元bytesと候補bytesを照合した。製品の入出力とLinuxのadapterは実装時の試験で確認する。

製品実装では、`sub/manifest` の `copy` 宣言に対し、`sub/copy` に `616263`、上位の `copy` に `616264` を置く例も確認する。上表と同じ宣言digestで照合し、`sub/copy` の3 bytesを返す。上位を根に渡す負対照はdigest不一致になる。
