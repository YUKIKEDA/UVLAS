# Case から指定されたチェックポイントの照合

[FD-08](../../specs/foundation/case-schema-fields.md#fd-08-初期状態と再開要求) の `Restart.manifest_sha256` を、[IO4-06](../../specs/foundation/io-execution-adoption-specification.md#io4-06-checkpoint) の指定ディレクトリ読戻しにつないだ。
この記録の内部入口は、型付き Case、呼出し側で解決した宣言のディレクトリ、マニフェストに割り当てた SourceId、準備側から渡された七つの指紋を受け取った。
現在の Case 入力は、取得時の読取り根を使う[相対取得の入口](case-checkpoint-relative-root.md)へ渡す。
同じハンドルで取得した `manifest.json` の SHA-256 を宣言と照合してから、JSON を読んで指紋・チャンクを照合する。
一致した全チャンクを独立して所有する。

## 検証

既存の flutter Case に `restart` を指定し、実ファイルの番号12を読む。
元のマニフェストと、末尾に空白1バイトを加えた版の SHA-256 は、先行試験で作った ext4 イメージから `debugfs` で取り出したファイルを Python の `hashlib.sha256` で計算した。
空白を加えても同じ JSON として読めるため、従来の直接読戻しは成功する。
元の指紋を指定した Case は `CheckpointCorrupt` を返す。
新しい指紋を指定した Case は読戻しに成功し、列挙上限が0なら `ResourceExceeded` を返す。
七つの指紋のうち材料だけを変えると `CheckpointIncompatible` を返す。
拒否後も先に取得した値42を読める。

## 実行結果

AlmaLinux 8.10 の GCC 14.2.1 Release で、既存の固定ビルドに変更した翻訳単位をリンクした。
ディレクトリ試験4件・136アサーション、マニフェスト試験4件・112アサーション、変更のない保存試験6件・92アサーションが通った。
宣言 SHA の照合だけを除いたビルドでは、空白追加版を元の指紋で読むと、拒否理由が `CheckpointCorrupt` から `ResourceExceeded` に変わった。
追加の試験はこの差を検出した。
[Linux の実行記録](case-checkpoint-manifest-digest/linux.log) [変異の記録](case-checkpoint-manifest-digest/mutation.log)

MSVC 19.52.37029.0 Release は対象3本のビルドと、ネイティブ保存試験を通った。
この環境の STL は `_MSVC_STL_UPDATE=202609L` で、現行の入力パーサーが採用する `202604L` とは版が異なる。
直接入力試験で `UnsupportedMemoryModel` を確認し、この環境ではディレクトリ・マニフェスト読込みの試験が失敗した。
Windows の読込み成功は認定済みの STL で別に確認する。
[Windows の実行記録](case-checkpoint-manifest-digest/windows.log) [入力試験の診断](case-checkpoint-manifest-digest/windows-parser-probe.log) [指紋と結果](im06-case-checkpoint-manifest-digest.json)

## 接続に必要な条件

Case の `Restart.directory` は非空の UTF-8 文字列として採用済みである。
この時点では呼出し側が用意したパスを受け取っていた。
ファイル Case と辞書 Case の相対パスは[後続の契約](../../specs/foundation/case-checkpoint-compatibility-contract.md)に従う。
準備済みモデルからの七指紋の導出と、[J1 の受理済み状態と探索履歴](../../specs/coupled/coupled-state-search-contract.md)をチャンクへ割り当てる規則は、実行・再開を接続する側で確定する。
同一環境での bitwise 一致と異環境の観測予算は、物理状態を復元する側で検証する。
