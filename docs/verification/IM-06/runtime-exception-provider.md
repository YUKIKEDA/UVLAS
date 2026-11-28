# 移設した Linux 成果物の例外ライブラリーを照合する

既存のランタイム監査は、生成元の記録と CLI・CPython 拡張の SHA-256 を照合してから両方を試す。Linux/glibc で `--linux-exception-provider-sha256` を指定した場合は、その移設先で同じ構文不正の JSON を読み、実際の例外ライブラリーも検査する。期待する SHA-256 は配布候補ごとの台帳から指定する。[使い方](../../build.md#非公開のランタイム用ディレクトリ)

構文不正の入力は `{"schema":x}` とした。CLI は終了コード2と byte 10 の構文診断、Python は `load_case` の `InvalidInput` を確認する。さらに glibc のロード記録で、対象の移設 CLI と拡張がそれぞれ `std::runtime_error(char const*)@GLIBCXX_3.4.21` に結び付けた提供元を特定する。各実行で提供元が一つに定まり、実ファイルの SHA-256 が期待値と一致したときだけ JSON を出力する。

## 実アーティファクトによる対照

| 環境 | 実行した Release のソース | 結果 |
| --- | --- | --- |
| AlmaLinux 8.10、glibc 2.28、CPython 3.12.16 | `a13bce0b1622b420a4b0a0083a0b8f3d9d70de66` | CLI と拡張の両方が `libstdc++.so.6.0.25` の同じバイト列を使用。付属 CPython の `stdlib/LICENSE.txt` も照合。 |
| Ubuntu 24.04.5、glibc 2.39、CPython 3.12.3 | `3a140d7bfbffef797a3fb7b29187dc458e50f53a` | 両方が `libstdc++.so.6.0.33` を使用。ディストリビューションの Python ライセンスを参照。 |
| Windows 11、CPython 3.12.8 | `63df6cc00d7180725370201dca3454ad35073cfd` | 従来の移設・所有権10試験に合格。Linux のオプション指定は拒否。 |

AlmaLinux の候補に別版の正規 `libstdc++.so.6.0.25` の SHA-256 を指定すると、実提供元の SHA-256 との違いを報告して終了コード1になった。

`3a140d7b` を同じ AlmaLinux で構築した成果物は、入力文書の確保モデルが未認定で、`UnsupportedMemoryModel` となり構文検査に進まない。この対照では診断の byte がなく、監査は構文経路未到達として終了コード1になった。両対照の標準出力に成功 JSON はなかった。

実測した監査ソースは `77dc9bd8b9e70cd21f179d9897610e1ec9b7ce9a`。Ruff の検査と書式、Pyright strict が通った。配布先の認定には、実際のビルドヘッダー、確保モデル、カーネル機能を引き続き組で確認する。[Linux 配布候補の認定境界](linux-release-runtime-qualification.md)
