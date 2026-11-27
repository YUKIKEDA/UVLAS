# 移設した CLI と Python で二形式の外部参照を照合する

同じ Release コミット `3a140d7bfbffef797a3fb7b29187dc458e50f53a` から CLI と Python 拡張を配置した。両環境で空の CPython 3.12 仮想環境に private wheel と NumPy 2.3.5 をオフライン導入し、`pip check` と Python 所有権の10試験が通った。移設先で読み込んだ実行物の SHA-256 を、生成元の記録と照合した。 [成果物と入力の指紋](im06-private-wheel-external-references.json)

| 実行環境 | 配置と導入の条件 |
| --- | --- |
| Windows 11 x64、MSVC 19.51.36260.0、CPython 3.12.8 | 同じホストの別ディレクトリへ移設し、開発コマンドを除いた `PATH` と空の仮想環境で実行した。タグは `cp312-cp312-win_amd64`。 |
| Ubuntu Base 24.04.5 x86-64、GCC 14.2、CPython 3.12.3 | 新しく展開した環境に実行時依存だけを配置し、ネットワークの経路がない名前空間で実行した。タグは `cp312-cp312-linux_x86_64`。要求シンボルの最大版は `GLIBC_2.38` と `GLIBCXX_3.4.32`。 |

## 同じ Case の成功と拒否

同じ回転平板と層を参照する実 BDF と実 MSH、その二つの対応表を一つの Case で宣言した。四辺形の元タグは両形式とも701。MSH の全厚250 mmと BDF の全厚0.25 mは、この Case の層厚に対応する。Case と4参照ファイルの SHA-256 は両OSで一致した。 [元の物理対照](external-shell-cross-format.md)

- 相互に対応する4ファイルを使うと、CLI `validate` は終了コード0、Python の `load_case` と `case_from_dict` は成功した。
- BDF のバイト列だけを変えると、CLI は `/model/geometry/2/deck` を示し、Python の `load_case` は `InvalidInput` を返した。
- 別 Geometry の Gmsh 対応表を作り、その内容に合わせて宣言済み SHA-256 も更新すると、CLI は `/model/geometry/1/sidecar` を示し、`load_case` は `InvalidInput` を返した。
- 元の4ファイルと Case を復元すると、CLI と `load_case` が再び受理した。

この記録の検証範囲は、移設した入口での外部参照の照合、拒否と復元、Case の所有である。一般配布ABI、独立した Windows 環境、物理求解と実状態からの再開はそれぞれの受入条件で判定する。

## 再実行

[先行の移設手順](private-wheel-rehearsal.md#再実行)を使う。台帳の SHA-256 を照合した `private-wheel-rehearsal-probe.py` を別ディレクトリの `probe.py` にコピーする。そこで[移設用の差分](private-wheel-external-references-probe.patch)を `git -c core.autocrlf=false apply --no-index --unidiff-zero <差分の絶対パス>` で適用する。続いて同じ `build` と `probe` を `probe.py` に対して実行する。差分は Release の固定コミットと、外部参照の対照を更新する。
