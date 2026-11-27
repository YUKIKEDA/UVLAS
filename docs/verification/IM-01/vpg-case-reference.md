# 凍結 V-PG から参照形状と固定支持を再生成する

[A5-02](../../requirements/practical-acceptance-plan.md#a5-02-部品部分系の必須入力と参照) の V-PG のうち、形状、厚さ、格子、材料定数、全厚固定支持と全厚接続を現行の `uvlas.case/1` に展開した。Case の構成は R4 の品質系列24件、明示係数系列4件、曲面接続系列2件、R5 後退テーパ系列4件である。元の保存順を保ち、各 Case を凍結 result の行番号まで追跡できる。

R4 のグラフは元のパラメータを十進値から有理数へ変換し、チャートの移動、幅、せん断を含めて展開した。任意係数系列は凍結 result の有理係数を直接使う。R5 のグラフは凍結ソースの構成関数からパラメータを読み、両側の直線で切った多角形を用いる。全厚は元の定数またはアフィン分布を使い、支持は `s=0` の辺に並進3成分の全厚拘束を置く。材料は凍結参照の単位密度と、同じ Lamé 定数を与える等方弾性値を使う。

曲面接続は凍結した二曲面の分割数を読み、左の `s=1` と右の `s=0` に、共有する `(t, ζ)` を渡す。両側の半セルの境界で接続範囲を分け、両側で対応するセルを指定する。外周辺は互いに逆向きに宣言し、全厚固定支持は左端に置く。

## 再生成と照合

```text
python -m tools.qa.vpg_case_reference generate out/vpg-reference
python -m tools.qa.vpg_case_reference verify out/vpg-reference
```

再生成先の親ディレクトリを先に作る。既存ディレクトリへの生成は拒否する。元 R&D は読み取り専用で、監査の固定ハッシュと、監査に記録された各 result・ソースのハッシュを確認する。出力の `manifest.json` は入力のハッシュと各 Case のハッシュを持つ。独立した出力先で再生成し、検証時には manifest だけでなく各 Case のバイト列も比べる。今回のマニフェスト SHA-256 は `a3464f54edfc40b1e2752f13bafa648c29a898faef877de7593b67db725e1fb4` だった。

Windows 11 の MSVC 19.51 と Ubuntu 24.04 の GCC 14.2 の Release ビルドで、34件すべてを `uvlas validate` と Python の `uvlas.load_case` で読み込んだ。両 OS で同じマニフェストを生成した。元 result の改変、result と監査のハッシュの同時付け替え、既存出力への再生成、生成後の Case 改変の4件が拒否された。

この再生の範囲は参照形状と支持・接続の Case 入力までである。曲面接続の物理資格は、凍結モデルを製品の型付きパッチから疎チャートまで通した別の試験で照合した。変位・応力・固有モードの比較では、凍結参照にある層別の状態と独立荷重を実装した Case を使う。再生結果の未載荷 Case とテンプレート由来の入力・出力設定は、その比較の合否判定に使わない。

形状と拘束の物理検証は [R4 の位置と Jacobian](pg-r4-reference-geometries.md)、[R5 の体積と仕事](pg-r5-swept-taper-reference.md)、[R4 の二曲面接続](pg-r4-curved-interface/verification.md)、[疎チャートの全厚比較](pg-chart-refinement-series/physical-trace-verification.md) に記録した。
