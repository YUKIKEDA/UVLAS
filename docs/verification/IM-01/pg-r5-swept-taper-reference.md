# R5 の後退・テーパ形状を製品の参照準備で照合する

[凍結したR5の第2版](../../../rd/structure/main-shell/paired-geometry/records/paired-shell-r5-domain-v2-result.json)の形状を製品の参照幾何・共通分割へ渡した。入力は、原試験が使う格子 $n=1,2$ と基準厚さ $h=0.01,0.0001$ の全4配置。結果JSONのSHA-256はR5監査が記録した `7552eb6b9cb4f9cebb07b31eeeb61abd3c392620ed9f5e897f1ba5ffe2ce7e43` と一致した。元実験、領域の関数、グラフ係数を展開する関数のハッシュも照合した。 [領域の原契約](../../../rd/structure/main-shell/paired-geometry/records/paired-shell-r5-domain-protocol.md)・[監査](../../../rd/structure/main-shell/paired-geometry/records/paired-shell-r5-audit.json)

共通分割が返す三角形の面積と重心から、面積 $4/5$、$s$ の一次モーメント $11/30$、$t$ の一次モーメント $131/300$ を厳密に検査した。全三角形で、有理の重心と正の法線深さを製品へ照会した。法線管の測度は製品のベクトル外積と独立に、グラフの第1・第2基本形式から平均曲率項とガウス曲率項を作って照合した。変厚の勾配は厚さ方向の接線と平行なため、測度には寄与しない。各三角形の面積と全厚に対応するJacobianを含む。

検証対象は、後退・テーパの共通分割、変厚、法線位置と体積測度。R5の保存場4件は索引・監査・結果の指紋を照合した。変位、応力、モードの製品比較では、保存場の実体と独立荷重を用いる。
