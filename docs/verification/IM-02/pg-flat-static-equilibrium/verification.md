# 平面シェルの非アフィン試行と積分モーメント

[静的 Case の契約](../../../specs/structure/pg-flat-static-equilibrium-contract.md)に用いる独立検証を記す。Python の `Fraction` で、固定厚と層内一定の StVK を評価した。

## 厚さの必要次数

全厚2.54 mm、長さ2 m、Young率200 GPa、Poisson比0の矩形で、軸変位を $u_Y=a Ls\zeta^2$ とした。$a=127/10^6$、$s\in[0,1]$、$\zeta\in[-1,1]$。固定端では変位が零、$\det F=1+a\zeta^2$ は全域で正である。

独立に展開すると Green ひずみエネルギーに厚さ8次の項が生じた。面内を解析積分して厚さだけを比べたとき、3点 Simpson 則ではエネルギーの相対差が1.875873%だった。正の有理数則は0〜9次のモーメントに厳密に一致し、エネルギーも解析積分と一致した。 [厚さ積分の対照](trial_quadrature.py)

## エネルギー、元力、接線

別の対照では、回転・並進した矩形を幅の異なる二つの区画を各2三角形に分け、幅の異なる2層で積分した。Poisson比は2/7。4形状の元 Cartesian 12列を使い、三つの非アフィン試行で正の変形を確かめた。

三角形のモーメントで評価したエネルギー、元力12成分、接線作用12成分は、Cartesian 座標で Green ひずみの不変量を微分した別計算と有理値で一致した。エネルギーが四次となる方向に沿った5点の厳密な差分でも二つの変分が一致した。層の分割と接線の双対性を照合した。接線応力作用から幾何項 $HS$ を除くと、三試行の各12成分で参照から外れた。対照で保持した係数の整数幅は最大134ビットだった。 [作用の独立検算](moment_witness.py)

## 求積の構成と再現

二進有理節点と補間重みを使い、0次から指定次数までの全モーメントを厳密に照合した。面内次数36と Duffy 変換後の次数37には37点と38点の正の則を得た。厚さ次数8には9点の正の則を得た。面内の重みの整数幅は最大820ビットだった。二つの面内の則と厚さ則を別々に使って Duffy 則の6,327個のモーメントを求め、すべて連続積分と照合した。係数の次数包絡からエネルギー、元力、接線一方向に必要な積の上限を求めた。 [モーメントと作業量の検算](support_and_rule.py) [正の有理数則の検算](positive_rational_rules.py)

以下の4本を Python 3.12 の標準ライブラリで実行する。通常実行と `-O` の出力がそれぞれ一致した。作用の対照出力の SHA-256 は `26c23a9dc22d58014b4510cfd56e4b63e159bbabf5218c1ca1f922e6f2903994`、作業量の対照出力は `307eacb4fc6b1e4a6de0ea1bba5b6e81441c9b7a0a0282bd976b9bea7efe967b`。

```bash
python3 -B docs/verification/IM-02/pg-flat-static-equilibrium/trial_quadrature.py
python3 -B docs/verification/IM-02/pg-flat-static-equilibrium/moment_witness.py
python3 -B docs/verification/IM-02/pg-flat-static-equilibrium/positive_rational_rules.py
python3 -B docs/verification/IM-02/pg-flat-static-equilibrium/support_and_rule.py
```

製品では、参照 seed から求めた変位・元自由行の残差・拘束反力を三つの Case 入力経路で比較する。三つの経路の結果は、独立な軸伸長の式と FEM の不確かさを含む受入幅で照合する。有限則と連続積分の誤差を分け、反力は同じ受理状態で求める。固有対はその平衡を次の入力に使う。
