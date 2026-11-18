# P1局所誘導核の独立参照

P1核の値と全13変数微分を、製品の境界積分式から独立した面積積分で照合するための参照である。19例について3基底×3成分の値と117微分を生成した。角度求積の細分化で観測した差は、固定受入許容差の最大約 `1.05e-5` 倍だった。近接点と正則化自己点の2例は、さらに256点・150桁まで精密化した。

この面積参照は、連続積分に対する厳密な誤差上界を与えない。製品が返す丸め誤差包絡は、同じ離散作用素の独立高精度評価で別に試験する。求積誤差と丸め誤差を混同しない。

## 面積積分を独立に評価する

生成器は [p1-area-reference-generator.py](p1-area-reference-generator.py)。製品コードや境界渦・sheetの式をimportしない。参照三角形を固定した重心座標の扇形に分け、半径方向の積分を解析的に処理し、角度方向だけを高精度Gauss求積する。摂動時も扇形の重心座標を固定し、幾何、法線、面積Jacobianを微分する。理想核の同一平面外点では、個々の扇形に特異点が生じる符号付き分割を避け、内部の重心から分割する。

半径変数を $$ \rho $$ とすると、分母は $$ D=A\rho^2+B\rho+C=A[(\rho+s)^2+t^2] $$、 $$ s=B/(2A) $$、 $$ t^2=C/A-s^2 $$ となる。変数 $$ u=\rho+s $$ に移し、次の原始関数を使う。

$$
\begin{aligned}
J_0&=\frac{u}{t^4r}-\frac{u^3}{3t^4r^3}, &J_1&=-\frac{1}{3r^3},\\
J_2&=\frac{u^3}{3t^2r^3}, &J_3&=-\frac{1}{r}+\frac{t^2}{3r^3},\\
J_4&=\operatorname{asinh}(u/\sqrt{t^2})-\frac{u}{r}-\frac{u^3}{3r^3},
\qquad r=\sqrt{u^2+t^2}.
\end{aligned}
$$

各 $$ J_j $$ は $$ u^j/(u^2+t^2)^{5/2} $$ の原始関数である。二項展開により次数0〜4の半径モーメントを構成する。生成時に40組のパラメータで原始関数の微分を独立に照合した。

入力の浮動小数点数はbinary64の正確な値としてmpmathへ渡す。面・辺・頂点queryは元の位置式から組み立てる。微分は13変数のcomplex-step、刻みは `1e-32` とする。

## 19例の範囲と観測した差

対象は理想核の近接点3距離、正則化した面・辺・頂点の自己点、skew三角形、遠方、理想核の同一平面外点・辺延長上、非平面四角形を構成する2子三角形である。各子は理想核と2種類の正則化半径を持つ。これらは局所核の参照であり、大域幾何、後流、連成、実務受入の認定を表さない。

角度求積を24、40、64、96点以降へ細分化し、前段との差を製品の固定受入許容差で割って比較した。19例の最大比は `1.043767975976733e-5`。90桁から110桁への変更による最大比は `1.489624075878338e-71` だった。値は65桁の文字列で保存する。元のR&D参照生成器は隔離コピーから2つの非特異な正対照だけに使い、面積積分とcomplex-stepの結果を対照した。旧R&Dのsourceと記録は変更していない。

[精密化生成器](p1-area-reference-refine.py)は、距離 `1e-5` の理想核近接点とskew三角形の正則化面自己点を再評価する。96→160→256点、125→150桁、complex-step `1e-32` と `1e-40` を比較した。角度細分化差の最大比は `3.7808509161009819799e-39`、精度変更差の最大比は `3.193805511280784e-107`、刻み変更差の最大比は `1.4625686398569393166e-53`。80桁の値を別fixtureに保存する。これらも観測した差であり、厳密上界へ読み替えない。

## 離散Gauss作用素の丸め包絡は別に調べる

[p1-gauss-reference-generator.py](p1-gauss-reference-generator.py)は標準ライブラリの整数とFractionだけでLegendre多項式、Sturm列、根の分離を構成する。8点と12点の全20根を310bitまで分離し、節点と重みを256bitの外向き区間へ量子化する。全40個の次数0〜 $$ 2n-1 $$ のモーメント包含、区間幅、重み削除の負対照を検査する。製品のNewton反復から独立した参照である。

遠方枝の丸め包絡は、製品と同じ12点×12点Gauss作用素を独立の高精度演算で評価して照合する。連続面積参照との差は、別の求積精度の試験で扱う。

## 別出力で再生成する

環境はPython 3.12.8、mpmath 1.3.0、NumPy 2.3.5。面積参照の2生成器は既存の出力を拒否する。Gauss生成器は上書きを拒否しないため、実行前に新しい空出力先であることを確認する。以下はリポジトリrootから実行し、`out/p1-reference-reproduction` は新しい空ディレクトリとする。

```powershell
$ref = 'docs/verification/IM-03'
uv run --with mpmath==1.3.0 --with numpy==2.3.5 --python 3.12.8 python -B -X utf8 "$ref/p1-area-reference-generator.py" --legacy-source "$ref/p1-area-reference-source/uvlm_induction_reference.py" --output out/p1-reference-reproduction/area.json
uv run --with mpmath==1.3.0 --with numpy==2.3.5 --python 3.12.8 python -B -X utf8 "$ref/p1-area-reference-refine.py" --generator "$ref/p1-area-reference-generator.py" --output out/p1-reference-reproduction/precision.json
python -B -X utf8 "$ref/p1-gauss-reference-generator.py" --output out/p1-reference-reproduction/gauss.json
```

各コマンドの終了コードを確認し、順に `tests/aero/fixtures/p1-area-reference-v1.json`、`p1-area-precision-v1.json`、`p1-gauss-legendre-v1.json` とSHA-256を照合する。今回、3種類とも別出力への再生成でbyte一致を確認した。sourceとfixtureのハッシュは製品検証manifestに収録する。最初の生成ログは [面積参照](p1-area-reference-generation.txt)と[精密化](p1-area-reference-precision.txt)に保存した。
