# IM-06 Nastran CONM2の重心と慣性

対象ソースは **1ce657590a8f17e0dcc3cc138eb52a496fda71c3**。[IO4-03](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran系の限定読込み) の認証済みデッキから、指定のCONM2一件を既に読んだGRIDに結び付ける。CIDが0または空欄ならX1–3は基本系の偏心、CIDが−1なら基本系の絶対重心座標として読む。UVLASでは絶対重心の三欄を必須とする。空欄の偏心成分と慣性成分は零として数値化し、元の空欄は結果に保持する。[MSC Nastran 2022.1 Quick Reference Guide, pp. 1558–1560](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf#page=1594)

## SI変換と基準点

質量はkg、GRID・重心・偏心はm、慣性はkg m²。GRIDの基本系位置には既存のCP包絡を使う。重心テンソル $J_C$ の対角は書かれたI11、I22、I33、非対角は積慣性I21、I31、I32の符号を一度反転した値である。慣性の基底はどちらのCIDでも基本Cartesian。元のカード、欄位置、符号付き十進値と元単位を保持する。

正質量を必須とする。重心慣性が全零なら点粒子として扱う。それ以外の重心慣性は、正定値と主慣性の三角不等式を厳密有理で検査する。後者には $\operatorname{tr}(J_C)I/2-J_C$ の主小行列式がすべて非負であることを使う。偏心 $d=C-G$ がある場合も、粒子の重心へ回転慣性を足す必要はない。[粒子・偏心剛体と慣性](../../specs/structure/connection-adopted-specification.md#5-粒子偏心剛体と慣性)

GRIDのまわりへ移した慣性は $J_G=J_C+m[(d\cdot d)I-dd^\mathsf{T}]$ で求める。平行軸移動は一度だけ適用する。重心慣性は厳密値として保存し、位置とGRIDまわりの慣性は有向丸めの閉区間として返す。参照宣言、指紋、要素IDで同じ元デッキのGRIDとCONM2を照合する。正のCIDの基底、Caseの物質点・質量台帳との対応、実際の構造自由度への組立ては後続の工程である。

## 独立対照

質量0.5 kgの四点を重心の両側 $\pm(0.1,0.2,0)$ m、$\pm(0,0.1,-0.3)$ mに置き、両符号の積慣性を持つ参照を作った。GRIDの偏心を $(0.3,-0.4,0.1)$ mとし、元デッキにはmm、gで入力した。同じ重心をCID 0の偏心とCID −1の絶対座標でそれぞれ指定した。両経路の質量、GRID位置、偏心、重心、重心慣性を独立値と照合した。

四点それぞれの速度から運動エネルギー $273/100$ Jを求めた。独立した並進・角加速度と仮想速度に対する質量作用の仕事率 $53/20$ Wも四点で直接足した。GRID基準の慣性から計算した包絡がそれぞれの有理値を含むことを確認した。零慣性の点粒子、正のCID、明示を欠く絶対重心、負質量、正定値または三角不等式を欠く慣性、異なる元デッキ、作業上限ちょうどと1不足も確認した。

## 実行結果

固定ソースの刻印と変更なしの状態を照合した。MSVC 19.51.36260 DebugとUbuntu 24.04.5のGCC 14.2.0 Debug（ASan・UBSan）で io.nastran_coordinate_cards、io.nastran_shell_cards、io.external_planar_pg の各3件が合格した。CONM2の絞り込み試験は両環境で4ケース・311アサーション。Linuxのclang-tidyとWindowsのclang-format 22.1.3を通過した。

Linuxの分離コピーで平行軸の加算を二重にすると、両CIDの運動エネルギーと仕事率の計4アサーションが失敗した。元の式に戻して、CONM2の4ケース・311アサーションが再び合格した。
