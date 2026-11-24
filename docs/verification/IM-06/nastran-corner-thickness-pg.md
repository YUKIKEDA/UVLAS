# IM-06 Nastran の角別厚さを Case の全厚と照合する

認証済みの CTRIA3 と CQUAD4 を、同じ Case が所有する平面 PG に対応付ける。Nastran の膜厚を Case の物理全厚として比較するのは、PSHELL の MID1・MID2・MID3 がいずれも同じ MAT1 を明示し、Case が対応する等方材料の一層を持つ場合に限る。各 MID の空欄はこの照合では受理しない。曲げ倍率は空欄または厳密に1、横せん断厚さ比と結合材料は空欄を要求する。既存の材料照合と面の位置・法線・被覆の認定も通過する必要がある。[FD-03 対応表の契約](../../specs/foundation/external-mesh-sidecar-contract.md)

## 元欄を保って SI に換算する

TFLAG が空欄または0なら、Ti を入力単位の長さとして SI に換算する。空欄の Ti には PSHELL の T を使う。TFLAG が1なら、明示 Ti に PSHELL の T を掛け、空欄を倍率1として扱う。角別指定に必要な T がないときは元面と欄を示して拒否する。成功結果は元の節点順、SI の厳密値、値を指定した欄を所有する。明示倍率には Ti と T 双方の欄を残す。[MSC Nastran 2022.1 Quick Reference Guide、CQUAD4 pp.1609–1611、CTRIA3 pp.1667–1669、PSHELL pp.2901–2902](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

## 面積が正となる範囲で重なる各 Case 厚さ面と比べる

Case の全厚は三角形ごとにアフィンである。元面の認定済みチャート写像に代入すると、元三角形では一次、元四辺形では双一次になる。外部厚さも同じ形状関数で表す。面積が正の共通部分で一致するには、各 Case 厚さ面の式を元節点に代入した値と、元面の全係数が一致しなければならない。ここで使う Case の式は比較対象の三角形の係数を表すもので、領域外の物理厚さを割り当てるものではない。線や点で接するだけの領域は照合しない。四辺形の元の双一次形状関数を使い、対角線は加えない。

非平行四辺形では、Case の $H=1+s+t$ m と元角の1、2、5/2、2 m が全域で一致した。元の節点順を反転しても値の取り違えはなかった。負対照は元角が1、1、1、2 mの単位正方形で、Case の対角線分割も四隅で同じ値を持つ。点 $(3/4,1/4)$ では元面の双一次厚さは17/16 m、Case は1 mになる。照合は差がある元節点・欄と Case の厚さ三角形を返した。

実デッキからの試験は、CTRIA3 の絶対厚さと CQUAD4 の厚さ割合を、Case 入力から作った所有 PG と比較した。両形式で欄を空けた角を含めた。倍率の変更を Ti と T の位置付きで検出し、T がない倍率指定と、膜・曲げで異なる材料番号を持つ PSHELL を拒否した。

この経路は均質な PSHELL の平面を対象とする。PCOMP の角別指定や、曲面 PG に対する変換時の厚さ場は別の照合を要する。

## 検証に使った構成

対象ソースは `182e67b3fac8943fe7f42b56e1e261e6204f8c86`。Windows 11 の MSVC 19.51.36260 Debug と Release、Ubuntu 24.04 の GCC 14.2.0 Debug（ASan・UBSan）で、`io.external_pg_regional_thickness`、`io.external_planar_pg`、`io.nastran_laminate_cards`、`io.reference_input` が各4件合格した。各構成で、領域厚さの独立試験は111アサーション、実ファイルの対応表試験は927アサーションだった。

隔離コピーで、共通部分の頂点に含まれる元節点だけを比較すると、四隅が一致する負対照を誤って受理した。元の実装へ戻すと、領域厚さ試験の111アサーションが再び合格した。変更した翻訳単位の clang-tidy と clang-format 21.1.2 も通過した。
