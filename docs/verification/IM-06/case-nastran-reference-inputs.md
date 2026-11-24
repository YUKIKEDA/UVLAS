# IM-06 認証済み Nastran の物理参照を一つにまとめる

認証済み対応表が同じ SI Case に属し、物性を読んだデッキが対応表の元デッキと一致することを確かめる。そのうえで元要素・Case 側 ID・対応表の宣言位置を保持する。[FD-03 対応表の契約](../../specs/foundation/external-mesh-sidecar-contract.md)

Shell には所有 PG と厚さ・積層・材料軸の照合結果を付ける。CBEAM/PBEAM の SI 端点・断面と基底、CONM2 の SI 質量・重心と慣性、RBE2 の参照点・方向は、対応表の対象 ID を添えた元デッキ側の物理量として保持する。対応を指定した元要素の順を保ち、全処理が成立したときにまとめて返す。

Case 梁の剛性や接続領域との物理的な等価性は、後段の準備で照合する。通常のファイル検査には構文・スキーマ・参照を確認する既存契約が適用される。[全 Shell の準備](case-nastran-planar-pg-set.md)

## 同じ実デッキから作った証拠と元診断

実デッキを認証し、二つの Shell のほかに、CBEAM 905 の長さ1 mと PBEAM の断面積2 m²、CONM2 906 の質量2 kg、RBE2 907 の従属点を照合した。梁・質量・接続はそれぞれの Case 名前空間の ID 1 に結び付けた。対応表とデッキの識別・ダイジェストを所有し、準備した梁がそのデッキを参照していることを確かめた。

RBE2 の ALPHA だけを非零にすると、既存の物理準備が返す理由と元デッキのバイト位置に、Case Joint ID と対応表の位置が添わった。その後も非零の ALPHA を含むデッキを保持し、同じ Case と作業領域で二つの Shell を再準備できた。梁・質量・接続の対応カード件数は三件ちょうどで成功し、一件不足は割当前に検出した。準備中に成功した割当は、後段で失敗した場合も共有 ResourcePlan に累積する。

## 検証に使った構成

対象ソースは `d8dc9510d5fc244ab36212381fa4e0d028e3978c`。Windows 11 の MSVC 19.51.36260 Debug と Release、Ubuntu 24.04 の GCC 14.2.0 Debug（ASan・UBSan）で、`io.nastran_shell_cards`、`io.nastran_pbeam_cards`、`io.nastran_laminate_cards`、`io.external_planar_pg`、`io.external_pg_regional_thickness`、`io.reference_input` が各6件合格した。各構成で、実ファイルの対応表試験は1696アサーションだった。変更した翻訳単位の clang-tidy と clang-format 21.1.2 も通過した。
