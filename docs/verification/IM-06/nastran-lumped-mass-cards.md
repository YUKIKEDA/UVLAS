# Nastran の CONM2 質量欄と参照系を読む

要素IDの重複を検査したデッキから CONM2 を読み、GRID と座標系の参照を元の欄と照合する。元デッキの宣言元と SHA-256、欄位置、正確な十進値を記録する。元のデッキは、ほかのカードを読む処理も引き続き使える。 [IO4-03](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran)

CID の空欄と 0 は基本系のオフセットで、正の CID は指定系の局所 Cartesian 成分を表す。CID が −1 のときは X1〜X3 が基本系における重心の絶対座標となる。慣性は重心まわりの I11、I21、I22、I31、I32、I33 を欄名に対応させて保持する。積慣性の符号と空欄を保持し、テンソルの非対角成分への符号反転、平行軸移動、単位換算、慣性行列の適格性の判定は後段で行う。 [MSC Nastran 2022.1 Quick Reference Guide、pp.1558–1560](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf#page=1594)

Windows の MSVC 19.51 Release と Linux の GCC 14.2 Release で、C++ の対象試験を各1件実行し、局所オフセットと絶対重心、慣性の六欄と元位置、予約欄・正の参照先の検査を確認した。確保または作業の上限が1不足しても元デッキを保持した。clang-tidy 21.1.1 は読取りと試験を検査した。GCC の分離コピーで I21 を I31 の欄から読むように変えると回帰が失敗し、復元後は通過した。
