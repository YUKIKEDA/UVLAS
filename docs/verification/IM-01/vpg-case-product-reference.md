# 凍結CaseからPGの物理参照を準備
<!-- verification: A5-02 FD-03 GEO-01 SH-D02 V-PG -->

[R4・R5の凍結Case](vpg-case-reference.md)を元のハッシュに照らして再生成し、JSONバイト列から型付きCaseとSI値を準備した。34 Caseに含まれる36 Shellを物理IDで選び、所有するPGグラフ、係数辞書、積分領域、回転・原点付き参照幾何へ渡した。生成物のマニフェストのSHA-256は `a3464f54edfc40b1e2752f13bafa648c29a898faef877de7593b67db725e1fb4`。

各Shellで採用格子と層IDを照合し、三角形の面積を合計した。矩形はチャート面積1、R5の切断領域は4/5となった。内部点で物理全厚を、一定厚さまたは凍結入力の $H=h(1+s/10+t/20)$ と有理数で照合した。[R4の物理幾何](pg-r4-domain-geometry.md)と[R5の切断領域](pg-r5-domain-geometry.md)は、層面・横境界を含む位置、法線、全厚、物理測度を独立式と照合する。

追加の二層Caseはmmを指定し、非零原点と軸を巡回する回転を持つ。所有結果は元のSI Caseを破棄した後も使える。参照点の座標を $x=3.5+\zeta/40,\ y=-4+2s,\ z=5+2t$ m と照合した。入力予算の枯渇と存在しないShellを拒否し、成功済みの参照識別子を保持した。

凍結系列はReleaseの `coupling.frozen_case_pg_reference` で、単位と所有・失敗境界は各構成の `coupling.case_pg_support_motion` で実行する。局所参照から先では、系列全体の条件付け、必須場、全拘束の組立とIM-02への引継ぎを確認する。
