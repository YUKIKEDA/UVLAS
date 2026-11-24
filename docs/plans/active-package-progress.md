# IM-01・IM-03・IM-06の進捗と残作業

実装中の3パッケージは **24／43段階（約56%）が完了**した。残り19段階のうち16段階が実装・検証、3段階が終了審査。下表は現在の成果と次の完了条件を示す。

| パッケージ | 完了／全体 | 残り | 統合済みの成果 | 次に完了させる成果と判定条件 |
|---|---:|---:|---|---|
| [IM-01 構造](#im-01-構造) | 11／12（92%） | 1 | 一つのCaseから、認定基底と有限求積点を保持する全パッチ、その全厚拘束と固定疎チャートを準備できる。凍結R4・R5の32配置で有限則の合同誤差を測り、必須場を厳密に照合した。 [Caseの準備仕様](../specs/structure/pg-case-chart-contract.md)・[32配置の照合](../verification/IM-01/pg-vpg-finite-conditioning/verification.md) | Caseで準備した全厚拘束を凍結接続の原行と照合し、同じパッチ・拘束・チャートをIM-02へ引き継ぐ。 |
| [IM-03 空力](#im-03-空力) | 7／13（54%） | 6 | 局所P1核、大域幾何、後縁所有、共有辺の誘導・全微分を実装した。非一様な既知流れを同じ試行の循環・荷重とその微分に使う。全指定運動と励振のイベントを共通時計で検査し、同じ時刻の片側で離散空力取付を評価する。型付き入力の再開時には現行ケースと保存した係数・時刻表を照合する。 | PGの物質点写像と循環時間微分を接続し、自由後流の時間更新とJ1の状態確定・再開を検証する。 |
| [IM-06 入出力](#im-06-入出力) | 6／18（33%） | 12 | 型付きschemaとSI変換、マニフェスト直下の全参照の取得、Gmshメッシュの読取り、チェックポイント容器を実装した。容器は既存番号を保護して保存し、検証を通った最新番号のチャンクを復元する。 [チェックポイントの検証](../verification/IM-06/checkpoint-manifest.md) | Gmshのsidecarで部材、材料、支持、接続、参照面への対応を定め、役割ごとの参照を解決する。実行再開にはJ1と探索の受理済み状態が必要。 |

割合は[作業分解](next-package-work-breakdown.md)の完了段階数を全段階数で割って求める。各段階の完了は、採用仕様、後続が受け取る成果、独立した検証根拠で確かめる。段階ごとに規模が異なるため、この数字から残りの工数は見積もれない。

## IM-01 構造

[#17](https://github.com/YUKIKEDA/UVLAS/issues/17) は固定疎チャートまで実装した。複数の指定支持を実際のケースから組み立て、R-02の全厚位置差を実空間で照合した。残る1段階はV-PGとIM-02への終了審査。

| 段階 | 状態 | 到達点と残作業 |
|---|---|---|
| 厳密代数 | 完了 | 有理数・疎多項式の準備。 [#91](https://github.com/YUKIKEDA/UVLAS/pull/91) |
| グラフとセル | 完了 | 有理グラフ、half-cell、矩形の活動map。 [#103](https://github.com/YUKIKEDA/UVLAS/pull/103) |
| 全原辞書と full V | 完了 | アンカー補完、必須場と並進の厳密な写像。 [#113](https://github.com/YUKIKEDA/UVLAS/pull/113) |
| 多角形と全積分分割 | 完了 | 被覆を検査した分割と資源計画。 [分割](https://github.com/YUKIKEDA/UVLAS/pull/128)・[資源計画](https://github.com/YUKIKEDA/UVLAS/pull/149) |
| 局所独立化 | 完了 | 局所rank、基底、写像と来歴を所有する。再利用時にグラフ・領域・活動面積・辞書を照合し、B/T/S/W/Zを全活動行で厳密に再検証する。全12系列を独立参照と照合した。 [準備の検証](../verification/IM-01/pg-rank-bounded-preparation/verification.md)・[再検証と資源](../verification/IM-01/pg-rank-bounded-preparation/validation-verification.md) |
| 参照幾何と微分 | 完了 | グラフ・辞書・共通分割を所有する入口から参照Jacobian、可変法線、実法線厚さ、高階・片側微分と物理測度を返す。全活動片の局所正則性を厳密に検査する。 [実装 `09e983f`](https://github.com/YUKIKEDA/UVLAS/commit/09e983ff5d63d9ab95e0b2c75883687b64f219b1)・[検証](../verification/IM-01/pg-reference-geometry/verification.md) |
| 保証付き条件付け | 完了 | 一定厚・区分一次厚さの平板を厳密に積分する。曲面は所有する正の積則で幾何内積を包絡し、合同変換と必須場を認定する。切断・変厚放物柱面は全130列・34,680点の製品変換を独立な連続内積へ適用し、最大行和誤差が $1.062\times10^{-23}$ 以下。15点則の連続誤差は基準の $10^{-20}$ を超えた。求積次数の適格性は対象形状と分割ごとに確認し、V-PGで実用系列を検証する。 [変換生成](../verification/IM-01/pg-conditioning-generated-transform.md)・[独立連続積分](../verification/IM-01/pg-curved-gram-continuum-v1/verification.md) |
| PreparedPatchと参照作用 | 完了 | 認定済み基底と物理参照を所有する。体積点、層界面の上下、外周側面で、参照位置、正の測度、形状値、Cartesian 勾配を取得し、成分別の誤差枠で Binary64／TwoWord に射影する。固定則の点・重みを保持し、同じ Workspace で体積位置と転置を作用させる。 [物理参照量の射影](../verification/IM-01/pg-prepared-reference-projection.md) |
| 全位置traceと拘束 | 完了 | 非整合面を共通細分し、厚さ段差、指定層と部分面、向きを検査する。原拘束と固定支持の指定値を保持し、厳密な重複行の尺度と双対仕事を照合した。時間依存の指定値と複数接続の閉路処理は固定ピボットの疎チャートで扱う。 [全厚拘束の検証](../verification/IM-01/pg-prepared-patch-trace.md) |
| 固定pivot疎chart | 完了 | 固定支持、物理拘束ID、元行順で代表を選ぶ。同じ因子で位置・接線・転置・元行反力を計算し、独立に作ったKKT解と照合した。Caseの複数Motionは各群の元行を個別認定し、共通時計で片側の位置・速度・加速度を評価する。R-02の2・3・4パッチでは元有理CSRによる補正後、全厚位置差を全域で $1.153\times10^{-16}$ m以下と確認した。 [複数Motion](../verification/IM-01/pg-structural-motion-groups.md)・[実空間の残差](../verification/IM-01/pg-chart-refinement-series/physical-trace-verification.md) |
| 早期R-02 | 完了 | 非整合を含む2・3・4パッチ、接続1・2・4本の系列で、原行の非零数、rank、L/U、準備・再利用の課金と常駐量を測った。4パッチの閉路はrank 248、因子23,503項。実項数ちょうどで準備に成功し、チャート識別子、rank、L/U各項数が一致した。1項不足では元拘束を保持して拒否した。実務規模の性能はA5-06で確認する。 [R-02の計測](../verification/IM-01/pg-sparse-chart-series.md) |
| 終了審査 | 残 | 凍結34 Caseの36 Shellから製品の物理参照を準備した。R4の28配置とR5の4配置で、位置、法線、全厚、層の片側と物理測度を独立式と照合した。32配置では指定有限則の合同誤差が $10^{-20}$ 以下となり、385,968成分で必須場の厳密一致を確認した。凍結した二曲面の整合・非整合接続には独立した全厚拘束の参照があり、Caseから同じ準備を通す検証に使う。IM-02へはパッチ、有限点と重み、原拘束およびチャートをまとめて渡す。 [Caseからの物理参照](../verification/IM-01/vpg-case-product-reference.md)・[R4の物理幾何](../verification/IM-01/pg-r4-domain-geometry.md)・[R5の切断境界と層](../verification/IM-01/pg-r5-domain-geometry.md)・[R4の二曲面接続](../verification/IM-01/pg-r4-curved-interface/verification.md)・[有限則と資源の照合](../verification/IM-01/pg-vpg-finite-conditioning/verification.md) |

## IM-03 空力

[#19](https://github.com/YUKIKEDA/UVLAS/issues/19) の終了審査まで6段階。共有辺・四辺形の誘導作用を循環計算へ渡し、既知流れと指定運動を共通時計で評価する。循環とその微分が、全荷重・後流へ進む共通の先行条件。

| 段階 | 状態 | 到達点と残作業 |
|---|---|---|
| 厳密な局所述語 | 完了 | 退化、内外、距離とfar判定。 [#90](https://github.com/YUKIKEDA/UVLAS/pull/90) |
| 所有幾何・ID・query | 完了 | 幾何snapshot、整数接続、厳密な物質位置式。 [#102](https://github.com/YUKIKEDA/UVLAS/pull/102) |
| 包絡と局所packet | 完了 | 位置・幾何包絡と局所適格性。 [#112](https://github.com/YUKIKEDA/UVLAS/pull/112) |
| P1核と全13変数の微分 | 完了 | 3基底の局所誘導核と全Jacobian。 [#129](https://github.com/YUKIKEDA/UVLAS/pull/129) |
| 共有辺・四辺形・全作用 | 完了 | 固定分割に従って全出現に散布し、核半径と強度トレースが一致する逆向きの解析辺を相殺する。値とJVPを独立面積分、VJPを双対内積で照合した。 [共有辺と四辺形の検証](../verification/IM-03/p1-shared-action-v1/verification.md) |
| 大域幾何と後縁所有 | 完了 | 大域接触・不適格の検査、alias運動学、後縁の放出所有。 [#153](https://github.com/YUKIKEDA/UVLAS/pull/153) |
| 既知流れと指定運動 | 完了 | 一様流・解析突風、指定移流、指定剛体運動を同じ厳密時刻と片側で評価する。全空力取付・構造支持の指定運動と過渡外力の周期節点や絶対イベントを確認し、節点をまたぐ区間を拒否する。離散空力取付では同じ姿勢から点の位置と速度を求める。再開時には保持した時刻表、既知流れ、後流モードと指定移流、選択した運動を現行入力と照合する。 [指定空力取付](../verification/IM-03/clocked-prescribed-aero-attachment-v1/verification.md)・[全イベントの検証](../verification/IM-03/all-motion-clock-events-v1/verification.md)・[型付き再開](../verification/IM-03/clocked-known-input-restore-v1/verification.md) |
| 循環と全微分 | 実装中 | 接触済みの離散照会点ごとに非一様な既知流れを評価し、実数の非貫通式を解く。同じ試行の原残差を再評価し、位置と既知流れの時刻への接線・転置を独立参照と照合した。PGの物質点写像、幾何と後流輸送も含む時間微分、接触切替え、複素更新が残る。 [非一様流での循環と微分の検証](../verification/IM-03/discrete-known-flow-circulation-v1/verification.md) |
| 全ベクトル荷重 | 実装中 | 同じ試行の循環と束縛辺の中点ごとの非一様な既知流れから、辺の力、節点力、合力、基準点モーメントを評価する。一つの準備済みPGでは、非貫通の位置・速度・接線を材料点で評価し、正準頂点の力を同じ形状行の転置で構造へ戻した。同一時計区間の固定刻みBDF2は循環の時間微分を三角形の荷重と現在時刻の接線・転置へ渡す。Caseの材料多項式を全パネルで調べ、層の片側を含む実PGの位置・速度・接線を求める。共有辺では面別荷重を保持し、指定した元頂点の形状行の転置で包絡つきの荷重を渡す。検査済みの元面所有表から対象Caseの全束縛面を照合し、元facetと局所頂点を材料行に結び付ける。登録別名は支持先が異なっても拘束後の行と指定運動の項が一致することを検査する。指定運動面と構造追従面の混在、J1の受理履歴が残る。 [全空力面の照合](../verification/IM-03/case-material-surface-model.md) [材料写像と元荷重の検証](../verification/IM-03/case-material-physical-v1.md) [曲面と荷重の双対仕事](../verification/IM-03/pg-aero-work-v1/verification.md)・[中点荷重と微分](../verification/IM-03/discrete-known-flow-loads-v1/verification.md) |
| 後流値・移流・試行 | 実装中 | 自由後流では束縛渦輪、先頭・旧後流、既知流れから粒子速度を評価する。既放出点を同じ物質IDと厳密な三時刻で照合し、自由・指定の固定刻みBDF2で位置残差と接線・転置を検証した。指定移流では候補位置も返す。新生点と再始動の既存点は、実際の親流速と試行流速を使う台形残差と接線・転置を持つ。新生点の後縁起点は同じ入力・時計・幾何の次のBDF2に引き継ぐ。全放出と退役の台帳、物理状態との一括確定、ファイル再開が残る。 [自由後流の粒子速度](../verification/IM-03/discrete-free-wake-velocity-v1/verification.md)・[固定刻みの後流輸送](../verification/IM-03/fixed-wake-transport-v1/verification.md) |
| J1状態への接続 | 残 | IM-04側とstep契約を共有し、BDF2と新生台形、輸送別の初期微分、受理済み循環の時間微分、全履歴、一回限りの状態確定を検証する。構造追従面の支持応答と、共通時計の型付き入力を含む受理済み履歴のファイル再開を接続する。強連成Newtonの実装はIM-04が担当 |
| 参照と数値予算 | 一部検証済み | 局所核に加え、有限翼2形状の非貫通残差と束縛循環第1列を独立参照と照合した。残る列、定常と非定常の荷重、収束、V-AERO・V-MULTIの全観測と数値予算を照合する。 [有限翼の独立参照](../verification/IM-03/finite-wing-case-v1/verification.md) |
| 費用と終了審査 | 残 | 全履歴の相互作用数、資源、並列時の結果一致とR-03を確認し、パッケージを審査する |

## IM-06 入出力

[#22](https://github.com/YUKIKEDA/UVLAS/issues/22) の終了審査まで12段階。取得した同じバイト列から、Nastranの基本系座標、支持方向、集中質量をSIで評価する。外部面から共通参照面への変換では、対応する面の全域にわたる位置と法線の誤差を包絡した。

| 段階 | 状態 | 到達点と残作業 |
|---|---|---|
| 厳格なscalar変換 | 完了 | 強いIDとIndexへの検査付き変換。 [#89](https://github.com/YUKIKEDA/UVLAS/pull/89) |
| parser資源計画 | 完了 | 固定JSON依存の例外経路を補修し、見積もりと受入条件を実装。 [依存補修](https://github.com/YUKIKEDA/UVLAS/pull/101)・[資源計画](https://github.com/YUKIKEDA/UVLAS/pull/104) |
| 所有文書 | 完了 | 元バイト列、digest、型、十進字句を保持。 [#104](https://github.com/YUKIKEDA/UVLAS/pull/104) |
| 全schema | 完了 | 全caseの型、入れ子、参照と相互条件の検証。 [#130](https://github.com/YUKIKEDA/UVLAS/pull/130) |
| SI・精度変換 | 完了 | 元の厳密値、SI倍率、成分別の精度管理。 [#154](https://github.com/YUKIKEDA/UVLAS/pull/154) |
| 参照snapshot | 実装中 | 宣言収集、パスの正規化、[単一ハンドルの照合](https://github.com/YUKIKEDA/UVLAS/pull/177)は統合済み。Windowsの根相対adapterは [#184](https://github.com/YUKIKEDA/UVLAS/pull/184)、Linuxの根相対adapterは [#195](https://github.com/YUKIKEDA/UVLAS/pull/195) で統合済み。[#202](https://github.com/YUKIKEDA/UVLAS/pull/202) で、マニフェストの元バイト列とdigestを単一のファイルハンドルから取得する部品を統合した。[#219](https://github.com/YUKIKEDA/UVLAS/pull/219) で、Fileは元マニフェストのハンドルと実体ID、Memoryは呼出し側が開いた根を検査済みケースとともに保持するようにした。[#229](https://github.com/YUKIKEDA/UVLAS/pull/229) で、参照ファイルを開く前に、宣言順と位置を保持して参照先を正規化・索引化する部品を統合した。[#236](https://github.com/YUKIKEDA/UVLAS/pull/236) で、指定された参照先1件の実体を内容の読み取り前に照合し、同じハンドルから取得した元バイト列を保持する部品を統合した。[#245](https://github.com/YUKIKEDA/UVLAS/pull/245) で、宣言順に依存をたどり、循環と共有経路の深さ上限を検査する内部処理を統合した。[#252](https://github.com/YUKIKEDA/UVLAS/pull/252) で、先頭の根宣言から取得したバイト列とダイジェストを照合し、依存走査へ渡す内部処理を統合した。マニフェスト直下の全8役割を同じ親の参照表で保持し、取得済みのバイト列からGmshメッシュを読める。宣言順に外部ファイルを取得し、元バイト列、依存辺、宣言位置、依存先を先に並べた完了順を所有する。形式ごとの子依存の抽出と、未認定の役割を含む入力の認定が残る。 [参照グラフの検証](../verification/IM-06/reference-graph.md) |
| 励振と既知入力 | 実装中 | 型付きSIケースから突風なしの一様流を準備し、定数・多項式・PCHIPの値と時間微分を評価する。指定した成分別の誤差予算でBinary64／TwoWordに投影した値を認定できる。参照snapshot、他の励振、J1への接続が残る |
| parametric入力とprepare | 残 | 型付きcaseを共通のRawModelへ渡す。物理準備には対応するIM-01／02／03の成果が必要 |
| Nastran | 実装中 | 認証済みデッキから基本系のGRID、面の元節点順・材料・厚さ・オフセット、SPCの方向、CONM2の重心と慣性を確認する。CBEAM／PBEAMは偏心後の梁の軸と長さ、各断面の断面量・物性をSIで保持する。RBE2は選択された成分の点運動と、独立GRIDへ戻る荷重・モーメントを評価する。対応表で元EID・MID・SPCADDをCaseへ結び付け、平面全域の対応を認定した後に一様全厚と明示材料値を照合する。一つのPCOMPが全域を覆う場合は層順・界面・明示材料値を法線の正逆で照合する。PSHELLは同じPIDの材料照合をまとめ、明示全厚を元面ごとの領域で比較する。CTRIA3の絶対厚さとCQUAD4の割合を含む角別指定は、元三角形の一次場または元四辺形の双一次場のまま、面積が正の共通領域でCaseに照合する。複数PCOMPの元面ごとに総厚を照合し、PIDと法線方向ごとに元層順・界面・材料を検証した。材料軸の方向、梁・剛体接続の全体組立てが残る。 [PSHELLの角別厚さ](../verification/IM-06/nastran-corner-thickness-pg.md) [PCOMPの部分領域](../verification/IM-06/nastran-pcomp-regional-pg.md) [PCOMPの層照合](../verification/IM-06/nastran-pcomp-pg-quantities.md) [認証済み対応表](../verification/IM-06/external-mesh-sidecar-role.md) [梁のSI入力](../verification/IM-06/nastran-si-beam-input.md)・[RBE2の仕事](../verification/IM-06/nastran-rbe2-projection.md)・[デッキの参照検査](../verification/IM-06/nastran-bulk-leaf-verification.md) |
| Gmsh | 一部実装 | 検査済みのスナップショットからMSH 4.1 ASCIIを読み、疎な節点タグ、所属、物理名、向きと元の座標字句を保持する。CLI は同じバイト列からメッシュ本体に子参照がないことを確かめる。明示した面集合を検査し、共有辺の向きと頂点のつながりを照合して面の連結成分と自由境界を返す。対応表から物理グループの元面をShellへ割り当て、平面では対応先となる同じCaseからPGを準備し、物理グループの明示全厚を共通細分の正面積領域で厳密照合する。非隣接面の自己交差、非一様な元厚さ、物理支持・接続が残る。 [Caseからの外部面照合](../verification/IM-06/case-external-pg-reference.md) [認証済み対応表](../verification/IM-06/external-mesh-sidecar-role.md) [入力契約](../specs/foundation/gmsh-mesh-input-contract.md)・[面位相](../verification/IM-06/gmsh-surface-topology.md) |
| 共通参照への変換 | 実装中 | Nastranの三角形・共面四角形と、明示した次数4以下のPG参照面の位置・法線誤差を全域で囲う。頂点が一致する曲面の内部差も検証した。SPCの方向、CONM2の慣性、梁断面とRBE2の作用をSIで照合した。Gmshの元双線形四角形でも位置・単位法線の差を囲う。GmshのfitとNastranの対応表で、全域の位置・法線上界をCase単位から換算した予算と照合し、観測IDを保持する。観測の独立評価、参照グラフの推定、変換前後の物理応答の照合が残る。 [変換予算の照合](../verification/IM-06/external-mesh-conversion-budget.md) [Gmshの曲面誤差](../verification/IM-06/gmsh-curved-pg-bounds.md) [曲面の誤差上界](../verification/IM-06/nastran-curved-pg-bounds.md)・[支持方向](../verification/IM-06/nastran-spc-component-directions.md)・[集中質量](../verification/IM-06/nastran-basic-mass-geometry.md) |
| 結果と診断 | 実装中 | 状態と来歴のJSON、確定観測のCSV、適格な可視化標本のVTU・PVDを内部で生成できる。VTU・PVDは層界面の両側と再始動区間を保つ。型検査・SI変換済みケースから、元入力、保持値、誤差半径、現時点で解決した既定値を出力し、その指紋を来歴に結び付けられる。物理回復から可視化標本を作り、物理準備で確定する設定とともに実行結果の保存へつなぐ工程が残る。 [解決済みケース](../verification/IM-06/resolved-case-document.md)・[観測CSV](../verification/IM-06/observation-csv.md)・[VTU・PVD](../verification/IM-06/visualization-series.md) |
| checkpoint容器 | 完了 | 型付きチャンクと指紋からマニフェストを作り、既存番号を保持したまま確定する。中断した保存を除外し、型、長さ、SHA-256、指紋が一致する最新番号を復元する。 [マニフェストと確定ディレクトリの検証](../verification/IM-06/checkpoint-manifest.md) |
| 実状態とresume | 残 | 全BDF2ステンシル、受理済みの微分・荷重・後流履歴、探索状態などを保存する。IM-04／05の実状態を使い、同版・同環境での再開結果のビット一致と非互換の拒否を確認する |
| CLI | 実装中 | `validate` は全スキーマと対応する外部参照を検査し、未対応の役割を拒否する。直下の外部ファイルは SHA-256 を照合し、MSH 4.1 ASCII の節とタグ、Gmsh・Nastran の版付き対応表を検査する。`inspect` は保存済み状態と来歴を照合し、来歴に指紋がある場合は解決済みケースも確認する。Python と検査処理を共有する。形式別の子依存の検査、prepare／run／resume、実行終了時の状態記録が残る。 [開発用CLI](../build.md#開発用cliの終了規約) |
| Python | 実装中 | ファイル・辞書の検査と結果確認に CLI と同じ C++ 処理を使い、検査済み入力と参照を所有する。物理準備・実行・再開、受理済み解析結果への接続が残る。 [入口の検証](../verification/IM-06/shared-case-frontends.md) |
| ランタイムと新規環境 | 実装中 | WindowsとLinuxで実wheelを空のPython環境へオフライン導入し、移設CLIと同じケース検査、外部参照の診断、所有権を確認した。Linuxは開発ツールなしのUbuntu Base、Windowsは開発PATHを外した環境で計測した。一般配布ABI、Windowsの独立した非開発環境、外部参照の検証成功と導入後の実行・再開が残る。 [実wheelの導入試験](../verification/IM-06/private-wheel-rehearsal.md) |
| 終了審査 | 残 | V-USE・V-STATEで3入力経路とCLI／Pythonの物理結果、両OSの導入・中断・再開を照合する |

schemaのLinux sanitizerの2試験は [#143](https://github.com/YUKIKEDA/UVLAS/pull/143) の修正を含むmainで再計測し、49.51秒と53.22秒で通過した。Linux sanitizerだけ300秒、Linux Debug・ReleaseとWindows Debug・Releaseは60秒の登録を確認した。 [#142](https://github.com/YUKIKEDA/UVLAS/issues/142) 数値表では同じ段階の子Issueを重ねて加えない。

## 初回実務適用までの位置

IM-00の共通基盤は終了審査を通過した。IM-01の成果をIM-02の材料・梁・接続・構造全体へ渡し、IM-02とIM-03の成果をIM-04の強連成へつなぐ。その後にIM-05の安定性・探索、IM-07の実務規模の数値・性能検証、IM-08の実験・実務受入がある。IM-06は並行しながら、実行と再開の段階で揃った物理状態を接続する。 [全体の依存とゲート](implementation-start-review.md#3-後続パッケージと依存)
