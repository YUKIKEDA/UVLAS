# IM-02・IM-03・IM-06の進捗と残作業

2026年11月28日更新。統合状況は `main` の `9a6d8ee3` 時点。完了数は統合済みの成果に基づく。

| パッケージ | 完了と残り | 次に通す検証 |
|---|---|---|
| [IM-02 構造](#im-02-構造) | 4成果群すべてに残作業がある。2 Pa載荷では、試作の4配置で接線の正定値性を認定し、最終候補の補正を全216行で照合した | Case入口から載荷平衡を確定する。同じ平衡点で反力・応力・質量・対象固有対を独立照合する |
| [IM-03 空力](#im-03-空力) | 13段階中7段階完了、残り6段階。翼格子・非定常荷重・小コア監査を統合済み | 自由後流の有限窓と観測量ごとの誤差予算を定める。循環・荷重・後流をJ1の同じ受理状態に結び付ける |
| [IM-06 入出力](#im-06-入出力) | 18段階中6段階完了、残り12段階。観測参照契約は統合済み。CLI／Pythonの共通検査は3入力経路を通し、レビュー中 | 確定した解析結果を参照と照合する。解析状態の保存・再開にはIM-04／05の受理状態が必要 |

段階と成果群は[作業分解](next-package-work-breakdown.md)に従う。規模が異なるため、完了数から残りの工数は見積もっていない。前工程のIM-01は12段階すべてが完了し、[終了審査](../verification/IM-01/exit-review.md)を通過した。

## 並行作業と先行条件

| パッケージ | 先行工程を待たずに進めている作業 | 先行工程が必要な検証 |
|---|---|---|
| IM-02 | 2 Pa載荷で、2回目のNewton補正まで認定した。候補更新で消費する演算数を実際の状態と補正から調べている。並行して、R-05矩形ケースの3入力経路で固定チャートと側面牽引の準備を通した | 載荷平衡を確定し、その状態で先端変位・支持反力・応力・質量・対象固有対を独立参照と照合する |
| IM-03 | 凍結ケースの6ステップで、有限後流の切捨てと保持点の位置差が循環・力に与える影響を分けて照合した。既知の2軌道間では三角形ごとの誤差を連続上界で囲えた | J1の受理状態で、時間発展する後流と構造応答を含めた累積誤差を検証する |
| IM-06 | 観測参照を検査するCLI／Pythonの共通処理でMSVCの警告を修正した。Linuxの観測参照とAPI試験3件が通過し、修正後のWindows CIは確認待ち。 [観測参照PR](https://github.com/YUKIKEDA/UVLAS/pull/420) | R-05の3経路の物理結果比較にはIM-02の載荷平衡が必要。実状態を保存して再開する試験にはIM-04／05の受理状態が必要 |

## IM-01 構造

[#17](https://github.com/YUKIKEDA/UVLAS/issues/17) は固定疎チャートまで実装した。複数の指定支持を実際のケースから組み立て、R-02の全厚位置差を実空間で照合した。V-PGの参照と拘束を検証し、IM-02への引継ぎ審査を通過した。

IM-02 と R-05 では、この参照パッチを使って Case の辺面力を[全層・全厚の側面へ組み立てる](../specs/structure/pg-case-lateral-traction-contract.md)。

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
| 終了審査 | 完了 | 凍結34 Caseの36 Shellから製品の物理参照を準備した。R4の28配置とR5の4配置で、位置、法線、全厚、層の片側と物理測度を独立式と照合した。32配置では指定有限則の合同誤差が $10^{-20}$ 以下となり、385,968成分で必須場の厳密一致を確認した。凍結した二曲面の整合・非整合接続をCaseから準備し、両OSで全元行の拘束、再構成、転置の仕事と反力を照合した。非整合の元行は288、rankは171。同じCase入口で時変支持を検査し、凍結曲面と平面対照で両側時刻の位置・速度・加速度10,368成分、原拘束と仕事を独立照合した。凍結R5の全4配置で、この同じ所有パッチの参照位置、形状と厳密な物理勾配、有限点と重み、元拘束を引継ぎ入口から照合した。選んだ有限則では、物理勾配のTwoWord射影が試験の絶対誤差枠を超える対照も確認した。Caseでは明示点列と正の積則を選べる。平面対照で独立な2,793点と支持チャートが一致した。曲面R5の17×17×3積則は全16,473点のGram準備と2点のTwoWord射影を終えた。凍結R5の全4配置・480点で非アフィン2場の尺度付き位置・物理勾配、全元行、$J$ と $J^\mathsf{T}$ の仕事を独立照合した。TwoWord作用の最大微分差は $2.61\times10^{-32}$。曲面の積則の物理作用への積分誤差はIM-02の判定に引き継ぐ。 [R5の物理参照](../verification/IM-01/pg-case-reference-handoff/physical-reference.md) 局所独立化の子Issueは終了し、オーナーレビューで引継ぎが承認された。 [終了審査](../verification/IM-01/exit-review.md) [Caseからの物理参照](../verification/IM-01/vpg-case-product-reference.md)・[R4の物理幾何](../verification/IM-01/pg-r4-domain-geometry.md)・[R5の切断境界と層](../verification/IM-01/pg-r5-domain-geometry.md)・[R4の二曲面接続](../verification/IM-01/pg-r4-curved-interface/verification.md)・[有限則と資源の照合](../verification/IM-01/pg-vpg-finite-conditioning/verification.md)・[Caseの全厚拘束](../verification/IM-01/pg-case-full-chart/verification.md)・[動く支持の全厚拘束](../verification/IM-01/pg-case-moving-chart/verification.md)・[R5の参照引継ぎ](../verification/IM-01/pg-case-reference-handoff/verification.md) |

## IM-02 構造

[#18](https://github.com/YUKIKEDA/UVLAS/issues/18) の4成果群を追跡する。IM-01で認定したCaseの参照基底と元拘束を、同じCaseの物理作用・求解・回復に使う。

| 成果群 | 到達点 | 残作業と判定条件 |
|---|---|---|
| Case積層からシェル物理作用 | Caseのシェル材料と局所の物理作用を統合した。 [#388](https://github.com/YUKIKEDA/UVLAS/pull/388) | 採用した全3D材料則、初期ひずみ、密度、積層と層界面の両側を通し、エネルギー・内力・接線・整合質量を独立照合する |
| 梁の力学と断面 | Nastran入力から梁の軸と断面量をSIで取得できる。 [梁入力の検証](../verification/IM-06/nastran-si-beam-input.md) | 採用SE(3)梁の内力・接線・質量・慣性を組み込み、V-BEAMで独立参照と照合する |
| 支持・接続・作用 | Caseの辺面力を全層・全厚の側面へ組み立てられる。 [#397](https://github.com/YUKIKEDA/UVLAS/pull/397) | 支持・非線形接続・圧力・粘性を同じ物理状態に接続する。V-JOINT／V-LOADで元行反力、仕事、必要微分を照合する |
| 同じ平衡での全体系 | 平面シェルの静的平衡の仕様を統合した。ローカル試作では無載荷対照と、2 Pa載荷の最終候補の残差・補正・支持レンチを確認した。 [#404](https://github.com/YUKIKEDA/UVLAS/pull/404) | Case入口から載荷平衡を確定する。同じ状態で元の力の釣り合い、反力、層別応力、切断面力、整合質量、対象固有対を照合する。V-PG／V-RECOVERは凍結R4の28配置とR5の4配置を含む |

2 Pa載荷の最初の試行を受理した状態で、8体積すべての接線を誤差付きで組み立てた。独立の整数演算で仮想仕事を積分し、後退誤差の認定が残っていた106行の右辺が厳密にゼロであることを確認した。その106成分をゼロとした補正候補は、全216行で絶対残差と成分別後退誤差の条件を満たした。

この補正を刻み1で適用した候補は、独立計算でポテンシャルの十分減少条件と全8体積の向きの保持を確認した。自由残差の最大絶対値は約 $1.88\times10^{-1}$ N から $4.05\times10^{-5}$ N に下がった。

更新後の接線で後続2回の補正を求め、独立計算で各回とも全216行の絶対残差と成分別後退誤差を確認した。いずれも刻み1で十分減少条件を満たし、全8体積の向きを保持した。自由残差の最大絶対値は約 $1.87\times10^{-12}$ Nを経て $5.20\times10^{-27}$ Nまで下がり、最後の指定点変位差は約 $1.09\times10^{-15}$ mだった。同じ候補で、独立に積分した元力327成分から支持レンチを回復し、外力の合力・モーメントと釣り合うことを確認した。

保存した4配置の接線全成分を、求積と丸めの誤差を含めて囲い、製品の疎行列検査で正定値性を認定した。誤差幅を広げた対照では認定を拒否した。最終候補で求め直した補正は最大約 $6.83\times10^{-32}$ m。接線作用を独立に積分したところ、補正式の最大絶対残差は約 $2.14\times10^{-41}$ N、成分別後退誤差の上界は約 $1.22\times10^{-14}$ で、全216行が試作の許容値を満たした。同じCaseの整合質量は、総質量0.1 kgと参照位置の一次・二次モーメントを平板の解析積分と厳密に照合した。

この載荷対照では、右辺が厳密に零だった106自由行すべてで、持ち上げた方向が厚さ中央面に対して反対称で、外力の仮想仕事も厳密に零だった。配置と方向の対称性を全8体積の係数で照合した。ローカルの構造実装では、元力から得た区間が零を含む場合にこの条件を調べ、証明できた右辺を厳密に零とする。力の区間を広げたアフィン載荷対照で、全216行のうち該当する106行に適用できることを確かめた。

後続試行の接線は、体積ごとに組み立てた誤差包絡を固定チャートで縮約し、補正の検査に渡す。ローカルのアフィン2 Pa載荷対照では、全8体積の接線を組み立て、縮約接線の正定値性と全216行の補正を認定した。

保存した最初の非アフィン載荷配置では、元力の包絡に256ビットの二進格子を使い、保存済みの接線包絡と組み合わせて全216行の補正を認定した。接線の計算でも、配置の勾配を丸めたときの誤差を接線包絡に加える経路を組み込んだ。実体積1つで非零42成分と零36成分を独立参照と照合した。

参照配置から始める2 Pa載荷の通常求解でも、組み立てた接線包絡を使って2回目の補正を認定した。処理回数の診断用上限を2,400万回とした実行は、補正後の候補を作る工程で資源上限に達した。

R-05の矩形について、パラメトリック、Gmsh、Nastranの3経路で固定チャート、材料、全厚の側面牽引を準備した。外部経路では元データと対応表を認証した。パラメトリック経路の通常求解を参照配置から実行し、2回目の補正と観測評価まで確認した。3経路で平衡を確定し、同じ確定状態を先端変位・支持反力の独立照合と固有振動へ渡す工程が残る。R-05に続いて、梁、非線形接続・圧力・粘性、およびV-PG／V-RECOVERの全配置を進める。

## IM-03 空力

[#19](https://github.com/YUKIKEDA/UVLAS/issues/19) の終了審査まで6段階。共有辺・四辺形の誘導作用を循環計算へ渡し、既知流れと指定運動を共通時計で評価する。循環とその微分が、全荷重・後流へ進む共通の先行条件。

水平翼と斜翼の正式格子を使った調和循環の検証を統合した。 [#396](https://github.com/YUKIKEDA/UVLAS/pull/396)

翼の前後縁から渦格子を準備する [#398](https://github.com/YUKIKEDA/UVLAS/pull/398)、非定常荷重の [#402](https://github.com/YUKIKEDA/UVLAS/pull/402)、小コア監査の [#409](https://github.com/YUKIKEDA/UVLAS/pull/409) は、固定容量バッファの修正を含めて `891037e0` に統合した。修正後のLinux Release回帰は各ブランチで通過した。既存3試験のDebug実行上限の変更も統合した。 [#416](https://github.com/YUKIKEDA/UVLAS/pull/416)

非定常荷重は、水平翼と斜翼の各3格子について、4種類の励振で翼別の力・モーメント48件を独立参照と照合した。

指定直線後流の試作では、128 mで打ち切った後流と半無限後流の差を評価した。遠方作用の式に左右の縦辺と横辺の相殺を入れ、水平翼・斜翼の各3格子、4励振の翼別鉛直力とピッチングモーメントの計96観測で、相対誤差を0.382%未満に囲った。自由後流の退役源と観測量ごとの誤差予算が残る。循環・荷重・後流の一括確定と再開はJ1の受理状態につなぐ。

凍結した自由後流の親状態を再生成し、配置と循環を固定して最古の1区間を落とした。面積分の解析上界を使い、非貫通点・保持物質点・束縛荷重点の8照会で凍結実装の境界式による直接誘導の差を囲んだ。次に、同じ指定翼運動で、自由移流のEuler検証経路を成長後流と2行の有限窓で各6回進めた。初回の力の差は、循環を固定して退役区間だけを除くと約 $1.93\times10^{-4}$ N、循環とその時間微分を更新すると約0.224 Nだった。2回目から共有物質点の位置も分かれた。採用時間積分で循環、移流、構造運動を通じた誤差を認定する作業が残る。

共有する後流形状が一致する初回では、凍結した1パネル対照の離散非貫通行列を有理数で解いた。退役による循環差は約 $6.17\times10^{-3}$ m²/sで、退役源の直接誘導から求めた上界は約 $2.39\times10^{-2}$ m²/sだった。同じ親履歴からBDF2で時間微分を評価し、鉛直の非定常力だけで約0.208 Nの差が生じることを独立に照合した。保持点が移動した後の行列変化と、構造応答を含む誤差伝播が残る。

同じ1パネル対照の全6回で、離散循環差を、退役源、保持源による右辺の差、先頭後流の形状変化による行列の差から再現した。6回目の循環差は約0.0215 m²/sで、保持源による右辺の差はそのうち約0.00238 m²/sだった。保存した2配置の間で三角形を補間し、各作用の差を囲ったときの循環差の上界は約0.0931 m²/s。BDF2の3時刻履歴から得た6回目の鉛直非定常力差は約0.0376 Nだった。構造応答を含む採用J1の軌道に沿って誤差を伝播させる工程が残る。

固定面・指定運動面で元面と後縁を所有するため、Case の前後縁から渦格子と非貫通点を作る写像を定めた。 [前縁と後縁の対応規約](../specs/aero/wing-edge-pairing-contract.md)

| 段階 | 状態 | 到達点と残作業 |
|---|---|---|
| 厳密な局所述語 | 完了 | 退化、内外、距離とfar判定。 [#90](https://github.com/YUKIKEDA/UVLAS/pull/90) |
| 所有幾何・ID・query | 完了 | 幾何snapshot、整数接続、厳密な物質位置式。 [#102](https://github.com/YUKIKEDA/UVLAS/pull/102) |
| 包絡と局所packet | 完了 | 位置・幾何包絡と局所適格性。 [#112](https://github.com/YUKIKEDA/UVLAS/pull/112) |
| P1核と全13変数の微分 | 完了 | 3基底の局所誘導核と全Jacobian。 [#129](https://github.com/YUKIKEDA/UVLAS/pull/129) |
| 共有辺・四辺形・全作用 | 完了 | 固定分割に従って全出現に散布し、核半径と強度トレースが一致する逆向きの解析辺を相殺する。値とJVPを独立面積分、VJPを双対内積で照合した。 [共有辺と四辺形の検証](../verification/IM-03/p1-shared-action-v1/verification.md) |
| 大域幾何と後縁所有 | 完了 | 大域接触・不適格の検査、alias運動学、後縁の放出所有。 [#153](https://github.com/YUKIKEDA/UVLAS/pull/153) |
| 既知流れと指定運動 | 完了 | 一様流・解析突風、指定移流、指定剛体運動を同じ厳密時刻と片側で評価する。全空力取付・構造支持の指定運動と過渡外力の周期節点や絶対イベントを確認し、節点をまたぐ区間を拒否する。離散空力取付では同じ姿勢から点の位置と速度を求める。再開時には保持した時刻表、既知流れ、後流モードと指定移流、選択した運動を現行入力と照合する。 [指定空力取付](../verification/IM-03/clocked-prescribed-aero-attachment-v1/verification.md)・[全イベントの検証](../verification/IM-03/all-motion-clock-events-v1/verification.md)・[型付き再開](../verification/IM-03/clocked-known-input-restore-v1/verification.md) |
| 循環と全微分 | 実装中 | 接触済みの離散照会点ごとに非一様な既知流れを評価し、実数の非貫通式を解く。同じ試行の原残差を再評価し、位置と既知流れの時刻への接線・転置を独立参照と照合した。循環の求解と接線・転置で照会ごとの作業領域を再利用する。水平翼M2/N4の循環・荷重の累積確保はGCC 14で35.0%減少した。 [作業領域と凍結有限翼の照合](../verification/IM-03/circulation-row-workspace.md) PGの物質点写像、幾何と後流輸送も含む時間微分、接触切替え、複素更新が残る。 [非一様流での循環と微分の検証](../verification/IM-03/discrete-known-flow-circulation-v1/verification.md) |
| 全ベクトル荷重 | 実装中 | 同じ試行の循環と束縛辺の中点ごとの非一様な既知流れから、辺の力、節点力、合力、基準点モーメントを評価する。一つの準備済みPGでは、非貫通の位置・速度・接線を材料点で評価し、正準頂点の力を同じ形状行の転置で構造へ戻した。同一時計区間の固定刻みBDF2は循環の時間微分を三角形の荷重と現在時刻の接線・転置へ渡す。Caseの材料多項式を全パネルで調べ、層の片側を含む実PGの位置・速度・接線を求める。共有辺では面別荷重を保持し、指定した元頂点の形状行の転置で包絡つきの荷重を渡す。検査済みの元面所有表から対象Caseの全束縛面を照合し、元facetと局所頂点を材料行に結び付ける。登録別名は支持先が異なっても拘束後の行と指定運動の項が一致することを検査する。荷重中点では自己辺の除外と後流寄与を保ちながら、照会間で作業領域を再利用する。 [有限翼の計測](../verification/IM-03/circulation-midpoint-workspace.md) 指定運動面と構造追従面の混在、J1の受理履歴が残る。 [全空力面の照合](../verification/IM-03/case-material-surface-model.md) [材料写像と元荷重の検証](../verification/IM-03/case-material-physical-v1.md) [曲面と荷重の双対仕事](../verification/IM-03/pg-aero-work-v1/verification.md)・[中点荷重と微分](../verification/IM-03/discrete-known-flow-loads-v1/verification.md) |
| 後流値・移流・試行 | 実装中 | 自由後流では束縛渦輪、先頭・旧後流、既知流れから粒子速度を評価する。既放出点を同じ物質IDと厳密な三時刻で照合し、自由・指定の固定刻みBDF2で位置残差と接線・転置を検証した。指定移流では候補位置も返す。新生点と再始動の既存点は、実際の親流速と試行流速を使う台形残差と接線・転置を持つ。新生点の後縁起点は同じ入力・時計・幾何の次のBDF2に引き継ぐ。成長後流は全正準後縁を放出して既存の物質キーを保持し、Case行上限で停止する。有限窓は固定刻みから行数を求め、試行で退役する最古行を返す。同じ物質帯で親の有向後縁から放出強度を引き継ぎ、保持と退役の各強度を両語と帯IDで所有する。一つのCase試行で新生点を台形、区間途中の保持点をBDF2へ送る。二面の試行で全行の残差と入力順のシードを使う転置を照合した。有限窓の保持距離・誤差の認定、物理状態との一括確定、ファイル再開が残る。 [自由後流の粒子速度](../verification/IM-03/discrete-free-wake-velocity-v1/verification.md)・[固定刻みの後流輸送](../verification/IM-03/fixed-wake-transport-v1/verification.md) |
| J1状態への接続 | 残 | IM-04側とstep契約を共有し、BDF2と新生台形、輸送別の初期微分、受理済み循環の時間微分、全履歴、一回限りの状態確定を検証する。構造追従面の支持応答と、共通時計の型付き入力を含む受理済み履歴のファイル再開を接続する。強連成Newtonの実装はIM-04が担当 |
| 参照と数値予算 | 一部検証済み | 凍結有限翼2形状の各5迎角で、有限コアの定常循環、非貫通の原残差、風軸の合力と基準点まわりのモーメント、固定20帯の荷重と両刻みの傾斜を独立参照と照合した。同じ格子の+1度で後流長32〜1024 mと正コア半径の差を分離し、同じ閉端の独立理想核とも照合した。凍結二翼の翼弦2列・翼幅12列では48×48の全係数と翼別荷重を独立参照と照合し、節点力の仕事が翼別合力・モーメントの仕事と一致した。同じ有限翼で辺の自己作用を除いた力、20帯の線分荷重とTrefftz抗力を独立比較した。水平翼・斜翼のM128/N192で、近傍抗力とTrefftz面抗力の差、独立な凍結極限との差はどちらも1%以内だった。格子・時刻・後流・コアを含む数値予算と非定常観測が残る。 [正式格子の抗力と原残差](../verification/IM-03/formal-wing-assembly-cost.md) [線分荷重と遠方抗力](../verification/IM-03/finite-wing-span-lines-v1/verification.md) [二翼の定常検証](../verification/IM-03/multiple-surface-steady-v1/verification.md) [有限翼の5迎角と20帯](../verification/IM-03/finite-wing-polars-v1/verification.md)・[後流長とコア半径](../verification/IM-03/finite-wing-wake-core-v1/verification.md) |
| 費用と終了審査 | 実装中 | 2・4・8行の自由後流で、粒子速度・微分・転置の結果指紋と相互作用数を維持しながら累積確保を減らした。遠方の共通ジェットを使った8行の輸送値・接線・転置は最短時間が15.63%短縮した。凍結参照との照合は両翼のM4/N12と水平翼M8/N24が通った。正式格子M128/N192では両翼の全幾何、接触・後縁所有、照会3件の作業域と仕事量の予報を確認した。水平翼・斜翼の正式格子では全24,576行の原残差と抗力の判定まで、それぞれ915.76秒と960.11秒で通った。実務規模での全履歴の費用、並列時の結果一致、終了審査が残る。 [正式格子の費用](../verification/IM-03/formal-wing-assembly-cost.md) [幅方向細分と自由後流の検証](../verification/IM-03/finite-wing-wide-admission.md)・[正式格子の事前見積もり](../verification/IM-03/finite-wing-interval-rounding.md) |

## IM-06 入出力

[#22](https://github.com/YUKIKEDA/UVLAS/issues/22) の終了審査まで12段階。取得した同じバイト列から、Nastranの基本系座標、支持方向、集中質量をSIで評価する。外部面から共通参照面への変換では、対応する面の全域にわたる位置と法線の誤差を包絡した。

Caseのファイル入力では親ディレクトリ、Python辞書入力では明示した参照ルートを起点に参照先を検査する。 [#403](https://github.com/YUKIKEDA/UVLAS/pull/403)

観測参照JSON v1の契約と例を統合した。共通Case検査への接続はレビュー中。Linux Releaseの実CLIとPythonでparametric・Gmsh・Nastranの3経路を通し、元バイト列、条件、記述子、単位、成分の不一致を拒否した。C++では所有した値と来歴、累積資源上限、取消後の再利用を確認した。確定した解析結果と参照の比較、prepare／run／resumeと両OSの新規環境での実行・再開の検証が残る。 [契約 #408](https://github.com/YUKIKEDA/UVLAS/pull/408)・[共通検査 #420](https://github.com/YUKIKEDA/UVLAS/pull/420)

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
| Nastran | 実装中 | 認証済みデッキから基本系のGRID、面の元節点順・材料・厚さ・オフセット、SPCの方向、CONM2の重心と慣性を確認する。CBEAM／PBEAMは偏心後の梁の軸と長さ、各断面の断面量・物性をSIで保持する。RBE2は選択された成分の点運動と、独立GRIDへ戻る荷重・モーメントを評価する。対応表で元EID・MID・SPCADDをCaseへ結び付け、平面全域の対応を認定した後に一様全厚と明示材料値を照合する。一つのPCOMPが全域を覆う場合は層順・界面・明示材料値を法線の正逆で照合する。PSHELLは同じPIDの材料照合をまとめ、明示全厚を元面ごとの領域で比較する。CTRIA3の絶対厚さとCQUAD4の割合を含む角別指定は、元三角形の一次場または元四辺形の双一次場のまま、面積が正の共通領域でCaseに照合する。複数PCOMPの元面ごとに総厚を照合し、PIDと法線方向ごとに元層順・界面・材料を検証した。一定全厚の平面では、元のTHETAまたは矩形座標系の軸と一定Case材料軸を、厳密な内積・外積の符号と角度で照合した。認証済み対応表で宣言した全Shellを同じCaseから準備し、厚さ・積層・軸を所有する結果を一括で返す。曲面、変厚、位置依存の材料軸と梁・剛体接続の全体組立てが残る。 [全Shellの検証](../verification/IM-06/case-nastran-planar-pg-set.md) [材料軸の厳密照合](../verification/IM-06/nastran-pg-material-angle.md) [PSHELLの角別厚さ](../verification/IM-06/nastran-corner-thickness-pg.md) [PCOMPの部分領域](../verification/IM-06/nastran-pcomp-regional-pg.md) [PCOMPの層照合](../verification/IM-06/nastran-pcomp-pg-quantities.md) [認証済み対応表](../verification/IM-06/external-mesh-sidecar-role.md) [梁のSI入力](../verification/IM-06/nastran-si-beam-input.md)・[RBE2の仕事](../verification/IM-06/nastran-rbe2-projection.md)・[デッキの参照検査](../verification/IM-06/nastran-bulk-leaf-verification.md) |
| Gmsh | 一部実装 | 検査済みのスナップショットからMSH 4.1 ASCIIを読み、疎な節点タグ、所属、物理名、向きと元の座標字句を保持する。CLI は同じバイト列からメッシュ本体に子参照がないことを確かめる。明示した面集合を検査し、共有辺の向きと頂点のつながりを照合して面の連結成分と自由境界を返す。認証済みの対応表から全Shellを一括準備する。平面と曲面の全域を同じCaseのfit予算で照合し、物理グループが指定する全厚は共通細分で厳密照合する。後段のShellの不一致は元面を指し、一括結果を返さない。元一次線と名前付き外周が全域で一致することを確認し、欠け、重複、別外周、内辺を拒否した。面のfit誤差が正のときも、元線とPG境界の一致を辺上で検査する。支持と接着接続の対応を照合し、元線と接する元面を、使用したPG参照とともに構造チャートへ渡す。`supplied` の参照ファイル、自己交差、非一様な元厚さが残る。 [外部メッシュからの構造チャート](../verification/IM-06/external-case-pg-chart.md) [Gmsh支持境界の照合](../verification/IM-06/gmsh-case-support-boundary.md) [Gmsh全Shellの実ファイル検証](../verification/IM-06/case-gmsh-pg-set.md) [Caseからの外部面照合](../verification/IM-06/case-external-pg-reference.md) [認証済み対応表](../verification/IM-06/external-mesh-sidecar-role.md) [入力契約](../specs/foundation/gmsh-mesh-input-contract.md)・[面位相](../verification/IM-06/gmsh-surface-topology.md) |
| 共通参照への変換 | 実装中 | Nastranの三角形・共面四角形と、明示した次数4以下のPG参照面の位置・法線誤差を全域で囲う。頂点が一致する曲面の内部差も検証した。SPCの方向、CONM2の慣性、梁断面とRBE2の作用をSIで照合した。Gmshの元双線形四角形でも位置・単位法線の差を囲う。GmshのfitとNastranの対応表で、全域の位置・法線上界をCase単位から換算した予算と照合し、観測IDを保持する。同じCaseを使った実ファイルの平面対照では、両形式とも上面積4 m²、体積1 m³で独立式に一致した。厚さを変えた対照は元要素と厚さの元欄を示して拒否した。一般の観測、参照グラフの推定、変換前後の物理応答の照合が残る。 [形式間比較](../verification/IM-06/external-shell-cross-format.md) [変換予算の照合](../verification/IM-06/external-mesh-conversion-budget.md) [Gmshの曲面誤差](../verification/IM-06/gmsh-curved-pg-bounds.md) [曲面の誤差上界](../verification/IM-06/nastran-curved-pg-bounds.md)・[支持方向](../verification/IM-06/nastran-spc-component-directions.md)・[集中質量](../verification/IM-06/nastran-basic-mass-geometry.md) |
| 結果と診断 | 実装中 | 状態と来歴のJSON、確定観測のCSV、適格な可視化標本のVTU・PVDを内部で生成できる。VTU・PVDは層界面の両側と再始動区間を保つ。型検査・SI変換済みケースから、元入力、保持値、誤差半径、現時点で解決した既定値を出力し、その指紋を来歴に結び付けられる。物理回復から可視化標本を作り、物理準備で確定する設定とともに実行結果の保存へつなぐ工程が残る。 [解決済みケース](../verification/IM-06/resolved-case-document.md)・[観測CSV](../verification/IM-06/observation-csv.md)・[VTU・PVD](../verification/IM-06/visualization-series.md) |
| checkpoint容器 | 完了 | 型付きチャンクと指紋からマニフェストを作り、既存番号を保持したまま確定する。中断した保存を除外し、型、長さ、SHA-256、指紋が一致する最新番号を復元する。 [マニフェストと確定ディレクトリの検証](../verification/IM-06/checkpoint-manifest.md) |
| 実状態とresume | 残 | 全BDF2ステンシル、受理済みの微分・荷重・後流履歴、探索状態などを保存する。IM-04／05の実状態を使い、同版・同環境での再開結果のビット一致と非互換の拒否を確認する |
| CLI | 実装中 | `validate` は全スキーマと対応する外部参照を検査し、未対応の役割を拒否する。直下の外部ファイルは SHA-256 を照合し、MSH 4.1 ASCII の節とタグ、Gmsh・Nastran の版付き対応表を検査する。`inspect` は保存済み状態と来歴を照合し、来歴に指紋がある場合は解決済みケースも確認する。Python と検査処理を共有する。形式別の子依存の検査、prepare／run／resume、実行終了時の状態記録が残る。 [開発用CLI](../build.md#開発用cliの終了規約) |
| Python | 実装中 | ファイル・辞書の検査と結果確認に CLI と同じ C++ 処理を使い、検査済み入力と参照を所有する。物理準備・実行・再開、受理済み解析結果への接続が残る。 [入口の検証](../verification/IM-06/shared-case-frontends.md) |
| ランタイムと新規環境 | 実装中 | WindowsとLinuxで実wheelを空のPython環境へオフライン導入した。移設CLIとPythonは実BDF・実MSHと両方の対応表を含む同一Caseを読み込み、改変や対応先の不一致を拒否し、復元後に再び成功した。Linuxは開発ツールなしのUbuntu Base、Windowsは開発PATHを外した環境で確認した。一般配布ABI、独立したWindows環境、導入後の求解と実状態からの再開が残る。 [実wheelの導入試験](../verification/IM-06/private-wheel-rehearsal.md)・[外部参照の成功と拒否](../verification/IM-06/private-wheel-external-references.md) |
| 終了審査 | 残 | V-USE・V-STATEで3入力経路とCLI／Pythonの物理結果、両OSの導入・中断・再開を照合する |

schemaのLinux sanitizerの2試験は [#143](https://github.com/YUKIKEDA/UVLAS/pull/143) の修正を含むmainで再計測し、49.51秒と53.22秒で通過した。Linux sanitizerだけ300秒、Linux Debug・ReleaseとWindows Debug・Releaseは60秒の登録を確認した。 [#142](https://github.com/YUKIKEDA/UVLAS/issues/142) 数値表では同じ段階の子Issueを重ねて加えない。

## 初回実務適用までの位置

IM-00の共通基盤は終了審査を通過した。IM-01の成果をIM-02の材料・梁・接続・構造全体へ渡し、IM-02とIM-03の成果をIM-04の強連成へつなぐ。その後にIM-05の安定性・探索、IM-07の実務規模の数値・性能検証、IM-08の実験・実務受入がある。IM-06は並行しながら、実行と再開の段階で揃った物理状態を接続する。 [全体の依存とゲート](implementation-start-review.md#3-後続パッケージと依存)
