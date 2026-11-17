# 全schemaの設計検査と実装への引継ぎ

本書はIssue #110の設計審査表である。[辞書](case-schema-fields.md)と[JSONの最小fragment集](case-schema-examples.json)を照合する。fragment集とroot全体のschema例を分けている。schema_examplesの5例は空力単体の定常/指定運動過渡、構造単体静的、連成過渡/フラッタに全rootを当てた設計例である。実ファイル・物理prepareを備えたcaseではない。これをCLIの有効case、数値通過例、実務入力として配布しない。

## DC-01 各variantの最小例

fragment名の `Type.tag` は辞書の対応variantである。単一型は型名を使う。JSONはnative integer、decimal number、string、boolを区別し、shapeを明示する。例のSHA-256は設計検証記録に保持する。全fragmentを同一caseへ結合してはならない。例えば各Material例のid=1は各局所例の定義であり、同一modelの重複IDを許す意味ではない。

| 対象 | fragment名 | 局所registry/条件 |
|---|---|---|
| root構成 | Units、Geometry.pg、Material.isotropic、Shell、Flow、Analysis.static、Observation.point、Scale | schema Caseの9必須keyに同じ型を配置する。numerics/resources/outputの全keyは辞書通り別途必須 |
| 元幾何 | Geometry.pg/surface_graph/nastran/gmsh、ReferenceGeometry.supplied/fit | importer参照のhashは字句例の64桁零であり、実ファイル一致を主張しない |
| 場 | ReferenceField.constant/polynomial/expression、TimeSeries.constant/polynomial/pchip/expression | shape nを各例で1とする。射影軸はn3、明示軸はn9 |
| 範囲外 | Outside.reject/hold/periodic | 周期例は値0の2節点tableと組合せ、一階微分0で整合 |
| 軸・運動 | Axis.projected/explicit、Motion | explicitは単位行列、projectedはflat graphのx方向。Motionは定数零でC2 |
| 材料 | Material.isotropic/orthotropic3d/effective3d、Applicability.unevaluated/declared | E=G=ρ=1、ν=0、SPD単位行列。実材料受入ではない |
| 層・シェル | Layer、Shell | ζ=[−1,1]、1層、MAT-D08の全3D。幾何厚さH=1なので実法線厚さ[−1/2,1/2] |
| 梁 | Beam、BeamCoordinate.arclength/scaled、Section.matrix/circle/rectangle | 2節点の位置0/1、自由反り。matrixはK6=I。circle/rectangleの明示せん断面積とtorsionを保持 |
| 接続 | Joint.adhesive/rigid_attachment/weld | 各対象IDと参照geometryを局所registryへ置く。adhesiveは同じ参照面のη写像を両側に明示。点列fitではない |
| 支持 | Support.constraint/trace_constraint/elastic | constraintは並進1成分、elasticは明示K/D。shell取付frameを必要とする回転成分はprepareに同じ物理接線を渡す |
| 質量 | Mass.particle/rigid、Ledger.add/replace | replace例は別Mass id2を局所registryに置く。粒子には慣性を足さない |
| 荷重 | Load.point_force/couple/gravity/traction/pressure/body_force/current_pressure | coupleはbeam対象、面/体積測度を区別。current_pressureは定義域を持つ零pressure式 |
| 空力 | AeroSurface、Gust.none/sine、Convection.free/prescribed、Retention.all/closed_window | sineの方向zとwavevector xで直交。ρ/基準量/εは明示正 |
| 解析 | Analysis.static/transient/flutter、BranchSeed.reference/prepared | transientは1区間[0,1],h1。flutterは速度[0,1]で物理coverage要求を保持し、成功を主張しない |
| 初期・再開 | Initial.equilibrium/cold_start/prepared_history、Perturbation.none/physical、Planes.binary64/two_word、FrameStatePayload、Restart | coldは動的配列のschema形状例でありprepare後の係数数適合は別。hi/loの微小値を保持 |
| 観測 | Observation.point/stress/section/reaction/aero/mode、Basis.world/reference_material/current_material/explicit、ComparisonReference.none/file | material basisはstressのlayerをregistryに置く。成分数・順序・次元は[観測の成分形状](case-schema-fields.md#観測の成分形状) |
| 尺度・誤差群 | Scale.uniform_components/explicit_layout、ErrorGroup.required/not_applicable/covered_by | xyz尺度の明示展開、prepared layoutの全列、結合群の単一代表を区別 |
| 数値抽出 | EigenPolicy.none/full/generalized/arnoldi、EigenSelector.largest_magnitude/shift_invert、Backward.not_requested/required | Arnoldi例はn3/nev1/ncv3、seed非零。全物理modeの被覆と呼ばない |
| 保存 | CheckpointOutput.off/on | onのpathは未作成。実OS保存試験の代替ではない |
| 式 | ExprNode.constant/variable/add/subtract/multiply/divide/negate/integer_power/sin/cos/exp/log/sqrt | node単体はoperand型の文脈が必要。以下の閉じた式例でtype/定義域を照合 |

### 限定式の単位と微分を独立照合する

| 入力式とdomain | 型/値/微分の独立期待 | 正常/負対照 |
|---|---|---|
| 参照chartのE0=ε s²、ε=1e−30 | 無次元、∂sE0=2εs、∂ssE0=2ε。t/ζ微分0 | εをbinary64丸め前に保持。時刻/current x端末の使用拒否 |
| 曲面の繊維方向w=(1,0,αs) | 参照無次元vector、MAT-D05で参照法線へ射影 | wを現在qから更新する入力を拒否。法線平行はprepare拒否 |
| 指定圧力p(X,t)=p0(1+X1/L) sin(ωt+φ) | p0 pressure、X1/L無次元、ω angle/T、φ angle。p_t=p0(1+X1/L)ω cos(...)（radへ統一後） | 単位違いのadd、log(pressure)をschema拒否 |
| 現位置圧力p(x,t)=p0 exp(x1/L) | ∂x1p=p/L、∂x1x1p=p/L²、他0 | 完全なδpを点のδxへ鎖則。q配列端末/外部関数名拒否 |
| deg系、ω=180deg/s、φ=0、t=0.5s | sin(π/2)=1、時間微分0。wavevectorもdeg/L→rad/mへ | degを1radと読む変異を後続SI試験で検出 |
| exp(log(2))、sqrt(4)、整数冪2³ | 2、2、8。logは正、sqrtの必要微分域は厳密正 | denominator domainが0を含む、log≤0、sqrt<0/微分域0をprepare拒否 |
| table/PCHIP支持運動 | 各区間は同じ位置多項式からv/aを生成 | C1で二階不連続なら明示event/左右値、無eventのC2要求を拒否 |

式例は採用能力の説明であり、製品の式評価/自動微分/区間演算を実行済みではない。これらの必要微分とdomainが先に定義されて初めて後続実装へ渡せる。

### root全体の適用条件

| schema_examples | model/flow/wakeの条件 | 検査する範囲 |
|---|---|---|
| aerodynamic_steady | surface_graphと固定参照attachment、構造8集合のうち物性/部材は空。Flow known/Wake active | 架空のシェル厚さ/剛性/質量/nを要求しない |
| aerodynamic_prescribed_transient | 同じ薄い面にt³の指定並進、t0のv/a零から平衡初期化 | 位置から同じ速度/加速度、構造支持速度の二重加算無し |
| structural_static | シェル/材料/全位置固定trace、aero_surfaces空。Flow/Wake/WakeCapacity none | 架空ρ/U/コア/後流行数を要求しない |
| coupled_transient | シェル物質写像へ付いた揚力面、既知流れ、整合平衡初期化 | 同じq/vを面位置/面速度/転置力へ。解析時刻[0,1] |
| coupled_flutter | 同じ枝/後流/step、Arnoldiと両physical coverage要求がtrue | falseによる実務受入の迂回を拒否。実nを得たprepareでseed/nev/ncvと全容量を再検査 |

これらは字句と相互条件を具体化する例であり、合成設定の物理精度・容量を保証しない。I1〜I5および公開R1の受入は[製品要求](../../requirements/product-requirements.md)の全証拠を必要とする。

## DC-02 全系列に適用する負例

各変異は正常fragmentの**1条件だけ**を変更する。同じ型の全列挙成員を試す。採用仕様の最初の1例だけで系列全体を検査したとしない。

| 変異 | schemaの期待 | 後続の境界 |
|---|---|---|
| root必須9keyを各1つ削除、任意3keyへnull、未知root key | InvalidInput | 全schemaを通らず終了 |
| model全8集合を各削除/非array化、各nested objectへ未知key | InvalidInput | 空集合は削除と区別 |
| 全variantのkindを未知/別variantへ変更、他variant専用key混入 | Unsupported/InvalidInputを理由付きで区別 | 既知kindだがkey違反はInvalidInput |
| Number→bool/string/null、Id→native integer/decimal、Count→1.0/負 | InvalidInput | 暗黙変換しない |
| uint64最大Id=18446744073709551615、最大+1 | 最大は保持、超過はInvalidInput | 浮動小数点を経由しない |
| index=INT64_MAX、INT64_MAX+1、負、範囲外参照 | 字句範囲/対象arrayの条件で判定 | checked arithmeticを先行 |
| 各Vec/Matのrow/columnを±1、row-majorをflat arrayへ | InvalidInput | 転置/平坦化を推測しない |
| 同じnamespaceの重複Id、別namespace同数字 | 前者拒否、後者許可 | SourceId同一宣言の再参照は許可 |
| dangling参照、型の違う参照、禁止参照のcycle | InvalidInput | 接続の物理閉路とは区別 |
| E/密度/全厚/長さ/尺度≤0、ν範囲外、layer gap/overlap | exact値条件はschema拒否 | 材料SPD/幾何detJ/拘束rankはprepare |
| PCHIP時刻重複/逆順/1点、period≤0、event逆順 | InvalidInput | 実係数/周期の微分整合はprepare |
| transient initial+restart両方/どちらも無し | InvalidInput | root任意keyだからと物理条件を省略しない |
| node operand未来/自己/範囲外、未使用node、未知op、型/次元不一致 | InvalidInput/Unsupported | callback名をTextのまま通さない |
| point.position予算長1、stress.first_pk予算長6、reaction予算長3 | InvalidInput | FD-09の3/9/6成分を全て要求 |
| uniform尺度のxyz長1、mixed constraintをuniform指定、explicit_layout短列/別digest | InvalidInput/準備版不一致 | uniformの意味付き展開と暗黙不足埋めを区別 |
| stress.greenとcontribution=elastic/viscousの組合せ | InvalidInput | Greenひずみを応力の加法分解と混同しない |
| 結合評価のcovered_by循環/連鎖/代表不在、covers重複/未被覆 | InvalidInput | 物理的な適用と予算所有を分離 |
| 全quota同値/1不足、巨大shape積overflow | 同値許可/不足失敗、入力不変 | 推定backend capとPlanを区別 |
| 1e−400、1+10^−30の元decimal、TwoWord low=10^−30 | 元値を保持 | SI/数値consumerの前に0や1へ落とさない |
| source pathにNUL、hash63/65桁/非hex、same SourceId別定義 | InvalidInput | 不存在/改変/脱出/symlinkはsnapshot処理 |
| formula評価に外部コード/任意q/乱数/loop/条件分岐 | Unsupported | 許可演算のDAGだけ |

NaN/Infinityは元文書が構文時点で拒否する。引用文字列はNumberへ変換しない。有限な元decimalのSI変換不能・丸め包絡失敗は、schemaのUnknownKindと混同しない。

## DC-03 採用仕様の漏れを調べる

| 必須内容 | 契約箇所 | 設計の判定対象 |
|---|---|---|
| IO4 root9+3/model8 | FD-01、CS-02/04 | 全key/空/排他/型付き参照 |
| SH graph x0/w/a/f、o/Q/Lref、Hζ/2、n、固定anchor、D | FD-03 | ζ±1で実厚さ±H/2、非単位width/skew、片側/全trace |
| MAT 3D全定数/有効則/E0/軸/層/粘性 | FD-02/04/12 | Mandel順、engineering変換、初期ひずみと初期配置の区別 |
| BM K6/質量/偏心/慣性/初期strain/反り | FD-04/08 | 2節点の位置/姿勢/速度、角速度angle/T、線密度と質量 |
| CN 接着/有限剛体取付/支持/粒子・剛体/質量台帳 | FD-05 | η写像と全厚対応、元拘束/転置反力、二重計上拒否 |
| SI 荷重測度/追従接線/物理減衰/回復 | FD-05/09/12 | pressure/traction/bodyの測度、軸、同じ状態の全微分 |
| AE/UI/UW 物理面/後縁/共有/核/移流/保持 | FD-06 | dual-quarter、同じ固定ε、W-P1、ID種別/放出/実履歴 |
| AE-07 known U/sine/PCHIP | FD-02/06 | 角周波数と位相rad、e·k=0、範囲外/イベント/導関数 |
| J1 branch/BDF2/初期整合/自由定常/探索 | FD-07/08 | 物理過去/算法補助の区別、各枝と全queue、未解像を正常化しない |
| N4 尺度/原残差/TwoWord/投影/固有cluster/coverage | FD-10 | backend収束と元問題資格を別に保持、未達/容量拒否 |
| IO4 provenance/観測/結果/checkpoint | FD-08/09/11 | 元bytes/hash、側/基底/測度/確定時刻、atomic保存要求 |

## DC-04 状態と引継ぎ

Issue #110は設計であり、schema実装/数値合格/物理受入の報告ではない。独立レビューが閉じた項目だけを設計検証記録へ記す。型辞書に未定義語や、必要fieldを「後で決める」とした箇所が残る場合、全schema実装開始条件を満たしたとしない。

次の実装Issueでは、辞書全variantの正常/負例を実際のtyped validatorへ接続する。独立oracleは手書きの型/shape/厳密値/参照関係であり、製品validator自身が再出力したJSONを正解にしない。SI、式評価/微分、snapshot/importer、物理prepare、checkpoint/CLI/Pythonはそれぞれ既存IM06の依存段階へ渡し、本設計の機械検査だけで完了しない。
