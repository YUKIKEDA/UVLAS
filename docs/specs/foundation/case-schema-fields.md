# uvlas.case/1 フィールド辞書

[全schema契約](case-schema-contract.md)の規範的な型辞書である。表のkeyは全許可keyであり、`?`以外は必須である。`List`の空可否は「空可」「1以上」で示す。行列はrow arrayで保持する。dの単位指数は`(length,mass,time,angle)`の順とし、`1`は無次元を表す。以下の型名はJSON上のkeyではない。

有限式入力の採用はIssue #110に記録した利用者決定Aに基づく。数値kernelの実装済みを示さない。

## FD-01 rootと単位

| 型 | key:型 | 条件 |
|---|---|---|
| Case | schema:Tag{uvlas.case/1}, units:Units, model:Model, flow:Flow, analysis:Analysis, observations:List<Observation>, numerics:Numerics, resources:Resources, output:Output, initial?:Initial, restart?:Restart, provenance?:Provenance | observationsは空可。初期/再開の相互条件はFD-08。provenance省略は補足空 |
| Units | length:Tag{m,mm,cm,km,in,ft}, mass:Tag{kg,g,t,lbm}, time:Tag{s,ms,min}, angle:Tag{rad,deg} | tは1000kg、lbmは0.45359237kg、inは0.0254m、ftは0.3048m。degはπ/180radであり有理数と偽らない。SI段階が変換誤差を認定 |
| Model | geometry:List<Geometry>, materials:List<Material>, shells:List<Shell>, beams:List<Beam>, joints:List<Joint>, supports:List<Support>, masses:List<Mass>, aero_surfaces:List<AeroSurface> | 全8配列を明示。各空可。analysisが消費する物理部分が空なら相互条件で拒否 |
| Provenance | description:Text, references:List<Text>, assumptions:List<Text> | 全リスト空可。来歴補足であり原bytes/hashを置換しない |

角度はroot angleによる物理角を明示し、無次元ひずみ/法線/回転行列へangle倍率を掛けない。圧力は`mass/(length*time²)`、密度は`mass/length³`。混合単位行列は成分の型で決める。独立した温度・電磁場・任意単位文字列を推測しない。

## FD-02 共通の物理点・場・時間系列

| 型 | key:型 | 条件 |
|---|---|---|
| Frame | origin:Vec<Q<L>,3>, rotation:Mat<Number,3,3> | 行列列が局所軸の世界成分。右手直交の認定はprepare |
| Target | kind:Tag{shell,beam,rigid}, id:対応Ref<Shell/Beam/Mass> | rigidはMass.kind=rigidだけ |
| MaterialPoint | target:Target, coordinates:Vec<Number,3>, side:Tag{interior,lower,upper} | shellは(s,t,ζ)、beamは(s,0,0)、rigidは物体座標の長さ。unused成分の零をexact検査。shellのζ/実厚さ変換は準備版を使用 |
| Region | target:Target, kind:Tag{volume,face,edge,point}, geometry:Ref<Geometry>, boundary:Index, side:Tag{interior,lower,upper} | pointはboundary=0。面/辺は幾何に定義したindex。未知index/target種別不整合拒否 |
| PolyTerm<d,k> | powers:Vec<Count,k>, coefficient:Q<d> | 座標を指定した無次元変数に正規化してから評価。同一powersの重複を拒否 |
| ReferenceField<d,n> constant | kind:Tag{constant}, value:Vec<Q<d>,n> | 成分順は使用箇所が指定 |
| ReferenceField<d,n> polynomial | kind:Tag{polynomial}, coordinates:Tag{chart,material,world_reference}, origin:Vec<Number,3>, lengths:Vec<Number,3>, components:Vec<List<PolyTerm<d,3>>,n>, error:Vec<Q<d>,n> | chartのorigin/lengthsは無次元、material/world_referenceはL。各length>0、error≥0。変数u=(座標−origin)/lengths。空termは厳密零。近似誤差を保持し、軸の直交性はprepare |
| TimeSeries<n> constant | kind:Tag{constant}, value:Vec<Number,n> | 使用側が各成分の次元を付ける |
| TimeSeries<n> polynomial | kind:Tag{polynomial}, origin:Q<T>, scale:Q<T>, coefficients:List<Vec<Number,n>>, domain:Vec<Q<T>,2> | scale>0、coefficients1以上、昇冪順。domain下端<上端、範囲外拒否 |
| TimeSeries<n> pchip | kind:Tag{pchip}, times:List<Q<T>>, values:List<Vec<Number,n>>, outside:Outside | 同じ長さ2以上、時刻厳密増加。AE-07の同じ補間係数をprepareで生成 |
| Outside reject/hold | kind:Tag{reject,hold} | 明示されたholdは端値保持 |
| Outside periodic | kind:Tag{periodic}, period:Q<T> | period>0、時刻幅と一致。値/一階微分の周期整合はprepare |
| Motion | origin:Frame, translation:TimeSeries<3>, rotation_vector:TimeSeries<3>, events:List<Q<T>> | 並進Lと回転angleは基準Frameの局所軸成分。同じSO(3)写像から位置・速度・加速度を生成。events厳密増加、空可。相対Logと区分C2の適格性はprepare |
| Axis projected | kind:Tag{projected}, direction:ReferenceField<1,3>, angle:ReferenceField<angle,1> | 法線射影と回転はMAT-D05。零/平行/退化をprepareで拒否 |
| Axis explicit | kind:Tag{explicit}, columns:ReferenceField<1,9> | row-major3×3の行列成分。SO(3)/連続性はprepare |

ReferenceFieldとTimeSeriesはFD-12の限定式variantも持つ。callback、任意コード、名前だけの外部providerを受けない。空力の状態依存荷重は既存couplingが同じmodelから作る。近似多項式を入力する場合は座標・出自・誤差を保持し、必要微分と物理予算をprepareで検査する。

MotionはPCHIPのC1だけを加速度のC2資格へ流用しない。支持運動の各table節点は明示イベントとし、左右の二階微分を同じ補間から生成して各区間で評価する。節点加速度が一意でない場合はsideを要求し、滑らかなC2運動と表示しない。C2を要求する連続区間で左右二階微分が一致しなければprepareで拒否する。gustにはAE-07のC1を要求する。

基準原点と姿勢を $c_0,Q_0$、局所並進と回転ベクトルを $d,\phi$ とすると、世界座標の原点と姿勢は $c=c_0+Q_0d$、$Q=Q_0\exp([\phi]_\times)$ である。PCHIPでは、同じ標本から導いた勾配を区間の三次Hermite係数に使う。[指定剛体運動の開始契約](../aero/prescribed-rigid-motion-start-contract.md)

PGの物理参照は $$ rhat=(x0+width*s+skew*t,t,f(s,t)) $$ 、 $$ X0=o+Lref Q(rhat+Hhat ζ nhat/2) $$ とする。thickness.valuesは物理全厚H、Hhat=H/Lref、ζ∈[−1,1]で実法線厚さ[−H/2,H/2]へ写す。アンカーは採用graph規約s=0であり自由入力ではない。

## FD-03 幾何と外部source

| 型 | key:型 | 条件 |
|---|---|---|
| Geometry pg | kind:Tag{pg}, id:Id<Geometry>, frame:Frame, reference_length:Q<L>, x0:Rational, width:Rational, skew:Rational, degree:Vec<Count,2>, coefficients:List<PgTerm>, domain:List<Vec<Rational,2>>, thickness:Thickness, n:Count, convention:Tag{pg-graph-v1}, boundaries:List<Boundary>, continuity:List<Continuity> | reference_length/width>0、総次数≤4。domain3頂点以上、反復終点なし。両方向共通n∈{1,2,4}、anchorはs=0で固定 |
| PgTerm | powers:Vec<Count,2>, coefficient:Rational | 重複powers拒否、宣言degree内。係数は無次元graph、frame/lengthで寸法化。零term可 |
| Thickness | vertices:List<Vec<Rational,2>>, triangles:List<Vec<Index,3>>, values:List<Q<L>> | verticesとvalues同数、3点以上。区分一次全厚、全values>0、index範囲内・三角形3index相異。幾何の被覆/正則性はprepare |
| Boundary | name:Text, vertices:List<Vec<Rational,2>> | 2点以上。順序付き境界。name一意 |
| Continuity | other:Ref<Geometry>, kind:Tag{smooth,fold}, own_boundary:Index, other_boundary:Index, orientation:Tag{same,reversed} | 自分自身への参照拒否。smoothは位置/法線/全厚の共通物理graphをprepareで検査 |
| Geometry surface_graph | kind:Tag{surface_graph}, id:Id<Geometry>, frame:Frame, reference_length:Q<L>, x0:Rational, width:Rational, skew:Rational, degree:Vec<Count,2>, coefficients:List<PgTerm>, domain:List<Vec<Rational,2>>, boundaries:List<Boundary>, continuity:List<Continuity> | 空力単体の薄い参照面。PGと同じ中央graph式、構造厚さ/n/基底/剛性/質量を要求しない。reference_length/width>0、総次数≤4 |
| Geometry nastran | kind:Tag{nastran}, id:Id<Geometry>, deck:FileRef<nastran-deck>, sidecar:FileRef<nastran-sidecar>, dialect:Tag{msc-2022.1}, mode:Tag{structural-extraction}, excluded_cards:List<Text> | IO4-03のallowlistだけ。未知剛性/質量/拘束cardを追加できない。sidecar単位必須 |
| Geometry gmsh | kind:Tag{gmsh}, id:Id<Geometry>, mesh:FileRef<gmsh-msh4.1-ascii>, sidecar:FileRef<gmsh-sidecar>, reference:ReferenceGeometry | 同じexternal sourceをPG DOFと解釈しない |
| ReferenceGeometry supplied | kind:Tag{supplied}, file:FileRef<reference-geometry> | 参照graph/部材/共有面と適合するsidecar |
| ReferenceGeometry fit | kind:Tag{fit}, degree:Count, distance_budget:Q<L>, normal_budget:Q<1>, conversion_observations:List<Ref<Observation>> | degree≤4、予算>0、観測1以上。有限点から未観測曲率を認定しない |

Nastran/Gmshのsidecar内容は本caseと同じ型辞書の幾何・材料・所有/対応定義を使用し、別の暗黙単位/材料補完を作らない。カード/mesh構文とsidecarのファイルschemaは各importer契約が閉じる。`SchemaCase`はFileRefの役割とcase内対応を検査し、ファイル内のcardを検査済みとは返さない。

## FD-04 材料・シェル・梁

| 型 | key:型 | 条件 |
|---|---|---|
| Material isotropic | kind:Tag{isotropic}, id:Id<Material>, young:Q<pressure>, poisson:Q<1>, density:Q<M/L³>, provenance:Provenance, applicability:Applicability | E>0、−1<ν<1/2、ρ>0 |
| Material orthotropic | kind:Tag{orthotropic3d}, id:Id<Material>, young:Vec<Q<pressure>,3>, shear:Vec<Q<pressure>,3>, poisson:Vec<Q<1>,3>, density:Q<M/L³>, provenance:Provenance, applicability:Applicability | E1/E2/E3、G23/G13/G12、ν12/ν13/ν23。正E/G/ρ、相反性compliance/SPDはprepare |
| Material effective | kind:Tag{effective3d}, id:Id<Material>, stiffness:Mat<Q<pressure>,6,6>, components:Tag{mandel,engineering}, density:Q<M/L³>, provenance:Provenance, applicability:Applicability | xx/yy/zz/yz/xz/xy順、engineeringのせん断倍率を変換。SPD/対称化記録はprepare |
| Applicability unevaluated | kind:Tag{unevaluated} | 数値適格/実材料受入を混同しない |
| Applicability declared | kind:Tag{declared}, strain_lower:Vec<Number,6>, strain_upper:Vec<Number,6>, reference:Text | 成分ごとlower≤upper。用途の受入状態は別台帳 |
| Layer | id:Id<Layer>, interval:Vec<Rational,2>, material:Ref<Material>, axis:Axis, initial_green:ReferenceField<1,6>, viscosity:Mat<Q<pressure*T>,6,6> | intervalは正規化ζ∈[−1,1]で下端<上端。層順で全厚を隙間/重複なしにexact被覆。粘性0を明示可。PSDはprepare |
| Shell | id:Id<Shell>, geometry:Ref<Geometry>, layers:List<Layer>, damping_components:Tag{mandel,engineering}, loads:List<Load> | layers1以上、loads空可。E0はGreen、Mandel順。normal-Schur/任意ABDの入力へ切替不可 |
| Beam | id:Id<Beam>, nodes:Vec<Frame,2>, coordinate:BeamCoordinate, section:Section, initial_strain:Vec<Number,6>, residual_energy:Q<energy/L>, viscosity:Mat<Number,6,6>, loads:List<Load> | 2節点SE(3)、ひずみ順Γx/Γy/Γz/κx/κy/κz。後3成分1/L。粘性混合単位はK6×T。loads空可 |
| BeamCoordinate arclength | kind:Tag{arclength}, length:Q<L> | length>0。弦長から推測しない |
| BeamCoordinate scaled | kind:Tag{scaled}, length:Q<L>, arc_scale:Q<1> | length/scale>0。BM仕様のK/m/初期twist変換 |
| Section matrix | kind:Tag{matrix}, stiffness:Mat<Number,6,6>, line_density:Q<M/L>, center:Vec<Q<L>,3>, inertia:Mat<Q<M*L>,3,3>, envelope:Vec<Q<L>,2>, warping:Tag{free}, recovery:Tag{unavailable}, source:Text | K6のΓ/κ混合単位。重心慣性/線密度を混同しない。剛性対称SPD/慣性物理条件はprepare |
| Section circle | kind:Tag{circle}, radius:Q<L>, young:Q<pressure>, poisson:Q<1>, density:Q<M/L³>, effective_shear:Vec<Q<L²>,2>, warping:Tag{free} | 全正、ν範囲。せん断面積を既定化しない |
| Section rectangle | kind:Tag{rectangle}, sides:Vec<Q<L>,2>, young:Q<pressure>, poisson:Q<1>, density:Q<M/L³>, effective_shear:Vec<Q<L²>,2>, torsion:Q<L⁴>, warping:Tag{free} | JをIy+Izで補完しない |

梁のmatrix剛性の(i,j)単位は、一般化ひずみの長さ指数をa=(0,0,0,−1,−1,−1)としてenergy/lengthからL^(a_i+a_j)を除いた単位である。`initial_strain`各成分も同じaを使う。適格性を満たすために微小慣性/剛性/粘性を追加しない。

## FD-05 接続・支持・質量・荷重

| 型 | key:型 | 条件 |
|---|---|---|
| Joint adhesive | kind:Tag{adhesive}, id:Id<Joint>, first:Region, second:Region, mapping:List<TraceMap>, split_at_attachment:Bool | 全位置traceの対応、mapping1以上。少数collocation代用の許可ではない |
| TraceMap | first_cell:Vec<Index,2>, second_cell:Vec<Index,2>, domain:List<Vec<Rational,2>>, first:Vec<List<MapTerm>,3>, second:Vec<List<MapTerm>,3>, orientation:Tag{same,reversed} | 共通η=(η1,η2)のdomainは3頂点以上の多角形。両写像は(s,t,ζ)全成分を持つ。domain被覆/物理一致はprepare |
| MapTerm | powers:Vec<Count,2>, coefficient:Rational | ξ_a(η)、ξ_b(η)の多項式係数。重複powers拒否。全domainのpullbackを使用し、点一致だけで全trace資格としない |
| Joint rigid_attachment | kind:Tag{rigid_attachment}, id:Id<Joint>, attached:Region, frame:Target, offset:Frame | 明示した有限範囲だけ剛体化 |
| Joint weld | kind:Tag{weld}, id:Id<Joint>, first:Target, second:Target, relative:Frame | SE(3)閉路整合はprepare |
| Support constraint | kind:Tag{constraint}, id:Id<Support>, target:MaterialPoint, components:List<Index>, motion:Motion | componentsは0..5、相異、1以上。回転を持たない点の回転指定拒否 |
| Support trace_constraint | kind:Tag{trace_constraint}, id:Id<Support>, region:Region, components:List<Index>, motion:Motion | componentsは0..2の並進全位置成分、相異/1以上。領域全体の位置traceを指定参照運動へ拘束し、少数pointで置換しない |
| Support elastic | kind:Tag{elastic}, id:Id<Support>, target:MaterialPoint, reference:Motion, stiffness:Mat<Number,6,6>, damping:Mat<Number,6,6>, frame:Frame | 相対並進/回転順、混合単位。対称PSD/仕事はprepare |
| Mass particle | kind:Tag{particle}, id:Id<Mass>, target:MaterialPoint, mass:Q<M>, ledger:Ledger | mass>0、回転慣性を創作しない |
| Mass rigid | kind:Tag{rigid}, id:Id<Mass>, target:MaterialPoint, mass:Q<M>, center:Vec<Q<L>,3>, inertia:Mat<Q<M*L²>,3,3>, frame:Frame, ledger:Ledger | mass>0、重心慣性。SPDと三角不等式はprepare |
| Ledger add | kind:Tag{add}, source:Text, owner:Target | 同一source/所有領域の二重計上はprepareで拒否 |
| Ledger replace | kind:Tag{replace}, source:Text, owner:Target, replaced:List<Ref<Mass>> | replaced1以上、自己/重複/循環参照拒否。分布質量置換には対象領域を別sourceと一緒に記録 |
| Load point_force | kind:Tag{point_force}, id:Id<Load>, target:MaterialPoint, frame:Tag{world,body}, series:TimeSeries<3> | 各成分force。偏心点/姿勢から全接線 |
| Load couple | kind:Tag{couple}, id:Id<Load>, target:MaterialPoint, frame:Tag{world,body}, series:TimeSeries<3> | 各成分force*L。回転のあるtargetのみ |
| Load gravity | kind:Tag{gravity}, id:Id<Load>, acceleration:TimeSeries<3> | 各成分L/T²。同じ質量台帳を消費 |
| Load traction | kind:Tag{traction}, id:Id<Load>, region:Region, area:Tag{reference,current}, frame:Tag{world,body}, spatial:ReferenceField<pressure,3>, multiplier:TimeSeries<1> | multiplier無次元。面/向きをregionに保持。current areaの接線を省かない |
| Load pressure | kind:Tag{pressure}, id:Id<Load>, region:Region, spatial:ReferenceField<pressure,1>, multiplier:TimeSeries<1> | 正圧は指定向きの圧縮側。current法線/面積をcoupling核から計算 |
| Load body | kind:Tag{body_force}, id:Id<Load>, region:Region, spatial:ReferenceField<force/L³,3>, multiplier:TimeSeries<1> | 参照体積測度 |

支持の(並進,回転)混合単位は各相対座標の次元を使ってenergyから導く。剛体frameのないtargetをbody座標の偶力へ暗黙変換しない。接触/摩擦/離脱/拘束反りは未知tagとして未対応を返す。物理接着グラフの閉路だけで拒否したり従属行を削除しない。

## FD-06 空力・既知流れ

| 型 | key:型 | 条件 |
|---|---|---|
| AeroSurface | id:Id<AeroSurface>, geometry:Ref<Geometry>, attachment:AeroAttachment, chord:List<Number>, span:List<Number>, leading:Index, trailing:Index, orientation:Tag{positive,reversed}, sharing:List<AeroSharing> | chord/span各2以上で厳密増加、後縁の物理index範囲検査。quarter/three-quarterを別入力で変更不可 |
| AeroAttachment fixed | kind:Tag{fixed_reference} | geometryの不変参照位置を使う。架空シェル不要 |
| AeroAttachment prescribed | kind:Tag{prescribed_motion}, motion:Motion | 参照geometryに同じSE(3)運動を作用。別の構造支持速度を加算しない |
| AeroAttachment structural | kind:Tag{structural}, host:Ref<Shell> | 同じ固定物質写像のq/vから位置/面速度。motionの重畳keyは禁止 |
| AeroSharing | other:Ref<AeroSurface>, own_boundary:Index, other_boundary:Index, orientation:Tag{same,reversed} | ID共有と近接を区別。incident適格性はprepare |
| Flow none | kind:Tag{none} | aero_surfacesが空の場合だけ。架空ρ/Uを必須にしない |
| Flow known | kind:Tag{known}, density:Q<M/L³>, uniform:TimeSeries<3>, gust:Gust, reference:FlowReference, events:List<Q<T>> | density>0、uniformはL/T。events厳密増加、空可 |
| Gust none | kind:Tag{none} | 意図した突風無し |
| Gust sine | kind:Tag{sine}, amplitude:TimeSeries<1>, direction:Vec<Number,3>, wavevector:Vec<Q<angle/L>,3>, angular_frequency:Q<angle/T>, phase:Q<angle>, origin:Vec<Q<L>,3>, time_origin:Q<T> | amplitude L/T、direction無次元。wavevectorは角度/長さ、angular_frequencyは角度/時刻（Hzではない）。phaseと共にradへ変換して同じsine式を評価。e·k=0/方向正規化はprepare |
| FlowReference | density:Q<M/L³>, speed:Q<L/T>, area:Q<L²>, chord:Q<L>, span:Q<L>, moment_origin:Vec<Q<L>,3> | ρ*>0、U*≥0、他正。U*=0の係数Unavailableを保持 |
| Wake none | kind:Tag{none} | aero_surfacesが空の場合だけ |
| Wake active | kind:Tag{active}, convection:Convection, retention:Retention, core_radius:Q<L>, model:Tag{wp1-stream-p0-span}, mapping:Tag{dual-quarter-v1} | core_radius>0。全束縛/全後流で共通。理想核の別経路を自動選択しない |
| Convection free | kind:Tag{free} | 全相互誘導+既知流れ |
| Convection prescribed | kind:Tag{prescribed}, velocity:TimeSeries<3> | L/T、指定移流であることを明示 |
| Retention all | kind:Tag{all} | 全物理履歴を保持、資源上限で停止 |
| Retention window | kind:Tag{closed_window}, duration:Q<T>, distance:Q<L>, error_budget:Q<1> | duration/distance>0、error>0。閉端/退役台帳と誤差認定を後続で保持 |

指定移流をfreeと同じ物理モデルとして扱わない。wakeのRは資源上限であり、窓の物理時間/距離と別項目である。

## FD-07 解析・枝・探索

| 型 | key:型 | 条件 |
|---|---|---|
| Analysis static | kind:Tag{static}, physics:Tag{aerodynamic,structural,coupled}, branch:Branch, speed:AnalysisSpeed, wake:Wake, loads:List<Load> | speed≥0、loads空可。同じ自由/指定定常後流を含む残差 |
| Analysis transient | kind:Tag{transient}, physics:Tag{aerodynamic,structural,coupled}, branch:Branch, intervals:List<TimeInterval>, wake:Wake, loads:List<Load> | intervals1以上、loads空可。隣接時刻をexact照合、イベントを跨がない |
| Analysis flutter | kind:Tag{flutter}, physics:Tag{coupled}, branches:List<Branch>, intervals:List<SpeedInterval>, step:Q<T>, wake:Wake, search:Search, loads:List<Load> | branches/速度区間各1以上、step>0。全状態Pの安定性/coverageを後続認定 |
| AnalysisSpeed none | kind:Tag{none} | physics=structuralだけ |
| AnalysisSpeed value | kind:Tag{value}, value:Q<L/T> | value≥0、aerodynamic/coupledだけ |
| Branch | id:Id<Branch>, seed:BranchSeed, continuation:Continuation | 枝IDは周波数sort番号ではない |
| BranchSeed reference | kind:Tag{reference} | 未荷重参照を推定に使い、整合平衡の成功とはしない |
| BranchSeed prepared | kind:Tag{prepared}, history:FileRef<preparation-history> | 実際の構造/後流/入力履歴を保持 |
| Continuation fixed | kind:Tag{fixed} | 速度継続を要求しない単点/過渡のseed |
| Continuation sweep | kind:Tag{sweep}, initial_increment:Q<L/T>, minimum_increment:Q<L/T>, maximum_increment:Q<L/T>, arclength:ArcLength, maximum_points:PositiveCount, maximum_retries:Count | 0<min≤initial≤max。失敗を擬似逆で隠さない |
| ArcLength | state_scales:List<Scale>, speed_scale:Q<L/T>, orientation:Tag{previous_tangent}, enable_at_limit:Bool | 全代表尺度正。列順/単位は対象layoutと照合 |
| TimeInterval | start:Q<T>, end:Q<T>, step:Q<T>, restart:Tag{initial,event,step_change,chart_change} | start<end、step>0、(end−start)/stepが整数。一般可変刻みBDF2を含めない |
| SpeedInterval | lower:Q<L/T>, upper:Q<L/T> | 0≤lower<upper、重複/隙間の扱いをsearch対象域としてexact保持 |
| Search | queue:Tag{lowest_interval_first}, speed_absolute:Q<L/T>, speed_relative:Q<1>, growth_budget:Q<1/T>, mode_association:Association, maximum_points:PositiveCount, maximum_retries:Count, coverage:Coverage | 許容差/予算≥0、速度尺度は別に正。完了点/未処理queueを後続で原子保存 |
| Association | physical_observations:List<Ref<Observation>>, mac_minimum:Q<1>, subspace_angle_maximum:Q<angle>, eigenvalue_absolute:Q<1/T>, eigenvalue_relative:Q<1> | 観測1以上。0≤MAC≤1、他≥0。近接根を周波数sortで同定しない |
| Coverage | branches:List<Ref<Branch>>, mode_observations:List<Ref<Observation>>, require_physical:Bool | 各1以上。flutterのcaseではtrue必須。代数抽出の成功でphysical資格を作らない |

physics=aerodynamicはaero_surfaces1以上、shells/beams/masses/joints/supports空を許し、架空の構造物性を要求しない。physics=structuralはaero_surfaces空、flow/wake/wake_rowsはnone。physics=coupledは構造部材と構造attachmentを持つ面が各1以上。flutterのCoverage.require_physicalとEigenRequirements.physical_coverageはtrue必須であり、falseで実務flutterの終了gateを迂回できない。内部の代数抽出APIで未認定packetを扱えることとは別である。構造単体に空力density/coreを要求しない。

`Analysis.loads`は部材をまたぐ重力/荷重を置く。Shell/BeamのloadsとLoadIdを全域で一意にし、同じ荷重を二重計上しない。staticでは初期条件の動的履歴を消費しない。flutterでは対象枝ごとに同じ支持/材料/モデル/後流政策を固定する。

## FD-08 初期状態と再開要求

| 型 | key:型 | 条件 |
|---|---|---|
| Initial equilibrium | kind:Tag{equilibrium}, branch:Ref<Branch>, perturbation:Perturbation | 同じ支持/流れ/荷重/後流の平衡から整合加速度等を作る |
| Initial cold | kind:Tag{cold_start}, state:InitialState, circulation:List<Q<L²/T>>, wake:InitialWake | 明示した冷始動。欠けた物理過去を零埋めしない |
| Initial prepared | kind:Tag{prepared_history}, file:FileRef<initial-history> | 実履歴と算法履歴の区別を内部状態validatorが検査 |
| Perturbation none | kind:Tag{none} | 摂動無し |
| Perturbation physical | kind:Tag{physical}, displacements:List<PointValue>, velocities:List<PointValue> | 各空可。同じ物理観測/位置場へ写し、係数番号を推測しない |
| PointValue | point:MaterialPoint, value:Vec<Number,3> | displacementはL、velocityはL/T。重複point拒否 |
| InitialState | shells:List<CoefficientState>, beams:List<BeamState>, rigids:List<RigidState> | 各空可、対象modelの全必要entityを一度ずつ含む。欠落は参照値で補わない |
| CoefficientState | shell:Ref<Shell>, reference:Digest, representation:Tag{absolute,reference_increment}, position:Planes, velocity:Planes | position L、velocity L/T。同じ準備基底n_s×3とのshapeはprepare後に照合 |
| Planes binary64 | kind:Tag{binary64}, high:List<Number> | decimalをまだ丸めない。返却要求とは別の入力表現 |
| Planes two_word | kind:Tag{two_word}, high:List<Number>, low:List<Number> | 同じ長さ。hi/loを先に加算しない |
| BeamState | beam:Ref<Beam>, nodes:Vec<FrameStatePayload,2> | model Beam.nodesと同じ順序で両端を保持。beam重複/欠落拒否 |
| RigidState | rigid:Ref<Mass>, state:FrameStatePayload | rigidはkind=rigidのMassのみ。重複/欠落拒否 |
| FrameStatePayload | position:Vec<Q<L>,3>, rotation:Mat<Number,3,3>, linear_velocity:Vec<Q<L/T>,3>, angular_velocity:Vec<Q<angle/T>,3> | 角速度はroot angle/timeからrad/sへ変換。回転行列は無次元でangle倍率を掛けない。SO(3)とチャートはprepare |
| InitialWake | surfaces:List<SurfaceWake> | 各面を一度ずつ。冷始動の空wakeは行数0として明示 |
| SurfaceWake | surface:Ref<AeroSurface>, rows:Count, points:List<Vec<Q<L>,3>>, coefficients:List<Q<L²/T>>, row_ids:List<Id<WakeRow>>, point_ids:List<Id<WakePoint>>, band_ids:List<Id<WakeBand>>, birth_times:List<Q<T>> | pointsはrows*(N+1)、coefficientsはrows*N、row_ids/birth_timesはrows、point_idsはpoints長、band_idsはN。N=span.size−1。全積をchecked照合。行/点/帯のIDは別namespace。重複物質ID/未来のbirth拒否。付着IDと放出IDは別 |
| Restart | directory:Text, manifest_sha256:Digest, compatibility:Tag{same_environment_bitwise,qualified_cross_environment}, budget:List<ObservationBudget> | directory非空/NUL不可。same_environmentはbudget空、crossは1以上。未知の互換性を正常にしない |

transientには `initial` と `restart` のどちらか一方が必須。static/flutterには `initial` を置かず、restartは探索/平衡保存状態の再開に限る。static/flutterのseedはBranchに持つ。restartは既存checkpoint manifestとchunkのschema/type/shape/bytes/hash/Model/参照/精度/時刻/ID/履歴を検査する要求であり、caseに不透明な任意state objectを埋め込まない。IO4-06の実状態形式は後続checkpoint契約が閉じる。schema段階はディレクトリを開かない。

## FD-09 観測・誤差配分

| 型 | key:型 | 条件 |
|---|---|---|
| Observation point | kind:Tag{point}, id:Id<Observation>, point:MaterialPoint, quantity:Tag{position,displacement,velocity,acceleration,force}, basis:Basis, comparison:Comparison, budget:ObservationBudget | quantityと基底/単位を一致させる |
| Observation stress | kind:Tag{stress}, id:Id<Observation>, point:MaterialPoint, layer:Ref<Layer>, interface:Tag{interior,lower,upper,both}, measure:Tag{green,second_pk,first_pk,cauchy}, contribution:Tag{elastic,viscous,total}, basis:Basis, comparison:Comparison, budget:ObservationBudget | layerのownerとpoint.shell一致。measure=greenはcontribution=totalだけを許可し、elastic/viscousを拒否。界面両側を平均しない |
| Observation section | kind:Tag{section}, id:Id<Observation>, region:Region, origin:Vec<Q<L>,3>, basis:Basis, quantity:Tag{force,moment}, comparison:Comparison, budget:ObservationBudget | 物理切断面の向きと基準点を保持 |
| Observation reaction | kind:Tag{reaction}, id:Id<Observation>, support:Ref<Support>, origin:Vec<Q<L>,3>, basis:Basis, comparison:Comparison, budget:ObservationBudget | 物理共役量。従属乗数の一意性を主張しない |
| Observation aero | kind:Tag{aero}, id:Id<Observation>, surfaces:List<Ref<AeroSurface>>, origin:Vec<Q<L>,3>, basis:Basis, quantity:Tag{force,moment,coefficient}, comparison:Comparison, budget:ObservationBudget | surfaces1以上/重複拒否。係数にFlowReferenceを使用 |
| Observation mode | kind:Tag{mode}, id:Id<Observation>, targets:List<Target>, points:List<MaterialPoint>, basis:Basis, quantity:Tag{frequency,growth,subspace}, comparison:Comparison, budget:ObservationBudget | targets/points各1以上。固有vector係数番号を物理IDにしない |
| Basis world/material | kind:Tag{world,reference_material,current_material} | materialは所有layerがある観測だけ |
| Basis explicit | kind:Tag{explicit}, rotation:Mat<Number,3,3> | SO(3)認定はprepare |
| Comparison | window:Vec<Q<T>,2>, statistic:Tag{instantaneous,history,peak,rms}, side:Tag{left,right,both}, reference:ComparisonReference | 下端≤上端。instantaneousは同値、他は正幅。不連続の左右を保つ |
| ComparisonReference none | kind:Tag{none} | 比較参照未設定、観測を受入済みとしない |
| ComparisonReference file | kind:Tag{file}, source:FileRef<observation-reference>, columns:List<Text>, uncertainty:List<Number> | columns1以上、uncertainty同数かつ≥0。基底/単位/条件の照合は受入段階 |
| ObservationBudget | absolute:List<Number>, relative:List<Number>, scale:List<Number>, groups:List<ErrorGroup> | 3配列は直後の「観測の成分形状」に定めるnと完全一致。absolute/relative≥0、scale>0。absolute/scaleの次元は各観測成分、relativeは無次元。固定製品予算より緩くしない |
| ErrorGroup not_applicable | status:Tag{not_applicable}, kind:Tag{geometry,structure_grid,aero_grid,time_wake_age,wake_length,core,quadrature,algebra,newton,eigen_range,search_width,load_mapping}, reason:Text | reason非空。適用しない物理/解析条件を明示、必須群を無効化するスイッチではない |
| ErrorGroup covered_by | status:Tag{covered_by}, kind:Tag{geometry,structure_grid,aero_grid,time_wake_age,wake_length,core,quadrature,algebra,newton,eigen_range,search_width,load_mapping}, representative:Tag{geometry,structure_grid,aero_grid,time_wake_age,wake_length,core,quadrature,algebra,newton,eigen_range,search_width,load_mapping} | この群は物理的に適用し、直接指定したrequired代表群が結合評価/予算計上を所有する。自己参照/covered_byの連鎖を拒否 |
| ErrorGroup required | status:Tag{required}, kind:Tag{geometry,structure_grid,aero_grid,time_wake_age,wake_length,core,quadrature,algebra,newton,eigen_range,search_width,load_mapping}, covers:List<Tag{同じ集合}>, absolute:List<Number>, relative:List<Number> | covers1以上かつ自分を含む。budgetの成分長一致/非負。二重計上・別caseの値を流用しない |

### 観測の成分形状

表のnは`ObservationBudget.absolute/relative/scale`の各長さである。`ErrorGroup.required`のabsolute/relativeも同じnとする。少数成分を全成分へ繰り返さない。点の単一成分だけが必要な場合もn成分の出力・予算を持ち、明示していない成分を検査済みと表示しない。

| Observationとquantity/measure | nと順序 | absolute/scaleの次元と定義 |
|---|---|---|
| point.position / displacement | 3、x/y/z | L。positionは現在位置、displacementは同じ参照点からの差 |
| point.velocity / acceleration | 3、x/y/z | L/T、L/T²。物理速度/物理加速度であり算法補助値ではない |
| point.force、section.force、aero.force | 3、x/y/z | force=M L/T²。指定basisの成分 |
| section.moment、aero.moment | 3、x/y/z | force×L。指定origin周りの右手正の偶力 |
| aero.coefficient | 3、x/y/z | 無次元の力係数、各力成分をFlowReferenceの同じ動圧×areaで除す。基準動圧0はUnavailable |
| stress.green | 6、xx/yy/zz/yz/xz/xy | 無次元の対称Greenひずみテンソル成分。せん断はEyz等でありengineeringの2倍やMandelの√2倍を含まない |
| stress.second_pk / cauchy | 6、xx/yy/zz/yz/xz/xy | pressure=M/(L T²)の対称応力テンソル成分。せん断成分に√2を掛けない |
| stress.first_pk | 9、xx/xy/xz/yx/yy/yz/zx/zy/zz | pressure。行は現在の力成分、列は参照面法線成分。対称と仮定して6成分へ落とさない |
| reaction | 6、Fx/Fy/Fz/Mx/My/Mz | 最初の3成分force、後3成分force×L。指定origin/basisに戻した物理反力・偶力であり拘束乗数の生配列ではない |
| mode.frequency | 1、frequency | 1/Tの周期周波数Hz。step mapでは同じモードの位相分岐を確定したlog(λ)/hの虚部の絶対値を2πで除す。角周波数をHzと混同しない |
| mode.growth | 1、growth | 1/T。step mapの位相分岐を固定したlog(λ)/hの実部 |
| mode.subspace | 1、maximum_principal_angle | angle。同じpoints順・targets・正の物理尺度で作った二つの対応部分空間の最大主角。単一vectorの位相差や係数の差ではない |

Greenひずみは配置から定まる単一テンソルであり、弾性応力/粘性応力の和へ分解しない。green/totalはその全Greenひずみを意味し、E−E0を無断で返さない。弾性/粘性/合計のcontributionはsecond_pk/first_pk/cauchyの同じ応力測度へ適用する。

`basis`は表のx/y/z軸を定める。応力・ひずみは同じmeasureのテンソルを指定基底へ変換してからこの順に並べる。第一PKの二つの添字には指定basisの現在側/参照側への同じ物理対応を保持する。材料軸の初期Green場の入力順/Mandel変換とは別の出力境界である。

interface=bothではlower/upperを別packetとして返し、**それぞれ**同じn成分予算を適用する。nを2倍にした連結値や平均値へ置換しない。subspaceの比較参照は同じ物理観測・次元・尺度を持つ必要がある。`ComparisonReference.none`では比較角の達成を認定できず、subspaceの比較値はUnavailableとする。一般化問題の正のλはsqrt(λ)/(2π)とし、負/零λを勝手に正へ射影しない。step mapのalias/対数枝はJ1の刻み細分化・連続mode追跡と照合する。固有振動数を出せない負λ/零mode、位相分岐未確定、対応未解像も理由付き未認定状態とし、0で埋めない。

各ObservationBudget.groupsとNumerics.error_groupsは全12群を各1回宣言する。物理的に適用する群はrequiredまたはcovered_by、適用しない群だけnot_applicableと理由を持つ。required代表群のcoversは自分自身と、その代表を直接参照するcovered_by群だけを重複なく列挙する。required/covered_byの物理適用と、requiredだけが持つ予算計上の所有を分ける。全適用群がちょうど一つのrequired代表へ対応することを検査し、自己参照/連鎖/循環/二重代表/未被覆/不一致coversを拒否する。代表をnot_applicableにしてはならない。各required群の絶対/相対配分だけを一回ずつ合計し、親budget以下とする。実際の誤差未評価をnot_applicableと呼ばない。

`Comparison.window`は静的caseでは同じ評価時刻の単一点とする。modeの物理周波数/growthは1/T、位相/部分空間角は無次元/angleの区別を保持する。切断面/応力measureが未対応ならゼロ出力へ置換しない。絶対/相対判定は固定の正尺度を使い、未知のゼロ近傍尺度を推測しない。

## FD-10 数値政策と求解上限

| 型 | key:型 | 条件 |
|---|---|---|
| Numerics | patch:PatchPolicy, quadrature:Quadrature, linear:LinearPolicy, newton:NewtonPolicy, precision:PrecisionPolicy, eigen:EigenPolicy, wake_rows:WakeCapacity, scales:List<Scale>, error_groups:List<ErrorGroup> | scales/error_groups各1以上。Rはwake_rowsで資源停止し、古いwakeを黙って削除しない |
| WakeCapacity none | kind:Tag{none} | aero_surfaces空だけ |
| WakeCapacity limited | kind:Tag{limited}, maximum_rows:PositiveCount | 全/窓とも実使用上限。到達時停止 |
| PatchPolicy none | kind:Tag{none} | shells空だけ。aeroのみの幾何をPG構造基底にしない |
| PatchPolicy active | kind:Tag{active}, maximum_degree:Count, maximum_partition:PositiveCount, geometry_error:Q<L>, trace_error:Q<L>, rank:Tag{exact_no_drop} | degree≤4、partitionは採用範囲、誤差>0。小方向削除の許可を置かない |
| Quadrature | structure:StructuralQuadrature, aero:AeroQuadrature | structure=noneはshells/beams空、aero=noneはaero_surfaces空の場合のみ |
| StructuralQuadrature none | kind:Tag{none} | 構造積分を消費しない |
| StructuralQuadrature active | kind:Tag{active}, surface_degree:PositiveCount, thickness_points:PositiveCount, independent_surface_degree:PositiveCount, independent_thickness_points:PositiveCount, absolute:Q<1>, relative:Q<1>, maximum_subdivisions:Count | 独立次数を主次数より大きくする。層/領域/背景cell分割を省かない。予算≥0 |
| AeroQuadrature none | kind:Tag{none} | 空力積分を消費しない |
| AeroQuadrature active | kind:Tag{wp1_induction_v1}, absolute:Q<1>, relative:Q<1> | 予算≥0。近接解析式/遠方12対8/距離4の採用評価法を保つ |
| Scale uniform_components | kind:Tag{uniform_components}, field:Tag{shell_position,beam_translation,beam_rotation,rigid_translation,rigid_rotation,structural_velocity,structural_angular_velocity,circulation,wake_position,wake_strength,constraint,force_balance,moment_balance,aero_residual,transport_residual,newton_correction,arclength_state}, targets:List<Target>, values:List<Number> | 下表で許可する物理成分の正尺度。対象の全対応DOFへ明示政策として展開し、係数数を利用者に推測させない |
| Scale physical_layout_policy | kind:Tag{physical_layout_policy}, field:Tag{mixed_blocks}, constraint:ConstraintScalePolicy | 全caseで高々1件。拘束/乗数/Newton/弧長を物理primary尺度から構成する初回政策。実layoutの推測を要求しない |
| ConstraintScalePolicy none | kind:Tag{none} | jointsとconstraint/trace_constraint支持が空、内部constraint/乗数block無し。弾性/減衰支持だけならnone可 |
| ConstraintScalePolicy physical | kind:Tag{physical}, position:ComponentScale3, orientation:ComponentScale3, multiplier_work_scale:Q<energy>, reaction_force:Vec<Q<force>,3>, reaction_moment:Vec<Q<force*L>,3> | 必要なposition/orientationを下の条件で明示、work/force/moment尺度は全正。元反力を内部乗数の尺度で代用しない |
| ComponentScale3 none | kind:Tag{none} | 対応する物理拘束成分を持たないときだけ |
| ComponentScale3 values | kind:Tag{values}, values:Vec<Number,3> | 全正。positionはL、orientationはangle。xyzまたは右局所回転xyz順 |
| Scale explicit_layout | kind:Tag{explicit_layout}, field:Tag{shell_position,beam_translation,beam_rotation,rigid_translation,rigid_rotation,structural_velocity,structural_angular_velocity,circulation,wake_position,wake_strength,constraint,force_balance,moment_balance,aero_residual,transport_residual,newton_correction,arclength_state}, targets:List<Target>, layout:Digest, dimensions:List<Vec<signed-native-integer,4>>, values:List<Number> | 同じ準備済layoutの全列順。両配列は同長、values全正、dimension指数int16範囲。digest/内容/列数の不一致を拒否 |
| LinearPolicy direct | backend:Tag{dense_lu,sparse_lu,sparse_llt}, accuracy:LinearAccuracy, maximum_refinements:Count, backend_memory:BackendMemory | refinemement≤5。LLTのSPD資格、元残差/制約/反力を省かない |
| LinearPolicy gmres | backend:Tag{gmres}, accuracy:LinearAccuracy, maximum_iterations:PositiveCount, restart:PositiveCount, maximum_applications:Count, maximum_adjoint:Count, maximum_preconditioner:Count, maximum_true_residual_checks:Count, backend_memory:BackendMemory | restartとmaxiterの実n条件はprepare。認定用作用も課金 |
| LinearAccuracy | absolute:Q<1>, relative:Q<1>, backward:Q<1>, norm:Tag{infinity,l2} | 各予算≥0。元残差と成分別後退誤差を別判定 |
| NewtonPolicy | maximum_iterations:PositiveCount, maximum_line_searches:PositiveCount, absolute:Q<1>, relative:Q<1>, correction:Q<1>, minimum_step:Q<1> | 予算≥0、0<minimum_step≤1。元残差の未達を補正normだけで受理しない |
| PrecisionPolicy | input:Tag{preserve}, retained:Tag{adaptive_binary64_two_word}, returned:Tag{retained,binary64_projection}, projection_absolute:Q<1>, projection_relative:Q<1> | 予算≥0。明示投影は実投影値を再資格。低位破棄不可 |
| BackendMemory | admission:Tag{estimate,require_proven_bound}, maximum_bytes:Count | 同じbackendの実際のadmissionへ渡す。推定をmalloc/RSSの上限と呼ばない |
| EigenPolicy none | kind:Tag{none} | flutterでは拒否 |
| EigenPolicy full | kind:Tag{full}, target:Tag{step_map}, eigenvalue_scale:Number, packet:Tag{spectrum,pairs,subspace}, requirements:EigenRequirements, limits:EigenLimits | 完全Schur/pair/要求clusterを既存製品入口へ |
| EigenPolicy generalized | kind:Tag{generalized}, target:Tag{structural_modes}, eigenvalue_scale:Q<1/T²>, requirements:EigenRequirements, limits:EigenLimits | K/M対称性とM正定値を原系で認定、負λを削除しない |
| EigenPolicy arnoldi | kind:Tag{arnoldi}, target:Tag{step_map}, eigenvalue_scale:Number, packet:Tag{pairs,subspace}, selector:EigenSelector, nev:PositiveCount, ncv:PositiveCount, seed:ArnoldiSeed, backend_tolerance:Q<1>, requirements:EigenRequirements, limits:EigenLimits | eigenvalue_scale>0、backend_tolerance>0。1≤nev≤n−2、nev+2≤ncv≤nは実n判明後prepareで検査。24/64を小nへclipしない。seedは明示非零 |
| Complex | real:Number, imaginary:Number | 実入力でもJSON boolean/arrayの再解釈なし |
| ArnoldiSeed uniform | kind:Tag{uniform_nonzero}, value:Complex | value≠0。実n全列へ同じ値を明示展開。nをcaseで推測しない |
| ArnoldiSeed explicit | kind:Tag{explicit_layout}, layout:Digest, values:List<Complex> | 同じ準備済layoutの全n列、少なくとも1値非零。全列/内容を照合 |
| EigenSelector largest | kind:Tag{largest_magnitude} | 元作用の固有値絶対値 |
| EigenSelector shift | kind:Tag{shift_invert}, shift:Complex | step mapの無次元λ。sλも正の無次元。shifted残差だけで認定しない |
| EigenRequirements | residual:Q<1>, residual_observation:Ref<Observation>, normalization:Q<1>, normalization_observation:Ref<Observation>, orthogonality:Q<1>, orthogonality_observation:Ref<Observation>, backward:Backward, association_absolute:Q<1>, association_relative:Q<1>, left:Bool, clusters:List<EigenCluster>, physical_coverage:Bool | 各予算≥0。Gram<1、実A共役対、左右別空間/associationは採用契約。clusters空可。flutterではphysical_coverage=true必須 |
| Backward none | kind:Tag{not_requested} | 成分別βの省略を原残差省略と混同しない |
| Backward required | kind:Tag{required}, budget:Q<1>, observation:Ref<Observation> | budget≥0。微小tailを閾値で消さない |
| EigenCluster | members:List<Index>, expected_dimension:PositiveCount | members相異/1以上、dimension一致。完全packet順とのassociationは後続抽出契約で照合。全Qを一部clusterの資格へ転用しない |
| EigenLimits | maximum_order:Count, maximum_results:Count, maximum_ncv:Count, maximum_backend_iterations:Count, maximum_attempts:Count, maximum_forward:Count, maximum_adjoint:Count, maximum_inverse:Count, backend_memory:BackendMemory | 全呼出し課金、枠超過前に止める。n/nev+1/複素対scratchのchecked計上 |

### 尺度の明示的な展開

`uniform_components`のvaluesは下表の成分数と完全一致する。同じfield/targetへ重複する尺度宣言を拒否する。targetsが空ならそのfieldの全対象、空でなければ列挙した対象の全対応成分を表す。対象種別不一致を拒否する。`explicit_layout`のtargetsにも同じ規則を使うが、列順と次元は準備済layoutの全列を指定する。

| uniformを許すfield | values順・長さ | 物理次元 |
|---|---|---|
| shell_position、beam_translation、rigid_translation | x/y/z、3 | L |
| beam_rotation、rigid_rotation | 右局所回転x/y/z、3 | angle。SI変換後rad |
| structural_velocity | 物理並進速度x/y/z、3 | L/T。角速度を同じ尺度で埋めない |
| structural_angular_velocity | 物理角速度x/y/z、3 | angle/T |
| circulation、wake_strength | scalar、1 | L²/T |
| wake_position | x/y/z、3 | L |
| force_balance | x/y/z、3 | force |
| moment_balance | x/y/z、3 | force×L |
| aero_residual | scalar、1 | L/Tの非貫通残差 |
| transport_residual | x/y/z、3 | L/Tの物質輸送残差 |

constraint/newton_correction/arclength_stateは準備済み混合layoutに依存するため`uniform_components`を拒否する。初回caseは`physical_layout_policy`を用いて下の有限規則から構成し、準備後の全列を指定する場合だけ`explicit_layout`を用いる。初回からlayout digest/係数数を推測させない。両経路を同じfieldに重複指定してはならない。

### 初回mixed layoutを物理尺度から閉じる

`physical_layout_policy`は以下の必須primary宣言と組み合わせる。各primaryの値は利用者が計算前に指定した代表物理量であり、結果・現在の解・ゼロ右辺から推定しない。

| 適用対象 | 必須のuniform_componentsまたは同じprepared explicit_layout |
|---|---|
| 各shell | shell_position、force_balance |
| 各beam | beam_translation、beam_rotation、force_balance、moment_balance |
| 各rigid Mass | rigid_translation、rigid_rotation、force_balance、moment_balance |
| 構造を含むtransient/flutter | structural_velocity、beam/rigidがあればstructural_angular_velocity |
| active aero | circulation、wake_position、wake_strength、aero_residual、transport_residual |
| jointsまたはconstraint/trace_constraint支持あり | ConstraintScalePolicy physical。並進位置拘束があればposition=values、姿勢拘束があればorientation=values。該当しない方だけnone |

拘束の原成分を $$ C_i $$ 、対応するcallerの正物理尺度を $$ c_i $$ とする。元componentの物理kind/axis/targetと $$ c_i $$ をprepareが保持し、**先に** $$ \widehat C_i=C_i/c_i $$ と無次元化する。その後の同次元の独立化行変換を $$ T $$ とすると、変換後は $$ C'_k=\sum_i T_{ki}\widehat C_i $$ 、正の行尺度は $$ s'_k=\sum_i|T_{ki}| $$ の上向き包絡とする。非零独立行だけを資格付きで持つため $$ s'_k>0 $$ が必要であり、表現不能/零行を小値cutで修復しない。mとradの値を直接加えない。変換誤差も元量の単位へ戻して伝播する。

無次元化拘束へ共役な内部乗数は仕事次元を持つ。内部乗数を使うKKT/補正/反力回復経路には明示した正の`multiplier_work_scale`を適用する。forward行変換の列包絡から逆dual尺度を推定しない。既存CNの双対写像を同じ元/変換系へ保存し、rank削減時には非一意な元乗数の個別尺度保証を生成しない。原拘束は $$ c_i $$ 、元力学残差はprimaryのforce/moment、物理反力はreaction_force/reaction_momentで、それぞれ独立に資格を確認する。変換後の小残差で代替しない。

各Newton未知量/原残差slotはprepareで物理由来（primary field、target、axis、必要なwork block）と同じ単位を持つ。Newtonのcorrection尺度はこの明示尺度を同じ実slot順へ連結し、状態/拘束変換には上記の無次元化と同じ対応を使う。未分類slotや不足primaryを内部既定1で埋めない。全slotが対応したことをsolve前に検査する。枝の弧長は同じ状態尺度 $$ s_{z,i} $$ とContinuationの正speed_scale $$ s_U $$ から $$ W_z=\operatorname{diag}(1/s_{z,i}^2),\ w_U=1/s_U^2 $$ と固定する。数値評価は先に尺度化し、二乗のoverflow/underflowを避ける。ArcLength.state_scalesの同じfield/targetはNumerics側の同じ宣言内容と一致させ、別の動的列数を要求しない。

Arnoldiのuniform seedも同じく、実nをprepare後に知ってから全n列を事前確保して展開し、実vector・列順・generation・内容を記録する。明示したnev/ncvを小nへclipせず、実nで成立しなければ失敗する。seed展開や上の尺度展開が物理mode被覆の資格を与えることはない。

prepareはuniform宣言から各対象の実列への展開表を作り、元の宣言、単位、対象、実列順、正尺度の全配列と内容identityを所有する。同じxyz値を各Cartesian係数群へ適用することは**指定した政策の展開**であり、短いexplicit_layout配列の不足埋めではない。基底/座標変換がある場合は同じ物理成分の意味を保持する対応を準備版に含め、宣言を無検査で別layoutへ流用しない。実solve前にはN4の全成分尺度/単位・内容同一性を照合する。未被覆の必要成分には既定1を補わず、欠落として診断する。

数値政策の配置は内部のencodingであり、ケースが既存の受入精度を緩める許可ではない。`Scale.field`は上の列挙だけを受理し、SI/prepareで実layout・寸法・成分数を全照合する。値を勝手にbroadcastしない。`EigenCluster.members`は明示した固定原問題の完全packetへの要求であり、風速を変えた物理mode追跡IDとして再利用しない。

`step_map`の固有値とshift、sλは無次元であり、物理growth/frequencyは採用J1に従いlog(λ)/hから得る。`structural_modes`のλとsλは1/T²、周波数は負λや零を勝手に切らず別状態で解釈する。association予算はλ/sλに対する無次元値である。実行時にFullSelectionを作るindexは、同じ完全packetのsorted physical eigenvalues順へ結び付ける。別風速/別packetへindexだけを流用しない。

## FD-11 資源・出力・保存

| 型 | key:型 | 条件 |
|---|---|---|
| Resources | profile:Tag{24gib,48gib,custom}, threads:PositiveCount, maximum_memory_bytes:Count, maximum_wall_seconds:Q<T>, maximum_steps:Count, maximum_total_iterations:Count, checkpoint_temporary_bytes:Count, input:InputQuotas | wall>0。profile名だけで実RSS適合と認定しない |
| InputQuotas | maximum_source_bytes:Count, maximum_depth:Count, maximum_nodes:Count, maximum_members:Count, maximum_decoded_string_bytes:Count, maximum_decimal_token_bytes:Count, maximum_total_text_bytes:Count, maximum_typed_bytes:Count, maximum_entities:Count, maximum_references:Count, maximum_scalars:Count, maximum_polynomial_terms:Count, maximum_table_samples:Count, maximum_expression_nodes:Count, maximum_expression_edges:Count, maximum_expression_depth:Count, maximum_expression_power:Count, maximum_expression_points:Count, maximum_diagnostic_path_bytes:Count | 原文書入口のlimitsと同じ値で再照合。文書に書いたquotaをparse前のホスト上限の代用にしない。実上限は両者のmin |
| Output | directory:Text, observations:OutputCadence, visualization:Visualization, checkpoint:CheckpointOutput, overwrite:Tag{reject_existing}, resume:Tag{explicit_only} | directory非空/NUL不可。勝手な上書き/暗黙resumeなし |
| OutputCadence | every_steps:PositiveCount, observations:List<Ref<Observation>> | 観測空可、重複拒否。確定境界で出力 |
| Visualization | enabled:Bool, every_steps:PositiveCount, layers:List<Ref<Layer>>, interfaces:Tag{both}, reference_geometry:Bool, current_geometry:Bool | layers空可。未認定応力を0で埋めない。標本格子とPG DOFを区別 |
| CheckpointOutput off | kind:Tag{off} | 保存しないことを明示 |
| CheckpointOutput on | kind:Tag{on}, every_steps:PositiveCount, directory:Text, retain:PositiveCount | unique temporary→flush/close/hash→manifest last→同FSrename。確定dir上書き不可。retainで現確定や利用中親を削除しない |

ホストが先に与えたInputLimits/SchemaLimitsは、入力文書自身の資源宣言によって拡張しない。24/48GiBは全caseの予算・外部backend・実測認定の対象であり、本schemaの字句検査を通っただけではPERF受入にならない。保存先の安全な作成/既存検査/renameとOS差は実行・checkpoint実装の責任である。

## FD-12 限定式入力

利用者決定A（Issue #110）により、算術・整数冪・sin/cos・exp/log・sqrtだけを持つ式DAGを追加する。定数、多項式、表、AE-07の解析gustは維持する。これは外部コードの実行許可ではない。

| 型 | key:型 | 条件 |
|---|---|---|
| ReferenceField<d,n> expression | kind:Tag{expression}, expression:Expression<n>, error:Vec<Q<d>,n> | domainはreferenceだけ、出力次元全てd。error≥0。材料値/軸/E0を試行状態から再生成しない |
| TimeSeries<n> expression | kind:Tag{expression}, expression:Expression<n> | domainはtimeだけ。出力次元は使用側の成分型 |
| Load current_expression | kind:Tag{current_pressure}, id:Id<Load>, region:Region, expression:Expression<1> | domainはcurrent_position_time、出力pressure。正圧の符号はpressureと同じ。x微分を点写像へ鎖則で引き戻す |
| Expression<n> | domain:ExpressionDomain, nodes:List<ExprNode>, outputs:Vec<Index,n> | nodes1以上、index=配列位置。各operand indexは自分より小さく、全出力index範囲内。未使用node拒否。循環/任意名前参照なし |
| ExpressionDomain reference | kind:Tag{reference}, coordinates:Tag{chart,material,world_reference}, lower:Vec<Number,3>, upper:Vec<Number,3> | chartは(s,t,ζ)無次元、他はL。lower≤upper、物理領域の被覆はprepare |
| ExpressionDomain time | kind:Tag{time}, lower:Q<T>, upper:Q<T>, events:List<Q<T>> | lower<upper、events厳密増加/範囲内。各イベント左右を別区間に分ける |
| ExpressionDomain current_position_time | kind:Tag{current_position_time}, position_lower:Vec<Q<L>,3>, position_upper:Vec<Q<L>,3>, time_lower:Q<T>, time_upper:Q<T>, events:List<Q<T>> | 各lower<upper。任意q配列/v/solver内部状態を端末にしない |
| ExprNode constant | op:Tag{constant}, value:Number, dimensions:Vec<signed-native-integer,4> | SI前のroot単位で解釈。dimension各指数はint16範囲、演算による和もchecked |
| ExprNode variable | op:Tag{variable}, axis:Index | reference0..2、time0、current_position_time0..3（x/y/z/t）。単位はdomainで定義した固定座標の物理単位 |
| ExprNode binary | op:Tag{add,subtract,multiply,divide}, first:Index, second:Index | add/subtractは同次元、multiply/divideは指数和/差 |
| ExprNode negate | op:Tag{negate}, operand:Index | 同次元 |
| ExprNode power | op:Tag{integer_power}, operand:Index, exponent:signed-native-integer | exponentはint32範囲。指数0は無次元1、負冪は零を含む定義域を拒否 |
| ExprNode analytic | op:Tag{sin,cos,exp,log,sqrt}, operand:Index | sin/cosはangle（SI radへ変換）、exp/logは無次元、sqrtは全次元指数偶数 |

各domainは、必要な全評価点と時刻を含む契約である。範囲外での評価、零を含む除算、logの非正、sqrtの負を拒否する。必要二階微分のあるsqrtは引数の正の下界を要求し、微分不能な零を有限微分としない。式の定義域は、nodeの区間包絡を使用してprepareで検査する。保守包絡で成立が示せなければ区間を明示的に再分割するか未解決を返し、clamp/微小値補完で続行しない。有限入力からoverflowした評価/微分は精度未解決、入力非有限はNonfiniteを区別する。

schema段階は全nodeの型・次元、acyclic、index、出力数、domainと必要変数の適合を検査する。SI段階は十進値/角度/次元を変換し、prepare段階は値と必要な一階/二階微分、定義域、連続性、誤差を検査する。次数のある微分は同じ式から生成する。利用者が別の速度・加速度式を渡して位置と不整合にしてはならない。

DAGのnode数/edge数/最大深さ/指数絶対値/評価点数は明示上限で確保前に検査する。資源計画はnode値、区間値、同時生存する全一階/二階微分、型情報と出力のbytesを含む。中間値の同一性はnodeの参照によるものであり、式文字列の任意eval、ファイル読込、再帰、ループ、条件分岐、外部symbol、乱数、任意関数、callbackを認めない。イベントによる区分はdomain表で明示する。
