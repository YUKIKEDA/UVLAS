# 数値基盤・具体型・主要処理の採用仕様

2026-11-06。版1.0。作業4の採用仕様。製品コードは未実装。
構造・空力・連成の数式は [構造](../structure/structural-adopted-specification.md)、[UVLM](../aero/uvlm-integrated-specification.md)、[J1](../coupled/coupled-adopted-specification.md) と初期微分追加規約を継承する。本書で物理式・品質目標を変更しない。

## N4-01 依存と再現可能な構築

| 用途 | 採用経路 | 制限・代替を検討する条件 |
|---|---|---|
| 密実／複素、疎行列、分解 | Eigen 5.0.1、ヘッダのみ、`EIGEN_MPL2_ONLY` | float64基準経路。対称性や正定値性を確認せずLLTへ渡さない |
| 非対称・複素固有値 | Spectra 1.2.0、作用素インターフェース | 部分スペクトルの収束を全モードの安定性証明としない |
| 精密参照、正確な線形拘束、SHA256 | Boost 1.92.0 Multiprecision／Hash2 | 有理数は参照係数・階数等の有限問題に限定。実数の非線形幾何が正確有理数になるとはしない |
| コンパイラ | MSVC 19.44以降の必要C++23機能、GCC 14.4.0 | 現在の構成時検査を継承。`std::expected`の必要機能も検査 |

公式配布のURLとSHA256を `build/numerical-deps/header-sources.lock.json` に固定する。Eigen／SpectraはMPL-2.0、BoostはBSL-1.0。自身のソースのApache-2.0と第三者の通知を別管理する。外部CMakeによる未固定版の自動取得を呼ばず、ロック済みヘッダを指定する。取得物・派生物・バイナリのハッシュを監査する。

基準経路はEigen内部の単一スレッド数値カーネル。実運用は独立な組立タスクを並列化し、可変Workspaceを共有しない。BLAS、疎直接法の別バックエンド、FMM等は未採用。性能計測で必要になれば同じ原作用素と誤差予算を満たすことを独立に検査して追加する。未検証の高速化を性能見積もりに織り込まない。

## N4-02 型・配列・単位

`Index`は符号付き64bit、ファイル上のIDは符号なし64bit。構造自由度、空力パネル、後流頂点、後流面、層、試行、モデル世代、数値世代は別の強い型とする。別型を整数が等しいだけで混用しない。外部32bit APIへ渡す場合は乗算を含む範囲検査後に変換する。

| 値 | 具体的な保持・意味 |
|---|---|
| `DenseReal/Complex` | 列優先、float64／complex<float64>、行列次元と先頭次元を保持 |
| `SparseReal` | CSC、値float64、row index／column offsetともint64、重複集約後に圧縮 |
| `PointSoA` | x／y／zの所有配列、各配列長と物理ID列を一致させる |
| `DofLayout` | 成分順、単位、尺度、patch／部材との対応、拘束前後の写像を所有 |
| `PrecisionBuffer` | `Binary64`または`TwoWord{hi,lo}`。hiとloは同形状、有限、合計が意味する値 |
| `ErrorCertificate` | 観測ID、予算、評価範囲、数値世代、評価表現、達成／未達／未評価 |
| `State` | 物理時刻、segment、確定番号、q/v、循環・全後流・履歴、現在の受理微分と荷重、探索状態 |

SIを内部単位とし、角度はrad。ANCF係数を節点回転と同一視しない。PG参照係数の寸法は採用基底の規約で固定し、梁の回転チャートとは別の成分タグを持つ。CSV・外部格子の単位変換は入口で一度だけ行い、変換と元単位を保存する。

空配列・0×0系は合法な恒等処理として定義する。空の拘束、空の後流、ゼロ右辺をNaNやゼロ除算で処理しない。寸法積・確保byte数はchecked arithmetic、確保前にResourcePlanと照合する。

## N4-03 所有権と確定

`Model`は検査・構築後に不変。`AcceptedState`と`TrialState`は別型。試行は親確定番号とModelFingerprintを保持する。処理は親状態をconstで受け、所有する候補と診断を返す。`Workspace`は呼出側所有・一度に一呼出しのみ使用し、返した状態や結果をその内部の一時配列へ参照させない。

借用`span<const T>`は同期呼出しの間だけ有効。Pythonへ返す結果は所有オブジェクトで、NumPy viewには所有者を保持させread-onlyとする。入力の非連続配列・型不一致は境界で明示コピーまたは拒否し、黙った単位変換をしない。計算中のGIL解放、キャンセルは確定境界で扱い、C++例外をPythonへ越境させない。

キャッシュキーはモデル、幾何、拘束チャート、疎パターン、数値行列、尺度、精度、後流政策、時刻区間、バックエンドの世代から作る。アドレスや配列長だけでは再利用しない。右辺のみ変更の場合も行列世代を照合する。質量の定数性は採用されたシェルブロックに限定し、梁・偏心慣性・チャートによる変化を消さない。

受理は全残差・幾何・時間・精度条件を満たす資格packetを生成してから、親番号と候補hashを再照合し一回だけ行う。失敗・キャンセル時の最後の確定状態は不変。二重確定・別親への確定・資格取得後の変更は拒否する。

## N4-04 大域拘束とPGの疎性

採用シェルの有界patchごとに参照係数、Ritz持上げ、正確線形制約の階数・包含を認定する。製品で全10万自由度の密RREF、密nullspace、密質量平方根を作らない。局所演算の密配列はpatchに限定し、大域接続は疎scatterとする。局所認定情報の再利用キーに参照幾何・材料・境界政策を含める。

J1の独立座標rと写像q=ψ(r,t)を保つ。チャートは固定pivotの独立／従属座標と疎拘束方程式 C(q,t)=0 により陰的に構成する。従属Jacobianの分解を共有し、ψ、J作用、Jᵀ作用、時間微分、二階／必要高階微分を陰関数の線形求解で評価する。Jやψの大域密行列を保存しない。Jᵀは同じ拘束Jacobianの転置求解で得る。

線形拘束の正確な重複は係数とトポロジーから除去し履歴を残す。非線形接続の階数不明・チャート特異・参照曲面の不適格は明示失敗。小特異値を恣意的に落とす擬似逆で正常化しない。チャートを変更する場合はJ1の全履歴変換／segment再始動規約に従い、現在状態だけを変換しない。

全域のKKTを解く経路は原問題との照合・反力復元に使用可能。ただしそれをrのBDF2からqのBDF2へ無断で置換する経路にしない。拘束仕事、反力の基底・符号・単位を同時に検査する。

## N4-05 線形・非線形求解

| 原系 | 基準解法 | 必須の照合 |
|---|---|---|
| UVLM非対称実行列 | pivot付き密LU | 非貫通の原作用素、共有辺・同じ核・精度 |
| 複素後縁更新 | 実LU＋複素小Schur／Woodbury | 先頭列とTの全依存、小Schurの条件、元の複素系 |
| 構造の正定値静的接線 | 対称性・正定値性確認後に疎LLT | 負曲率・剛体modeを隠さない |
| 不定KKT・一般動的接線 | Eigen SparseLU、COLAMD、64bit CSC | pivot・有限値・元残差、拘束・反力 |
| 強連成Newton | 全ブロックJVP、密空力LU／疎構造の固定ブロック前処理＋再始動GMRES | 線形精度を全残差評価から認定。失敗時は明示全ブロック小ケースと照合 |

GMRESは同一線形求解中に前処理を固定する。更新は新しい外側Newton反復で行う。再始動幅64、最大512作用を初期資源設定とし、これは物理受入精度の代用ではない。停止して残差が未達なら未収束。逐次緩和だけで強連成の方式を置き換えない。

Newtonは最大30反復、補正line search最大12回、初期Armijo係数1e−4。幾何不適格候補は親に影響させず短縮する。残差だけでなく尺度付き補正と観測変化も評価する。30回で失敗なら原因と最後の確定状態を返す。刻みを変更するならJ1固定h区間の再構築へ戻り、過去の曲がった後流を無検証に再標本化しない。

残差は原系の成分別物理尺度で評価する。非貫通は1e−11、拘束位置／速度はそれぞれ代表長／代表速度で1e−11、力釣合いは代表力で1e−10を初期の上限とする。ゼロ右辺でもModelから先に定めた正の尺度を使用する。観測の解法予算を超える場合はさらに厳しくする。小さい残差を前進誤差の保証としない。

各観測の全数値予算に対し、反復・求解／固有抽出の合計はその5%を初期割当とする。これは全数値1%に5%を足す意味ではなく、残りの空間・時間・後流・コア等と合算する。既存のより厳しいケース仕様を緩和しない。

## N4-06 精度の保持と消費

[線形精度契約](linear-solve-accuracy-contract.md)を継承。Aの評価誤差、格納Aに対する求解、保持解、返却投影、実際の消費側を別々に認定する。既存の4096格子のfloat64投影未達を成功へ変更しない。

基準LUはfloat64。原残差は誤差なしtwo-productと補償加算で評価し、最大5回の反復改良を行う。通常の相殺は精密参照と照合する。停滞時はtwo-word保持解へ昇格し、低位補正を次の補正で失わない。低位を捨てる`to_double()`は明示投影操作で、投影後の原残差・必要観測・消費者の適合性を再検査する。

TwoWord値を受け取る行列作用、内積、幾何・渦核・荷重・履歴更新は精度政策を引数に取り、hi/loの両方を消費する専用経路を実装対象とする。Boost binary100を独立参照および小規模fallbackに用いる。大規模を無条件に任意精度密分解するfallbackは設けない。未実装または資源枠内で認定できない消費経路は`PrecisionUnsupported/Unresolved`で拒否する。内部だけ精密に解いてfloat64へ丸め、正常な時刻状態へ確定してはならない。

精度の昇格は配列の型とcheckpointに記録し、下流のsupported precisionを検査する。全対応を実装し、既存の投影未達例を受入試験へ移植するまでは、その製品経路を完成扱いにしない。精密値が必要な条件の製品性能は別に測定する。

## N4-07 固有値・安定性

小規模の完全スペクトルはEigenの一般実／複素および対称一般化問題で照合する。大規模はSpectraの一般Arnoldi、同一連成ステップPの作用、必要領域のshift-invertを採用する。全状態の密Pを作らない。

P作用は各方向について、同じ確定平衡・履歴・Jnewの分解を使って implicit step derivative を解く。随伴は同じ行列の転置／共役転置求解から構成し、FDを大規模製品の主作用にしない。

$$
P=C_p-C_xJ_{new}^{-1}R_p,\qquad
P^*y=C_p^*y-R_p^*J_{new}^{-*}C_x^*y.
$$

初期の探索はnev=24、ncv=64とし、元残差・再始動／ncv倍増・shift領域追加・枝／部分空間追跡で対象の被覆を評価する。資源枠が許す範囲で拡大し、被覆が未認定なら安定／臨界値を返さず未解決にする。初期24個だけで全翼の安定を断定しない。左／右固有vector、固有値感度と残差による臨界区間の不確かさも保存する。

元固有残差は1e−10以下を初期上限とし、臨界速度／周波数の割当予算から必要な厳しさを追加する。重根・近接modeは部分空間で扱う。負の構造固有値、静的発散、算法寄生根、wake-only root、フラッタをJ1の規約で分類する。

## N4-08 処理I/Fとエラー

以下は実装用の内部契約であり、公開ABIを固定しない。各返値は所有する結果と診断。`expected<Result, Failure>`を基本とし、Failureにはcode、処理、モデル／親／試行世代、問題成分・ID、必要／実測資源、最後の確定状態の識別を含める。

| 処理 | 入力 → 所有する返値 | 変更・不変条件 |
|---|---|---|
| `prepare_model` | RawModel,Units,Policies → Model,PreparationReport | 参照basis／幾何／材料／拘束rankを認定。生入力を変更しない |
| `evaluate_structure` | Model,TrialKinematics,DerivativeRequest → EnergyForceMass,Recovery,JVP/VJP | 全材料・慣性・予ひずみ・支持を含む。指定状態以外のcacheを混用しない |
| `evaluate_aero` | Model,TrialGeometry,Wake,Flow,TransportDerivative → AeroPacket | 核・非貫通・循環・全荷重・必要微分を同じ世代にする |
| `initialize` | Model,InitialData,TransportPolicy → QualifiedInitialState | 初期付加質量・境界弱流束、物理／算法履歴を区別 |
| `evaluate_step` | Model,Parent,Unknowns,h → ResidualPacket,Linearization | 全r/u/Γ/freeX、同じ荷重微分・誕生を評価 |
| `solve_linear` | Operator,Rhs,Scales,Precision,Limits → Solution,Certificate | 原系／保持／返却／消費者を検査。親状態を変更しない |
| `advance_trial` | Model,Parent,StepPolicy → QualifiedTrial | Newton失敗で部分状態を確定しない |
| `commit` | Parent,QualifiedTrial → AcceptedState | 親番号／hash照合、一回限り |
| `linearize_equilibrium` | Model,QualifiedEquilibrium,Chart,History → StepOperator,AdjointOperator | 同じ離散平衡・物理IDと年齢slotの対応 |
| `search_flutter` | Model,SearchPlan,SearchCheckpoint → SearchResult,Checkpoint | mode／枝／全queue、臨界不確かさ、完了済み条件を保存 |
| `recover` | Model,AcceptedState,ObservationPlan → OwnedResults | 全荷重と受理済みDΓを使用。BDF1へ再計算しない |

Failure.codeは`InvalidInput, Unsupported, GeometryInvalid, RankUnresolved, ChartSingular, Nonfinite, FactorizationFailed, LinearUnconverged, NonlinearUnconverged, EigenCoverageUnresolved, PrecisionUnsupported, PrecisionUnresolved, ResourceExceeded, Cancelled, StaleGeneration, CheckpointCorrupt, CheckpointIncompatible`を区別する。ライブラリ内部例外・bad_allocは入口で変換し、失敗を空の正常結果で代替しない。

## N4-09 容量と性能の責任境界

5,000空力未知量の実密行列1枚は200MB。原AとLUと補償作業を別配列で持つ場合、その枚数を計上する。構造10万自由度は疎値＋64bit index＋分解fill＋拘束チャート＋局所patchを実際のgraphで数える。100kの三重対角試験や2D合成疎行列をPG実構造の性能合格にしない。

後流の状態数はpanel数とは別。各線形化で同じJnew分解を再利用し、Arnoldi basisは全状態次元×ncvで見積もる。自由後流の直接相互作用は二乗コストを隠さず記録する。指定後流の利用はその対象ケースの感度・物理適合性で認定し、性能のためだけに無検証で置き換えない。

32GB搭載・プロセス24GiBと64GB搭載・プロセス48GiBの二つのResourceProfileを用意し、どちらも入出力・checkpoint・固有抽出・fallbackを含め事前確保を検査する。性能の最終受入機と上限は作業5の代表ケースで確定する。現機の約48GiB搭載を64GB機での性能達成と呼ばない。

8〜12時間要求を満たす見積もりは速度点数、非線形反復、線形作用回数、後流状態数、準備・I/Oを積み上げる。数値基盤の容量試験だけから一晩達成を宣言しない。実装前は方式・容量上の成立条件と最初の性能ゲートを固定し、製品性能の実測完了は実装後とする。

## N4-10 実装移植時の必須検査

両OSの密実／複素更新、不定KKT、負の一般化固有値、100k疎、一般複素Arnoldi、実PG／自由後流の全履歴行列、two-product／有理rank／区間平方根／SHAの基準試験を保持する。さらに製品化時に、疎チャートとJ/Jᵀの随伴一致、TwoWord全消費者、世代不一致・予算超過・失敗時不変、並列Workspace、Python寿命、実格子I/O、再開、代表graph性能を追加する。これらの将来試験を今回の通過済み項目に数えない。

公式取得元： [Eigen](https://libeigen.gitlab.io/)、[Spectra](https://spectralib.org/)、[Boost 1.92.0](https://www.boost.org/releases/1.92.0/)。採用版の実ソースと取得ロックを仕様の根拠とし、最新版へ黙って追従しない。
