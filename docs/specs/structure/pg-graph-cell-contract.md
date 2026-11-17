# PGグラフ・セルの所有と資源

Issue #94の実装補足。[PG係数契約](pg-coefficient-contract.md)2.2/2.3節の式・領域・次数と、[exact代数の資源契約](pg-exact-algebra-resource-contract.md)を変更しない。対象は無次元有理グラフ、全half-cell、明示片側、有理評価と矩形活動map。全辞書、多角形、Q/o/L_ref/厚さ、rank、条件付け、求積、物理参照資格は親 #17 の後続である。

## 所有と片側

`pg_graph.hpp`のPgGraphRequestは同期呼出し中だけx0/w/a/fを借用する。PgGraphはそれらを明示cloneし、n、版、全cellを所有する。出力point、3成分cell多項式、3成分値、矩形mapも所有型で、rvalueのviewは禁止する。viewはownerの移動・代入・破棄まで。移動元の演算/viewはInvalidInput。呼出しはconst入力を変更せず、Planは共有して累積計上する。可変cache/共有Workspace/user callbackを持たないため、別workerは別Planで同じconstグラフを使える。

`PgSide::Interior`は各cellの厳密内部だけに使う。境界j/(2n)のLowerはj−1番cell、Upperはj番cell。0はUpper、1はLowerだけが内向きである。両軸を独立に検査する。from_localもu/vが−1ならUpper、+1ならLowerを要求し、指定cellと一致することを確認する。未知tag、範囲外、矛盾するsideはInvalidInput。locatorは選んだcellとexactなs/t/u/vを返し、境界係数を平均しない。

cell_polynomialsは局所3成分を所有して返す。evaluateは同じsourceのpointだけを使い、chart微分は元ChartST多項式を直接微分・評価する。CellUV微分はchart値へ既定の逆尺度を掛ける。この独立した2経路をcell係数からの直接評価とも照合する。高階orderは各軸の非負・既知tag・入力の検査を先に行う。全成分が次数超過で零なら尺度の冪を生成せずexact零を返し、巨大な正orderを冪のoverflowで誤拒否しない。非零の場合の合計orderは高々4なので、尺度の分母もnative整数の検査付き計算で足りる。

anchor_atはs=0、s側Upper、cs=0からt側を明示して評価し、tのchart微分だけを返す。他cellの左端との相違を検査する。物理座標へ写す入口ではない。

## 内容identityの正規形

GraphContentIdentityはfactoryが所有する実内容からSHA-256で生成する。callerからhash/成功フラグを受ける入口は設けない。型名やobjectの生byte、host limb幅、十進の元字句をhashしない。

正規streamはdomain文字列`uvlas.pg.graph.v1`、版、x0/w/a、f、n、実cell列、anchor定義、UnitSquare領域の順とする。整数長・指数・件数はbig-endian uint64。各有理数は分子sign byte、絶対分子のbyte数と最小big-endian magnitude、正分母のbyte数とmagnitudeを持つ。零はsign0/分子長0/分母1。fはcanonical項数と(a,b,係数)の標準順。実cell列は(cs,ct,linear-index)と両軸の境界numerator/共通denominator2n、尺度4nを標準順に含める。anchorはcs0/u−1と既定版、領域はUnitSquareの版と境界0/1を含める。定数文字列もbyte長を前置する。sourceのfixed limbをbyte単位に読むstreamingであり、別の十進文字列/vectorを作らない。既存Boost::hash2を使い新依存を足さない。

pointとmapはsourceグラフの独立snapshotも所有する。evaluate/validate_forはdigestだけでなく、x0/w/a/f、n、版、実cell列の内容をexact比較する。同じshapeの別graphはStaleGeneration。mapの生成は全cellと共通exact面積から実被覆を確認するprivate factoryだけが行い、任意cell一覧・boolmapを受けない。

このidentityはgraph内容だけである。後段の完全参照identityは、この所有内容にQ/o/L_ref/厚さ/対象領域等の実入力を合成して生成する。本型を完全参照資格へ暗黙変換しない。

## 一つのcallを計画する

Limitsは既存ExactLimits、maximum_cells、stop_tokenを明示する。nから得る全cell数と入力bits/次数/版/sideを最初に検査する。外部受入は全callで1回行い、内部のexact呼出しごとに同じcapを再利用して保持結果を予算外へ出さない。RequireProvenBoundは未証明のためUnsupported、推定/Plan不足とchecked overflowはResourceExceeded、取消はCancelled。内部Failureの詳細を保持する。

計画passは数値を確保せず、既存exact実装の有理演算bit式、polyのmerge_bound、各raw項数を伝播する。Vは返却係数のbit上限、Bはactive exact呼出しの代数作業上限として区別する。quarticのcell代入は $$ U=\sum_{(a,b)\in f}2^{a+b}\le129 $$ 、各raw係数は $$ V_f+20 $$ 以下、既存実装の作業上限は $$ U(V_f+21) $$ 以下。線形成分は既存4則の式を順に適用する。点の変換は入力bit上限Pから高々P+7であり、native16の整数pair生成に必要なB=9も最大へ含める。内部native定数の字数も生成前に検査する。評価は次数超過を先に除き、微分の係数増加を高々12bit、点の冪を高々4Pとして同じ既存式へ渡す。矩形areaのコピーだけは全map snapshotのH_inとcopy専用2整数/resizeを数え、存在しないGCDを加えない。graph全体を複製するcallは、そのclone計画も含める。これらの最大をBとする。

旧契約のW/C/G/D/K/Rを使い、H_inは全借用入力の実capacityを合計する。Lは下表の保持有理数slotの上限、Uはactive exact呼出しのraw項数上限である。

$$
E_{call}=H_{in}+h[(2L+16+6U)C+\max(G,D,K)+R].
$$

active exact出力をL側にも数える保守的な重複を許した推定であり、実heap/RSSや新たなproven boundではない。予測をFailure.resourcesの実測値へ入れない。内部exact呼出しのB/Uがこの計画を超えないことを実装と照合する。zero短絡でactive exact処理が無い場合は、そのscratchを0とする。

| call | Lの有限な保持集合 |
|---|---|
| prepare/clone | x0/w/aとf高々15項の所有出力：18有理数 |
| locate/from_local | source snapshot18＋座標4＋変換中scalar5：27 |
| cell_polynomials | 出力3成分高々20項＋affine置換4項＋scalar11：35 |
| evaluate | 微分2段の高々30項＋出力3＋scalar8＋尺度1：42 |
| anchor_at | 保持するpointのsnapshot/座標22＋evaluate42：64 |
| 矩形map | snapshot18＋面積/総面積2＋native件数/単位1の2：22 |

Planの事前比較は、実装が呼ぶ全allocateの累積上限を $$ \sum_j n_j\,sizeof(T_j) $$ で求め、remaining_bytes以下かを確保前に確認する。予約用のdummy配列や課金を足さず、実allocateだけをPlanへ計上する。privateのRawTerm実sizeofを得るpureなlayout accessorを既存代数のprivate header/sourceへ追加する（算術と呼出し順は変更しない）。S_RをRationalStorage、S_PをPolynomialStorage、S_TをExactTerm、S_rawをRawTermのsizeofとする。

poly clone/make/微分の上限Q(m)は $$ m(S_{raw}+S_T+S_R)+S_P $$ 。cell代入はglobal raw Uと各外項の1+2+…+2^(a+b)を含め、raw累積は $$ 3U-|f| $$ 。それに返却高々15項の $$ 15(S_T+S_R)+S_P $$ を加える。locateはgraph clone・point storageと高々9回のscalar storage、cell抽出は高々11回scalar、affine/linearのQ(2)を3回・Q(3)を1回と上記代入を数える。evaluateは値storage・高々8回scalar・微分Q(|f|)を高々2回。anchorはlocateとevaluateの累積、mapはgraph clone・map storage・高々4回scalar。各call自身のdata storageと全cell列も含める。零係数等で実量が小さくなることはあり、計画上限と実allocated差分は別に返す。

## 固定literalと開始審査

全n/全cellの標準順、境界/中点の正逆、両側とcorner、非単位w=4/5・a=3/20を検査する。flatと総次数4の手書き多項式について位置・両軸mixed/非零4階・次数超過の巨大orderを独立有理literalへ許容差0で比較する。localとchartの尺度差、cs>0の左端とanchorの差を含める。内部sideの根拠はcell indexと±1のlocal座標で確認し、まだ作っていない区分PG辞書の不連続微分を検証したとは言わない。

正規streamの小literalは独立canonical bytes/SHAへ照合し、同値有理表現/移動/cloneの同一性、同shapeのx0/w/a/f/n/版変更の拒否を確認する。矩形は4/16/64cell、各面積1/(4n²)、総面積1とfactory-only型を検査する。入力/旧結果の不変、取消、Plan/外部予算の同値上限と1不足、bad_allocも含める。producer byte同一の隔離control後に片側/尺度/hash入力欠落/活動漏れ等の意味ある変異を検出する。

本書/APIの所有、片側、canonical stream、有限slotと累積payloadをrootが審査してから数値bodyへ進む。通常具体化の範囲であり、新たな物理方式/所有者判断はない。source固定後に対象Debug/Releaseを各1回へまとめ、既存代数は算術を変更しない限り重複実行しない。新しい性能27標本の採取はこの閉包の完了条件に足さない。
