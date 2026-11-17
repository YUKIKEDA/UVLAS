# exact代数の型・資源・検証の補足契約

2026-11-17。Issue #87、[PG係数契約](pg-coefficient-contract.md)2.1/5/6/7節の内部具体化。[主シェル採用](main-shell-adopted-specification.md)の空間・範囲・精度を変えない。対象は参照準備用の有理数と二変数疎多項式であり、PG辞書、活動域、条件付け、PreparedPatch、物理状態の認定ではない。本書の固定Boost sourceと生存領域の対応は独立レビュー済みである。

## 所有する型と演算を閉じる

`ExactRational`は既約p/q、正のq、零0/1を所有する。整数pairは符号付きASCII十進整数だけを受け、leading zeroとplusを値として正規化する。有限十進入口は符号、少なくとも1数字、小数点とe/E指数を受け、`.5`と`5.`も有効とする。空白、区切り文字、hex、非有限表記を読み替えない。Binary64入口は語を復号し、十進0.1と二進0.1を区別する。整数pair・十進・二進で同じ値ならcanonicalな分子分母が一致する。

結果はmove-onlyで、明示cloneだけが予算付きコピーを行う。移動元はexact zero。`ExactPolynomial2<ChartST/CellUV>`は座標を型で区別し、異なる座標同士の加減を許さない。合成だけは変換先の型を明示して置換する。Boost型・遅延式を公開APIへ出さない。公開関数は利用者向けの英語ドキュメントを持ち、内部仕様IDだけを説明にしない。

多項式はa昇順/b昇順、重複exact集約、零項除去の項列を所有する。空列が零で、零のdegreeは不在。保存するa、b、a+bは全てIndex範囲内と事前に確定する。各軸次数と総次数の上限は入力・生成するraw項に適用し、相殺後だけ小さくなることを理由に中間の上限を無視しない。

有理数4則・比較、多項式の加減積・scalar倍・合成・明示軸の高階微分・点評価・下端付き積分を閉じる。未知axis、負指数／微分階数、分母零はInvalidInput、非有限Binary64はNonfinite。積分は指定下端で零になる定数を含め、0階微分は所有clone、次数を超える微分はexact zero。合成の各有効項について、変換先各軸／総次数を $$ a\deg(P)+b\deg(Q) $$ のchecked算術で確保前に検査する。零置換で消える外項は全入力検査後に短絡し、不要な巨大pow／raw項数を生成しない。

結果・項・診断は自分の領域を所有する。readonly borrowはownerの移動／代入／破棄まで、一時ownerからのborrowは禁止。状態／Modelを入力しないこの部品は、存在しない親世代を0で作らない。既存`core::Result`と`guard_entry`、Operation::PrepareModelを使い、内側Failure詳細とResourcePlanの累積payload規約を保つ。Failure時は入力と旧結果を変更しない。

## 数学bitと確保byteを分ける

Limitsはinput字数、集約後項数、raw項数、各軸／総次数、返却整数bit、作業整数bit、外部メモリpolicyを明示する。0を無制限としない。通常値・作業値・項数・指数・checked byteの不足はResourceExceeded。exact小項を切り捨てない。

有理数加減の中間分子bit上限は $$ \max(n_1+d_2,n_2+d_1)+1 $$ 、分母は $$ d_1+d_2 $$ 。乗算は $$ n_1+n_2,d_1+d_2 $$ 、除算は $$ n_1+d_2,d_1+n_2 $$ 。数学的な式と元operandの最大bitから作業bit Bを定める。multi-termの集約ではraw項数Uと各係数bit上限bから $$ U(b+1) $$ の保守的上限を使う。分子／分母を片方だけ数えない。

以下のbyte forecastはBoost 1.92の可変長`cpp_int`と`rational_adaptor`を、et_off・fresh非alias出力・有限lifetimeで使う経路に限定する。fixed整数や任意の式templateへ一般化しない。係数4・5・6はsourceの領域式から導き、任意倍率を置かない。

$$
w=\mathrm{CHAR\_BIT}\,\mathrm{sizeof(limb)},\quad
h=\mathrm{sizeof(limb)},\quad W=\lceil B/w\rceil+2.
$$

+2はcarry・長除算・Lehmerの余分なlimbの作業枠である。Bをalgebraic_work_bits、word丸めとguardを含むWwをplanned_work_bitsとして別に返し、max_work_bitsには後者を比較する。word幅32/64により受入計画が異なり得るが、exact値のcapや物理精度は変更しない。実測した最大bit数と呼ばない。Iを固定backendのinline limb数とし、C(W)=0、R(W)=0（W≤I）、C(W)=4W、R(W)=W（W>I）とする。Cは新規整数slotのheap capacity、Rは1回のresize中に新領域と同時に残る旧heapである。既存operandはsizeでなく実capacityからH_inを得る。相殺で零になっても過去の大capacityを隠さない。

| nested処理 | 追加limb-payloadのforecast |
|---|---|
| 非alias除算 | D=2C。remainder rとlong-division t |
| GCD | G=6(W+1)+5C。共有poolと、swap後にowningとなり得るU/V/t3個、Lehmer fallbackと除算tの2個 |
| variable Karatsuba | K=5W（W≥60）、他0。再帰ごとにheapを追加せず、最外shared作業域を数える。実operandの短絡で狭める実装も可 |

外側の整数slot数をNとし、一般の有理演算は $$ H_{in}+h[NC+\max(G,D,K)+R] $$ byteを比較する。加減N=7（gcd,t1〜t4と出力2）、乗除N=6（gcd_left/right,t1/t2と出力2）、正規化factoryN=6（parsed pair、出力pair、g/t）である。正規化だけならnested=max(G,D)、比較N=2ならnested=K、cloneN=2ならnested=0を使う。known pathに存在しないGCDを必須量として表示しない。

polyの親処理では同時にglobal raw出力・現在展開・置換先の3配列が残り得るので、係数2整数×3配列で6U slotを数える。固定slot16は累積、第1／第2power、active powerのresult/baseの計10整数と、nested rational乗算の出力2＋helper4の6整数を合算したもの。nested G/D/Kは別の一時最大値として加える。sortは指数比較だけを行い、normalize/公開時にこれより大きいlive集合を持たない。scalar処理ごとにcapをresetして蓄積raw係数を予算外へ出さない。

| 固定source／関数 | 領域の根拠 |
|---|---|
| cpp_int.hpp:205〜224、366〜393 `cpp_int_base::resize` | Iの定義、max(4*old,requested)の増設、新確保→copy→旧解放 |
| cpp_int.hpp:396〜485 | normalizeはsizeだけを縮め、copy/move/swapによりcapacityの残留・移譲がある |
| cpp_int/multiply.hpp:87〜89、225〜266、391〜450 | cutoff40limb、5s、300limb未満のstack／以上のheap、非alias出力条件 |
| cpp_int/misc.hpp:678〜803、882〜1075、1126〜1228 | GCDの6t pool、Lehmer aliases、fallback後のowner移譲とlocal t |
| cpp_int/divide.hpp:17〜145、185〜223、455〜465 | r/t、y.size()+shift+1のresize、alias時の追加copy |
| rational_adaptor.hpp:59〜88、611〜695、759〜819、986〜1008、254〜290 | normalize／加減／乗除／比較の外側live整数 |
| cpp_int/cpp_int_config.hpp:104〜141 | 構成でlimb32/64、decimal block9/18桁が異なる |
| cpp_int.hpp:1752〜1927、cpp_int/bitwise.hpp:330〜375/464〜527 | decimal block演算とshift。直接文字列代入のoctal解釈は使わない |

整数十進は検査済みdigitsのblock累積で構築し、元入力にNULを要求せず、別のtext heapを作らない。複数桁積と加算はnative limb内の値から行う。十進指数の10の冪は反復二乗の有限slotを使い、不要なpower scratchを破棄してから正規化する。zeroでも字句／指数の検査を省かない。Binary64は整数仮数とshiftで復号する。p/q全体のBoost文字列parserや、混在backend演算による非canonical0/kへ依存しない。

## 受入policyと実測を混ぜない

`ExactEstimateAdmission{bytes}`は上記scopeの同時limb-payload forecastをcap以下で明示的に受け入れる。allocator管理領域、CRT例外構築、stack、他の生存Model、プロセス基礎メモリまでのhard capではない。`ExactRequireProvenBound`のproven値は初版では不在で、確保前にUnsupportedを返す。forecast overflow／不足はResourceExceeded、bad_allocも同じcodeで返す。推定byteをFailure.resourcesのobservedへ入れず、不明な実量は省略する。

所有storage／term／raw配列の実sizeof×個数はResourcePlanへ確保前に課金し、外部整数limbと二重保証しない。text scratchが将来必要になれば、それも所有payload側へ計上し、外部式への追加を0とする。実装でcopyやalias／文字列／遅延式を増やす場合は計画を更新し、そのまま現在のslot式へ押し込まない。

空polyは入力・axis・階数・limitsタグの検査後にno-scratchへ短絡し、B=planned work=0とする。公開返値は生成callの外部予測byte、代数的bit上限B、guard付きplanned work bit、raw項数を所有する。これらは観測された最大値でない。timingはEnabled時のPreparation scopeだけを記録し、Disabledは時計読取りを行わず、数学値のcanonical byteを変えない。exact演算は参照準備用の可変長確保を伴い、反復物理kernelのheap0達成を主張しない。

## 独立参照と最初の終了条件

[PG契約](pg-coefficient-contract.md)のliteral系列を、版付き`tests/structure/fixtures/exact-polynomial-v1.json`へ持つ。原式・入力・手計算canonical分子／分母・標準項列が正本で、製品factoryで期待値を生成しない。実際の数値は演算でも使うprivate storage accessorから観測し、手書き文字列へ許容差0で比較する。製品に試験専用分岐を追加しない。

十進／Binary64 0.1、2^53超、最小subnormal／最大finite、非零下端積分、両軸高階微分、重複／相殺、合成次数の事前拒否、cross-frame型拒否を覆う。大整数copyとowned storageの確保失敗は独立test processの標準new置換で注入し、事前の不正／予算拒否はnew試行0と旧値不変を確認する。置換hookは製品へlinkしない。

最終固定版について対象Debug/Release、独立変異controlと検出、対象format/tidy、fixture hash／来歴、計測scopeと未実行OSを記録する。既存の数値suite／凍結R&Dを再実行せず、Linuxの未実行値を成功としない。親Issue #17の残るPG全能力とR-02は閉じない。
