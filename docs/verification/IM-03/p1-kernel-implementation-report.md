# P1局所誘導核の実装と検証

採用した局所P1核を、3基底×3成分と全13変数微分まで実装した。対象は元の局所幾何packetを消費する内部部品である。共有四角形の全16変数、JVP/VJP、大域幾何、循環、後流、荷重、連成と実務受入は既存IM03-0の後続に残る。

製品の数値bodyは `995b170c1bf417cf0583f7b6b90b664b79e46d5e`、primitive境界試験の追加は `0f8fd4195774a3878510090acf075e7cea4a753c` に固定した。後者で変えたのは試験1件だけで、製品bodyと既存試験の内容は保持している。以下はこの2版の実行範囲を分けた記録であり、QA上の数値認定や製品規模の性能を主張しない。

## 元の位置場と全微分を所有する

`P1Basis` は元の幾何owner、query、facetと散布descriptor、contextを保持し、9基底成分と117微分を所有する。成分の単位はそれぞれ $$ \mathrm{m}^{-1} $$ と $$ \mathrm{m}^{-2} $$、微分順序はquery位置3、頂点位置9、半径1である。元の両語と物質位置式をq格子で再構成してから評価し、query位置を二語へ投影しない。

近接枝は境界渦線と一定面内渦層を併用し、理想核の同一平面外点では法線微分を保持する。半径零の微分は採用した右微分である。遠方枝はDuffy変換した正方形上のGauss12²を返し、Gauss8²との差と両則の丸め包絡を全基底・全列で検査する。求積差の診断と丸め誤差包絡は異なる型にした。

固定256bitの外向き区間とchecked1024bit整数を用い、sqrtの257候補、log/atanの128項と剰余を具体化した。入力の既存checked32768bit整数については、MSVCのBoost内部Karatsubaが固定backendでもheap scratchを確保することをnative stackで確認した。32bit limbの学校式積に置き換え、3200/6400/19000bitの密な入力を、別の可変長整数へ厳密対照した。積・既存桁・carryの和が $$ 2^{64}-1 $$ 以下になる根拠を実装に残した。

## 採用契約7節の全行に試験を対応させた

| 採用系列 | 製品の試験と参照 |
|---|---|
| 全基底・近接 | `P1 all 126 components agree...`。独立面積分＋ambient complex-stepの19例で全126成分。`P1 refined near and self...` の2例は全成分を再対照し、精密化参照中心の包絡内包含も観測 |
| 理想外側極限 | 同19例のcoplanar exteriorとline extension。法線方向を含む全13列と半径零の右微分を照合 |
| 自己・incidence | `P1 orientation permutations...` の全6置換×7支持で全基底・全列を評価。非平面quadの両子は19例の独立参照へ対照。LG-F01〜14の理想面上・非incident接触等の入口拒否と両子descriptorは、変更していない既存local geometry packet試験を再利用 |
| 尺度・低語 | `P1 huge common origins...` で巨大原点とlow、尺度 $$ 2^{\pm900} $$、通常doubleで半径二乗がunderflowする例。回転と尺度1e−3/1e3は `P1 geometric covariance...` |
| 分岐・保存 | far等号と両隣接値、移動・回転・尺度、全6置換の全126成分、一定強度の独立Simpson線積分による渦輪復元。全7物質支持×9位置方向×9基底成分の有限差分に、 $$ u_{\xi_i}+\beta_i u_p $$ の全出現鎖則を照合 |
| 演算包絡 | 高精度primitive対照に加え、`P1 primitive dyadic boundaries...` でsqrt零・平方数、指数−9000〜9000、logの微小入力と2冪境界隣接、atanの象限境界。厳密有理interval、subgrid、全幅等号と直下、全基底・全列のfar予算境界も別に検査 |
| 所有・資源 | owner寿命・移動元・入力不変、quota等号/1不足、取消、全実確保故障点のdetails保持、allocator再入、独立workspace並列の全126成分bit一致。準備済み評価は結果1確保のみ。round-to-nearest、FTZ/DAZ、確保callbackによる環境変更も拒否 |

19例の連続面積参照は固定の相対1e−8または絶対1e−10で照合した。代数試験は1e−10、幾何共変と鎖則の比較は1e−8を用いた。参照の細分化・桁数・刻み変更と独立性は [独立参照報告](p1-independent-reference-report.md) にある。

精密化2例で参照中心が丸め包絡に入ったことは観測結果であり、残る連続積分誤差の厳密上界ではない。遠方の丸め包絡は同じ固定Gauss12²を100桁の別実装とcomplex-stepで評価し、全126成分の包含を確かめた。Gaussの20根・重みと全40モーメントは、Fraction/Sturm参照へ製品の実テーブルを照合した。この検査と連続面積分の離散誤差検査を混同しない。

鎖則試験の初稿ではEdge10の期待重みを入力facetの向きで書いたため失敗した。採用済みの昇順PointId規約に従い、期待重みを訂正して全支持を再実行した。製品側の位置式は変更していない。

## clean版で対象D/Rと変異を実行した

環境はWindows AMD64、MSVC 19.51.36260.0、SSE2、`/fp:strict`、Boost 1.92.0。元依存prefixを読み取り、専用Ninja Multi-Config buildを用いた。両版・両構成のbuild-infoは `uncommitted_changes=false` で保存した。

| 固定版と対象 | Debug | Release |
|---|---|---|
| 995b170のP1通常/resourceの2試験 | 2/2成功、28.23秒＋4.08秒 | 2/2成功、3.15秒＋0.51秒 |
| 0f8fd41の追加primitive境界1件 | 66 assertions成功 | 66 assertions成功 |

変更した製品/試験4 TUのclang-tidy 21.1.1と、最後に編集した各TUの追検査で対象診断は0件。第三者headerからの抑制診断総数を、製品の警告として数えていない。変更C++のclang-format 21.1.2 dry-runも成功した。LinuxとsanitizerはPRの対象CIへ委ね、この記録では成功と扱わない。無関係なsuiteや全体CIの待機は行っていない。

995b170の固定複製で無変異controlを成功させた後、低語除去、面内項欠落、半径seed欠落、外側法線微分欠落、境界向き反転、far差無視、半径/全幅混同、incident自己全零、絶対値微分混入、尺度seed誤りの10変異を実行した。全てcompile/link成功後にCatch2が数値試験の失敗を返した。最初の尺度変異は未使用引数でcompile拒否されたため検出件数に含めず、同じ欠陥を引数を参照する形で再実行した。全幅とfarの変異は、製品が実際に使うprivate qualifierへ入れ、そのliteral境界試験で検出した。追加のprimitive境界試験は変異対象bodyを変えていないため、10変異を重複実行していない。

元の生ログはoutに保持し、公開artifactではworkspaceと依存prefixのパスを置換した。試験出力と終了値は保持している。再実行用 [隔離手順](p1-kernel-implementation-v1/run-isolated-mutations.ps1) はcompile databaseから依存prefixを読み取り、新しい出力先でのみ動かす。

## 計画の再利用と資源の境界

式DAGは非確保のcounting passで数え、checkedな全call workとNode/Jet payloadの計画を予約した後に構築する。最大のimmutable shapeをworkspace準備時に保持するため、各評価で729パターン×2の計画を繰り返さない。実payloadはworkspace 3,048,032 bytes、結果6,960 bytesである。外部backendの28,672 bytesは別のforecastであり、allocator overheadやRSSの証明ではない。証明済み上限を要求された経路はUnsupportedにする。

変更前後の固定したDebug計測コピーで、同じ入力を交互3巡、各巡3回、許可された最小番号の論理CPUへ固定して実行した。他担当の重い処理は計測中に停止した。near評価の中央値は226.731→83.9404 ms、farは2025.8929→1907.7752 msだった。nearの計画部分の中央値143.745 msを、準備側へ移した効果を観測した。全18組で126成分の中心両語・半径・指数・stepはbit一致した。元の計測snapshot、各9標本の最小/中央値/最大とphase別内訳はartifactに保存した。

規約に従い、同じ固定コピーを `/O2 /MD /fp:strict /utf-8` でも別々にビルドし、一コア固定・他担当の重い処理停止の下で交互3巡を追加した。Releaseのnearは最良19.8127→11.5759 ms、中央値21.4694→12.3865 ms、観測範囲19.8127〜22.346→11.5759〜13.8172 msだった。farは最良259.6783→253.7131 ms、中央値266.7216→260.0498 ms、範囲259.6783〜273.9994→253.7131〜268.9146 msで重なり、総時間の短縮を認定しない。Releaseでも全18組の126成分の両語・半径・指数・stepはbit一致した。

計測コピーの当時の宣言に合わせ、beforeはstop_tokenの値渡し、afterは参照渡しのheaderを保存した。いずれも数値式は同じで、最終clean版そのものを再計測したものではない。製品全体の性能や実規模のR-03受入を表さない。farの多数回区間Jet評価の費用は残るため、共有/大域作用を実装する後続gateで、要求規模に対する時間・メモリ予算と併せて審査する。

独立レビューでは、所有・元位置場・13seed、外側極限、区間primitiveと剰余、全call計画、details伝播、環境検査、独立参照の分離を確認した。指摘されたdetails欠落・再入code・繰返し計画は本実装で訂正した。採用式からの逸脱はない。
