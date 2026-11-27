# 局所幾何packetの包絡・資格契約

Issue #109、IM03-3。[局所幾何契約](local-geometry-contract.md)5〜10節と[所有geometry/query契約](owned-geometry-query-contract.md)のA2を具体化する。UI-D01/03/04/05/06、UW-D07、N4-02/03/06/08/09を維持する。本書は実装開始用の型・演算・資源境界であり、実行結果は検証記録へ分ける。

## 1. 入力・所有結果・資格

内部入口は `prepare_local_packet(const OwnedGeometry&, const ResolvedQuery&, FacetId, LocalAccuracy, const GeometryLimits&, LocalPacketWorkspace&, ResourcePlan&, stop_token) -> AeroResult<LocalGeometryPacket>` とする。現在配置を評価し、参照配置の非退化はA1の所有snapshotに結び付ける。外部公開数学APIとcore共有失敗型は増やさない。

`LocalAccuracy` は有限非負の幅要求 `tau` と明示した `Binary64 / TwoWord` の返却表現を持ち、無意味なdefaultを設けない。出力表現の選択は元位置式を投影する操作ではない。両語・物質式を最後まで厳密に消費し、選択した中心表現への丸めを包絡へ含める。要求幅を満たせない場合はPrecisionUnresolvedとする。

`ScaledEnclosure` は有限中心（表現variant）、有限非負の絶対誤差半径、符号付き2冪指数を持ち、真値を `2^exponent * (center ± radius)` へ包含する。TwoWord中心は両語の厳密和である。零は指数0へ正規化する。factory以外に有効な包絡を構築させず、指数付き表現を普通のdoubleへ暗黙変換しない。

packetはsnapshotとqueryの不変ownerを保持し、Contextと内容identity、対象Facet/PhysicalFace/Band、元QueryId/定義/正規支持、半径と精度要求を所有する。導出量は位置3、法線3、面積ベクトル3、長さ1、距離1、半径1の12包絡である。far関係、最近実体、物理資格、元位置式の微分散布descriptor、work/資源報告を固定要素として持つ。借用getterはlvalue限定とし、移動元は無効として拒否する。

物理資格は `variant<VerifiedIncidence, Detached>` とする。VerifiedIncidenceの頂点/辺/面は正規支持と対象facetの接続からのみ作る。同じ座標であることをincidence証明にしない。支持が対象facetに属する場合は厳密な位置関係とも照合し、不一致はGeometryInvalidとする。支持が属さず、点が有限閉三角形上なら非incident接触としてGeometryInvalid。面外または同一平面の有限面外はDetachedである。VerifiedIncidenceかつ半径の厳密和が零なら理想自己としてGeometryInvalidとする。

最近実体はA0の全最小候補maskから、頂点→辺→面の次元順と正規ID順で選ぶ。A0の閉集合maskは保持し、物理支持と混同しない。packetは常に `GlobalGeometryNotChecked` を所有する。13/16変数のdescriptorはA1の元query定義・物質座標とfacet scatterを保持し、第三重みを二語へ先丸めしない。核のJacobian値はまだ生成しない。

## 2. 評価式と平方根の包含

A1と同じ $$ q=2^{-2148} $$ の厳密整数点を使う。差をe,f,r、整数外積をc、整数二乗面積をJ2とする。位置成分、法線成分、面積成分、長さ、距離、半径は、符号付き $$ s\sqrt{N/D}\,2^E $$ （N非負、D正）に揃えて評価できる。

| 量 | s、N、D、E |
|---|---|
| 位置成分I | sign(I), I², 1, −2148 |
| 法線成分c_i | sign(c_i), c_i², J2, 0 |
| 面積成分c_i/2 | sign(c_i), c_i², 1, −4297 |
| 長さ | +, longest_edge_squared, 1, −2148 |
| 距離 | +, A0 distance numerator, denominator, −2148 |
| 半径R（2^-1074単位の厳密和） | +, R², 1, −1074 |

位置・面積・半径は元の符号付きdyadic整数から直接中心と半径を作る。上位53bitを最初の語へ、その厳密な残差の上位53bitを低語へ移し、残った整数を半径へ上向き変換する。二語間の指数差を106bit以内へ制限しない。巨大原点と小さい局所差の両語を、256fraction-bitの連続仮数へ先丸めしない。Binary64要求では低語相当も半径へ含める。共通指数はhigh/lowをともに有限に保持できる範囲から選び、指数0で両語を保持できる場合は不要に正規化しない。低語または半径の残差が選択指数のsubnormal格子より小さい場合も、厳密残差を最小正数以上の上向き半径へ含める。零へ丸めて正確としない。幅を満たせなければPrecisionUnresolvedとする。

N=0は厳密零を返す。平方根を必要とする非零量では $$ k=\lfloor\log_2(N/D)\rfloor $$ をbit長と整数比較から求め、 $$ t=\lfloor k/2\rfloor $$ とする。負の奇数は負方向へ丸める。正規化した $$ z=(N/D)2^{-2t} $$ は[1,4)に入る。

$$ A/B=(N/D)2^{512-2t} $$ となるよう、シフトの正負に応じてNまたはDだけを左シフトする。A/Bに残る共通の2冪因子は、両整数のtrailing zeroの小さい方だけを左右から除く。正の同じ因子を消すので、以下の全比較と最終等号は不変である。257bit整数sを上位bitから逐次構成し、各候補に対して $$ s^2B\le A $$ を厳密比較する。257回後のsは $$ s/2^{256}\le\sqrt z<(s+1)/2^{256} $$ を与える。等号 $$ s^2B=A $$ なら上下端を同じにする。平方根の符号を比較する許容差はない。各候補の判定は、先に整数商 $q=\lfloor A/B\rfloor$ を求めて $s^2\le q$ と比較しても同じになる。元の等号は $s^2=q$ かつ剰余が零の場合に成立する。

この有理区間から指定表現の有限中心h（またはh+l）を作る。中心との差と区間半幅を厳密な整数差で評価し、上向きに丸めた有限double半径へ含める。指数E+tを別保持し、中間scalbnや平方根doubleのunderflow/overflowで真値を失わない。符号が負なら区間端を反転する。256fraction-bitは内部の固定精度であり、利用者のtauを緩める根拠ではない。包含できても幅超過なら失敗する。tau=0では幅が厳密零の場合だけ成功する。

## 3. 幅認定と整数容量

幅は半径の2倍である。位置・距離・長さの幅をell、面積をell²、法線を1で正規化してtau以下を要求する。半径も長さと同じ尺度で包絡する。比較は元の正の整数ell²と各出力半径・tauの厳密dyadic値を使い、必要なら両辺を二乗する。近似したellの下端を零へ丸めて比を捏造しない。

A0の表より、点成分3173bit未満、c成分6354bit以内、J2は12708bit以内、距離Nは最大19062bit、Dは12708bit以内である。平方根正規化のシフトを含むA/Bは、一般に `max(bit(N),bit(D))+514` 以内、候補s²との積はさらに514bitを加えても20090bit以内である。保守容量21000bitは既存checked32768bitに収まる。反復中も元N/Dを無制限に積み重ねない。

幅認定では各dyadic量をsignificandと指数に分け、正の整数同士のbit長を先に比較する。同じbit長へ整列した場合だけ交差積を実体化する。出力半径とtauは各53bit significandのdoubleだが、実装では2^-1074単位の整数へ復号するため、各整数は最大2098bitである。幅認定の積はこの復号整数を使い、共通の十進/浮動小数点近似へ置き換えない。二語中心は語ごとに53bitであり、両語の間の指数差を連続106bitへ制限しない。元dyadic中心の厳密な残差は元入力の3173bit域と出力語の範囲で比較する。半径をr、共通指数をEとすると、length比較は $$ 4r^2 2^{2E}\le\tau^2\ell^2 $$ 、area比較は $$ 4r^2 2^{2E}\le\tau^2\ell^4 $$ 、法線比較は $$ 2r2^E\le\tau $$ である。コードも幅2rを認定し、半径rだけを認定しない。area比較右辺はtau復号整数の平方とell²の平方を掛け、2×2098+2×6354=16904bit以内である。左辺は半径復号整数の平方と係数4で4198bit以内となる。長さ比較右辺は10550bit以内、法線比較は2098bitの整数同士である。比較の指数整列も大きい方のbit長以内なので、保守容量21000bitを超えない。指数差をそのまま巨大左シフトへ変えず、整数容量を超えるshiftを要求しない。checked演算のoverflowはPrecisionUnresolvedとして返し、有限包絡を捏造しない。

## 4. payload・作業・外部確保

LocalPacketWorkspaceは既存GeometryWorkspaceと同じlease境界を使い、追加の固定 `RootScratch` をResourcePlanで準備する。RootScratchは23個のExactInteger（N,D,A,B,root,candidate,square,product,value,consumed,remainder,radius,temporary、e/f/cの各3成分、J2）を持つ。配列数は固定で、byteは実sizeofを使う。packet結果は別の `PlannedBuffer<LocalPacketStorage>(1)` を持つ。A0の所有分類結果も同じ呼出しのresult_planへ一時的に確保するため、一回の累積payloadは `ExactTriangleResult::payload_bytes()+LocalGeometryPacket::payload_bytes()` である。分類結果を破棄してもPlanの累積計上は戻さない。snapshot/query所有者の保持と固定owner headerは既存A1と同じprocess overheadとして扱う。

既存A1の単独入口はgeometry予算とpredicate予算を別に検査する。A2ではGeometryLimits.max_stepsを一回全体の上限、PredicateLimits.max_stepsを述語subcapとする。A1実課金後の残量とsubcapの小さい方をA0へ渡し、checked合算した実課金を包絡へ引き継ぐ。A2の制御失敗では内側cause/detailsを保持し、未完了の幾何原因を表すtyped problemを付けない。

数値準備の全体stepは、既存A1/A0の実課金に加え、12量それぞれについて平方根候補257回、変換1回、幅認定1回を各1stepとする。厳密零でも同じ上限予算を入口で非課金admissionし、実行した候補/変換/認定だけ実課金する。追加上限は3108stepである。各候補前に取消とchecked累積stepを確認する。payload不足・work不足・取消の既存優先関係を維持し、未訪問の物理不備を推測しない。

32768bit/et_offの積は逐次実行し、既存の固定Boost外部scratch forecast `7L limb = 28672byte` を再利用する。除算・平方根をBoostの別algorithmへ委譲せず、上記bit試行と比較で評価する。新しい同時積・動的精度型・再帰階層を重ねない。allocator overhead/RSSのhard capとは呼ばず、proven boundなしのstrict要求は既存と同じUnsupportedとする。結果確保と外部scratchのbad_allocは既存cause/detailsを保ち、leaseを解放する。

任意の計測要求は数学のaccuracyと分離したmetadataである。disabledを既定とし、時計読みを行わない。enabledでは既存StageTimingで分類、接触資格、包絡、および包絡内のroot/幅認定を記録する。root/幅は包絡時間へ内包されるため合計へ重複加算しない。計測の有無で中心/半径/資格/実stepを変えない。準備段階の観測を物理kernelや実務全体の性能へ一般化しない。

## 5. 検証と後続境界

独立oracleはtest-only cpp_rationalで元語と物質式から幾何量を構成し、製品の平方根routineを共有しない整数平方根囲いを使う。真値包含に加え、固定tau以下の幅、正規ID/符号/最小実体を確認する。LG-F01〜14、全query/支持、両quad子、全6頂点置換、alias同値別分解、2^±900、巨大原点＋低語を含める。通常tau=1e-10、共変1e-8、厳密分類は許容差0である。半径なら予算内だが幅2rでは予算外となる境界、幅と予算が等しい境界、低語/subnormal残差の上向き包含を独立対照にする。

独立workspace並列、同じworkspace再入、context/content変更、旧packet不変、移動元、payload/work上限、取消、結果/Boostの確保失敗を対象にする。隔離byte同一controlと低語・包絡・接触・自己ID・global資格の変異から実検出を保存する。影響するA0/A1と新targetを固定sourceのD/Rでまとめ、旧記録は上書きしない。

P1核の解析式/遠方求積/速度/全微分、大域交差/TE、ANCF写像、循環/荷重/後流の資格は後続である。後続はこのpacketの元位置式・全包絡と同じcore/contextを消費し、包絡を無視して中心だけへ置換しない。
