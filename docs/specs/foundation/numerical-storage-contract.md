# IM00-4の数値保持型と境界の内部契約

N4-02、Issue #8を具体化する。対象は保持・配列配置・layoutの所有と検査であり、物理kernel、拘束式、解法は追加しない。SI単位の変換は入口の責務を継承し、この型で黙って変換しない。

## 共通の所有と失敗

製品型は `uvlas::numerics` に置く。すべての所有配列はResourcePlanを通す固定長のPlannedBufferで、各型はコピー禁止・移動可能とする。移動元は空になる。所有者の移動・代入・破棄で借用viewの保証が終わる。rvalueからの配列借用はコンパイル時に禁止する。IDとメタデータの個別取得は値で返し、範囲外をInvalidInputとする。

factory名のcopyは明示コピーを意味する。生入力と親状態を書き換えず、返値は入力配列の寿命に依存しない。形状・範囲・配置の不一致はInvalidInput、非有限値と集約中の非有限化はNonfinite、表現不能な寸法と資源不足はResourceExceededとする。各診断に呼出側のOperationを保持する。

複数配列のfactoryが途中で失敗した場合、既に成功した個々の配列確保の計上はResourcePlanの累積規約に従って保持する。破棄した準備配列の計上を巻き戻して予算を複製しない。ソートやID検査の作業配列も計上する。

## 密配列

DenseRealはdouble、DenseComplexはcomplex<double>。shapeはrows、columns、leading_dimensionをIndexで保持する。列優先の成分は `column * leading_dimension + row` に位置する。leading_dimensionはrows以上とし、列間のpaddingも所有する。

非空の物理配列長はleading_dimension×columnsで、積を検査してから入力長と照合する。0×N、N×0は0要素として扱い、途中の大きな積を作らない。paddingを含む入力全体の実部・虚部を検査する。所有後の書込みも可能だが、利用側は外部へ渡す前にvalidateで非有限値を検出する。内部Eigen境界もこの検査を行う。

## CSC疎配列

SparseRealは値double、row indexとcolumn offsetはint64。canonicalなCSCのコピーでは、offset列がcolumns+1要素、先頭0、末尾nnz、単調非減少であることを確認する。各列のrow indexは範囲内かつ厳密昇順とし、重複を拒否する。columns+1自体のoverflowも検査する。0列では静的な0をoffsetの唯一の要素として返し、payloadを確保しない。

tripletからの構築は、列・行・元の入力位置をキーに並べ替える。等しい列・行の加算順は元の入力順に固定する。重複はdoubleのNeumaier補償加算で集約し、途中または結果が非有限なら拒否する。集約後にrow／value／offsetを圧縮する。合計0になった要素も構造的な位置として残す。これは精度昇格の代替ではなく、TwoWord保持とその消費はIM00-5以後の責務である。

## 点と自由度layout

PointSoAはx／y／zと物理PointIdを入力順に所有する。4列の長さを一致させ、全座標の有限性とIDの一意性を検査する。IDの0とuint64最大値は合法。点IDを添字や他の種類のIDへ暗黙変換しない。

DofLayoutは全座標の順序、DofId、ComponentId、座標種別、SIの7基本次元指数、正で有限な成分尺度を所有する。座標種別は位置、ANCF勾配、梁の回転chart、その他のscalarを区別する。ANCF勾配を節点回転と同一視しない。patch／memberとの対応は複数許し、各全座標に少なくとも一つの対応を要する。同一の対応は重複とし拒否する。対応列は全座標index、所有対象の種類、物理IDの順に整列する。

拘束前後のlayoutは、固定pivotの独立座標の選択と順序を所有する。independent_to_fullは独立座標順の全座標indexで、重複・負値・範囲外を拒否する。full_to_independentはその逆対応をoptionalで返す。選ばれなかった全座標は従属であり、値が0という意味ではない。

物理値の非線形写像 $$ q=\psi(r,t) $$ とJ作用は拘束モデルが評価する。このlayoutの添字対応を固定係数の疎行列や物理値の変換関数へ読み替えない。独立座標の列が空でも、全座標がすべて従属のlayoutとして保持する。

## Eigenとの受け渡し

Eigenを使うadapterは `src/numerics/eigen_views.hpp` の非公開ヘッダに置く。公開ヘッダへ依存ライブラリの型を出さない。製品の所有Denseは、寸法・leading dimension・有限性を検査してからOuterStride付きのread-only Mapで借用する。製品のSparseはint64 CSCをそのままread-only Mapで借用する。いずれもコピーしない。

Eigenからの借用は元のMatrix／Mapからscalar、row/column major、inner/outer stride、要求寸法、有限性を検査する。Refへの変換時に一時配列を作ってしまう経路を使わない。外部配列のpaddingの割当終端を検証できないため、元のEigen配列から生spanを返す入口はpacked列優先に限定する。所有Denseのpaddingを持つMapの提供は別に対応する。非packed・row-major・float型・要求shape違いを黙ったコピーや型変換で受理しない。temporaryからの借用は削除したoverloadで拒否する。

## 独立した確認

| 契約 | 確認 |
|---|---|
| 密の列優先・padding・complex | 既知の成分位置、入力変更からの独立、虚部、空形状、paddingの非有限 |
| CSCの圧縮と重複 | 手計算のoffset/index、厳密整数の相殺和、明示0、異常pattern、集約overflow |
| 点の対応 | 入力寿命、全軸の非有限、長さ違い、IDの重複と両端 |
| 自由度layout | 非単位尺度、次元指数、座標種別、複数patch/member、独立順序と逆対応、従属、異常入力 |
| Eigen境界 | pointer一致、既知成分、元のrow-major／float／stride／shapeの拒否、非有限、temporary借用不能 |

これらの確認は保持型の実装適合に限る。物理ソルバーの精度、数値認定、実務規模の性能、実験妥当性、全用途の受入を意味しない。
