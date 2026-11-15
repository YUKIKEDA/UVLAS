# Python結合：所有権と失敗の契約

IM00-10は、IO4-02のPython境界を、実装済みの密行列保持と状態・資格・確定の製品処理へ接続する。根拠はN4-03、N4-08と[状態契約](state-contract.md)、[保持契約](numerical-storage-contract.md)。解析ソルバー、物理状態の配列形式、ケースのJSON schema、ファイルのrun/resume、wheel配布は定義しない。数値処理と物理表現は各consumerの実装範囲、ケースのschema・run/resume・wheelはIM-06の範囲に残る。

## 入口と返却値

開発用拡張は `_uvlas`。`DenseReal`、`DenseComplex`、`Model`、`TrialState`、`AcceptedState`、`EvaluationReport`と`QualificationPacket`を、対応するC++製品処理へつなぐ。試作の `compute(n)` やJSON readerは移植しない。

製品処理は `Result` を返す。`ok` が真なら `value`、偽なら `failure` に所有する値が入る。失敗の `value` へのアクセスはPythonの `ValueError`。`Failure` は17種類の名前を区別し、operationと、存在するID・資源量・最後の確定状態を保持する。未知の診断は `None` とし、IDゼロで代用しない。

Pythonの呼出し構文・クラス型が違う場合は `TypeError`、Python/NumPyオブジェクトの確保に失敗した場合はPython自身の例外があり得る。これを製品の計算失敗と混同しない。C++の処理入口の例外は `ResourceExceeded` または `InvalidInput` に変換する。ID・世代・byte予算を表すPython整数はdoubleへ変換せず、負数、64bit上限超過、bool、非整数を `InvalidInput` で拒否する。

## 配列をいつコピーするか

| 入口・出力 | 配置・型 | 所有と変更 |
|---|---|---|
| `DenseReal.copy` | native float64、2次元、Fortran連続、必要なalignment | GIL保持中に製品の所有配列へコピー。非有限は拒否 |
| `DenseComplex.copy` | native complex128、同上 | 実部・虚部とも有限。精度・単位・型を暗黙変換しない |
| Model定義・Trial payload入力 | uint8、1次元、C連続 | GIL保持中に製品がコピーする。物理float配列と解釈しない |
| 密行列view、Model定義、Accepted payload | 同じ型と配置 | 読取り専用のzero-copy view。baseが不変な所有者を保持する |
| `TrialState.payload_copy(plan)` | 独立したuint8配列 | 可変Trialへの借用を返さず、planで予算照合した `PlannedBuffer` へ明示コピー |

配列入力への暗黙のdtype変換、C/F配置変換、非連続配列の詰め直しは行わない。必要なら呼出し側がNumPyで明示変換する。空の配列も所有者をbaseに持つ。入力を保持するのは同期のコピー完了までで、その後のPython側の変更は所有済み結果を変えない。コピー中に別のnative writerが同じメモリを書き換えないことは呼出し側の責任である。

viewはWRITEABLEを解除するだけでなく、baseの所有者がwritable bufferを公開しない。`setflags(write=True)`、派生view、memoryviewから書込みを再開できない。Trialは置換・キャンセル・確定できるためlive viewを公開しない。結果やPlanを破棄しても、返却済みのviewと独立snapshotは保持される。

## GIL、同時使用、成功の公開

入力の検査とコピーをGIL保持中に終え、所有済みC++配列の検査、候補hash、qualify、commitの間だけGILを解放する。その区間にはPython callback、Python属性操作、Python参照カウントの変更を入れない。

ResourcePlan、Trial、資格packetには所有者ごとにmutexを持たせ、各所有者のgetterと変更処理で共用する。残量・有効性の取得も含め、全入口をtry-lockで拒否できるようにする。複数所有者は容量3以下を検査した上で `std::less<std::mutex*>` の全順序に挿入ソートし、重複を除いて非blocking取得する。途中で取れなければ、取得済みlockをその場で解放して `InvalidInput` を返す。GIL保持中にworkerの終了をmutexで待たない。単一の通常GIL interpreterだけを対象とし、free-threaded CPythonと複数interpreterは認定しない。

commitはPythonのResult、成功所有者、holderを先に確保する。native commitの成功後には `noexcept` なAcceptedStateのmove/emplaceだけを残す。置換・cancelもPython Resultとその参照を変更前に準備する。包装用の確保失敗で候補だけが消費され、Accepted結果が失われる順序を作らない。

Planは成功配列の累積予算であり、配列の破棄で返金しない。Python包装の失敗で既に作られた配列が破棄された場合も、この意味を維持する。Pythonオブジェクトheader、固定長metadata、アロケータoverheadはpayload予算に含まれない。独立snapshotの可変長native payloadは必ず計上する。

## 評価の意味を変えない

`TrialState.target()` は評価を始める前に取得する。`EvaluationReport.create(model, evaluated_target, bounds)` は、そのtargetと、呼出し側が計算した4種類の保守的上界を保存する。順序は残差・幾何・時間・精度。`None` は未評価。上界の妥当性は物理評価器の責任で、このbindingは支配方程式を評価しない。完了した評価を別候補のtargetへ付け替えない。

qualifyとcommitは実候補のhashと全targetをC++製品処理で照合する。別候補、古いrevision、未評価・未達、二重確定は拒否する。成功・失敗・cancelのいずれでもAccepted親は変更しない。

## 決定と試験

CTestの `python.ownership` は `tests/python/test_ownership.py` を実行する。期待値は固定のNumPy配列と明示したID・世代・状態遷移から作り、主実装で再生成しない。バイト比較に許容差を加えない。

| 決定 | 確認と負対照 |
|---|---|
| 所有・readonly | real/complexと4空形状、元入力の変更、所有者GC、derived view・memoryview、WRITEABLE再有効化の拒否 |
| 厳密な入口 | dtype、非連続、逆stride、rank、endian、alignment、非有限、None、Python型、uint64全範囲 |
| 失敗の保存 | 17カテゴリ・12operationの固定export表、資源診断と未知情報、4カテゴリそれぞれの未評価・未達 |
| 状態 | 別候補・変更後のreport、古いpacket、二重確定、snapshotの独立、親不変、予算失敗とcancel |
| GIL・排他 | 実候補のhash中にPythonを実行し、同じTrialのgetterと変更を非blocking拒否。元入力を変更してもnative候補は不変 |

試作 `rd/foundation/numerics/experiments/cpp-io-binding` からは、ownerをNumPy baseに保持する形とGIL区間の考え方を引き継ぐ。試作の成功ログを製品の検証結果には使わず、製品の処理と所有者を上の試験で確認し直す。
