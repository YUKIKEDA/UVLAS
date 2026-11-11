# 実務ケース・受入計画

2026-11-06。作業5の統合計画。製品の受入結果ではない。[製品要求](product-requirements.md)、[確定精度](acceptance-criteria.md)、採用設計1〜4を継承する。

## A5-01 判定単位

一つのケース版は、入力を生成するソースと引数、参照データ、比較する物理量・座標・時刻、許容差、数値誤差群、資源枠を含む。入力と参照を結果を見る前にhash固定する。外部モデルの校正に使ったデータを独立validationへ重複計上しない。

状態は`Specified / SourceAvailable / ReferenceQualified / ImplementationVerified / NumericallyQualified / ExperimentQualified / Accepted`を別のフラグで持つ。順番の違う種類の証拠を一つの合格へ潰さない。sourceが入手可能でもmeasurement uncertaintyが未評価ならExperimentQualifiedではない。部品の移植が通っても、実務リリースは全必須群の受入まで行わない。

比較の基準と数値予算は下表のまま維持する。実験比較では同じ単位の比較差と根拠付き不確かさを加え、許容値以下を要求する。測定不確かさがない場合はゼロとせず判定保留。異なる誤差を相殺したり、独立性不明の項を二乗和にしない。

| 量 | 製品目標 | 全数値予算 |
|---|---:|---:|
| 合力・moment・静変位 | 5% | 1% |
| ねじれ | 0.25° | 0.05° |
| 対象固有振動数 | 2% | 0.4% |
| フラッタ速度／周波数 | 5% | 1% |
| 一次調和振幅／位相 | 5%／3° | 1%／0.6° |
| 同じ連続体モデルの層別応力／ひずみ | 5% | 1% |

以下の追加指標はこの目標を緩和しない。参照がゼロに近い変位は代表長L、力はF0、momentはF0L、応力は材料弾性尺度E0で絶対量を定める。数値移植試験のゼロ量はそれぞれ1e−10倍を初期上限とし、既存fixtureのより厳しい検査を維持する。実験でノイズ床以下の量は事前に識別し、相対誤差を恣意的な分母で小さくしない。

## A5-02 部品・部分系の必須入力と参照

入力の正本は以下に列挙する**凍結した試験ソース内のモデル構成関数・定数・引数**と、対応resultのsources／artifactsである。同じ数値を別文書へ手写しして二つの正本を作らない。製品移植時にfixture adapterでこれらをcase.jsonへ展開し、原ソースhashと全パラメータを記録する。旧resultの上書きや、元の試験をそのまま製品の合格結果へ転用しない。

| ケース | 要求・入力の正本／参照 | 実装後の必須比較 |
|---|---|---|
| V-PG | GEO-01、STR-01、MAT-01。[R4](verification/paired-reference-r4-results.md)、[後退テーパR5](verification/paired-shell-r5-domain-v2-result.json)の曲面・厚さ・歪み・支持・独立体積／全面力 | 同じ全3Dエネルギー・質量・トレース、格子／求積別の変位・仕事・stress・mode。製造右辺だけで独立荷重を代替しない |
| V-MAT | MAT-01/02、OUT-02。[層kernel](verification/material-layer-kernel-result.json)、[shell統合](verification/material-shell-integration-v2-result.json)、[層粘性](verification/structural-layer-damping-result.json) | 同じ積層・材料軸・初期ひずみ、独立材料点／厚さ積分、内力・接線・散逸・層両側。材料定数未指定を推測しない |
| V-BEAM | STR-03/04。[梁物理v4](verification/beam-se3-physical-v4-result.json)、[断面](verification/beam-section-interface-result.json)。独立Cosserat ODE・Timoshenko伝達ODE・断面Poisson | 純曲げ、大回転先端力、3厚さ、対象mode、偏心質量・慣性・原点変更。polar momentでSaint-Venant torsionを代用しない |
| V-JOINT | STR-02/03/04、BC-01。[接続](verification/connection-adoption-results.md)、[実梁shell](verification/connection-beam-shell-result.json)、[全荷重](verification/structural-combined-result.json) | 実KKT／独立消去、全厚さトレース、交差・段差・剛体／梁frame、反力、wrench、全接線と同じM/K/C |
| V-LOAD | BC-01、LOAD-01。[既知入力](verification/uvlm-known-input-result.json)、[構造時系列／支持](verification/structural-history-state-result.json) | 時刻／座標／補間、指定運動、gust、支持の仕事と散逸、範囲外／重複の拒否 |
| V-RECOVER | OUT-01/02。[全回復](verification/structural-recovery-result.json) | 応力測度・材料軸・層片側、切断面wrenchと基準点変更、各観測の1%数値予算。未提供の梁回復を創作しない |
| V-AERO | PHY-01、GEO-01。[UVLM終了](verification/uvlm-adoption-completion.md)、[最終12観測](verification/high-frequency-final-result.json)、[抗力](verification/drag-common-limit-result.json) | 同じ有限格子／コア解、理想核の独立参照とその不確かさ、20帯・合力／moment・局所荷重。旧未達も回帰する |
| V-MULTI | PHY-02。[面間調和](verification/multiple-surface-harmonic-refined-result.json)、[大域写像](verification/uvlm-global-binding-result.json) | 独立2翼コード、同時相互作用、共有／TE所有、仕事、別幾何／別核の混同拒否 |
| V-J1 | SOL-01/02/03。[実PG全自由後流](verification/pg-uvlm-free-step-result.json)、[BDF2応答](verification/coupled-bdf2-prescribed-result.json)、[J1終了](verification/j1-adoption-completion.md) | 同じ離散平衡、全未知量JVP、全履歴P、非線形微小応答、二次収束、初期境界微分。2方向の照合を全翼精度と呼ばない |
| V-SEARCH | SOL-03、RUN-02。[枝](verification/static-branch-continuation-result.json)、[探索](verification/flutter-search-control-result.json) | 折返し、mode順序交差、cluster、成長率の曖昧さ、2臨界区間、queue完全再開。実翼のmode被覆は別に認定 |
| V-STATE | RUN-01/02、QUAL-03/04。[実連成再開](verification/coupled-checkpoint-result.json)、[IO4](io-execution-adoption-specification.md) | 同一環境bitwise、異OSの許容差、全現在微分と荷重・全履歴、破損・欠落・世代・二重確定・中断保存 |
| V-USE | IO-01/02/03、API-01/02、ENV-01/02、DIST-01/02。[N4／IO4](numerical-foundation-adoption-specification.md) | 手定義の同じ物理モデル→Nastran／Gmsh→変換量、CLI／Python、初期化／再開、fresh環境のsource／CLI／wheel |

表の全行を初回実務適用の必須ゲートにする。部品で式を確かめ、以下の実務統合ケースで組合せと長時間運用を確認する。対応する製品コード・結果が存在するまではImplementationVerifiedにしない。

## A5-03 実験候補S10の具体化

正本は [Zenodo v1](https://zenodo.org/records/19114860)。取得済み：LE／TEのNastran、解析済みmode／figure、readme、Run_Numbers.xlsx、raw ZIPの全126目録、代表4原データ。[数値・単位・所在の診断](verification/s10-reference-metadata-result.json)。

| ケース | 事前に特定した入力・参照 | 判定条件／現在の不足 |
|---|---|---|
| X-S10-GVT-LE | S10_LE.nas、GVT4004のver1〜7。解析済み公表周波数4.31／28.9／37.5Hz、32markerの座標と3mode | 同じfixture・付加物・材料の版、mode対応を照合。周波数2%、数値0.4%。3公表modeだけで全安定性mode被覆としない。系統不確かさとモデル校正履歴は未確定 |
| X-S10-GVT-TE | S10_TE.nas、GVT4003のver1〜6。公表4.27／27.4／38.7Hz | LEと同様。片方を材料／支持校正に使った場合は校正対象と独立検証対象を明記 |
| X-S10-STATIC | TE、Run4012、root AoA3°、名目U20/30/40/50/60/70m/s。最初の20m/s原データを取得 | 名目値でなく記録した動圧・角度・密度の整合した入力。32marker変位・断面ねじれ。支持座標・marker番号・未変形baseline・校正・付着流範囲が認定されるまで参照採用保留 |
| X-S10-FLUTTER | LE、Run4030、root AoA0°、FPS_LE_II、速度掃引。原データ97,063点の1000Hz系と23,296点のcamera系を取得 | 切離し前、微小振幅、定常性・掃引率・成長率の不確かさ、同じ変形平衡・modeを認定。生データの振動開始時刻をそのまま線形臨界としない。速度／周波数5%、数値1% |
| X-S10-DECAY | Run_NumbersのWire_pull対応表を正本に、原ファイル・構成を選定して固定 | GVTや静的ファイルを解放応答へ読み替えない。初期変位・解放・支持・減衰、入力履歴・同期と不確かさの取得が必要。現時点は未認定 |

動圧はreadmeに従いraw値×98,100でPaへ変換する。取得した静的例では11列目とDynP_1000hzが一致し、速度との内部整合から得る密度は約1.225kg/m³。しかし、これは実測密度や圧力較正の不確かさの独立証拠ではない。

readmeの系重量11g、ラン表の10g表記、deckの分散質量は対象物の範囲が同じと確認できていない。差を誤植と断定せずfixture／fuse／payloadを対応付ける。strainのx位置は取得した静的MATに存在せず、ZIP目録にもsensor_loc.mはない。yからxや材料軸を推測して局所ひずみの実験合格にしない。

**S10はSourceAvailableであり、ExperimentQualifiedではない。** この判定を明示したまま採用方式の製品移植を進められるが、実務受入／リリースの代替にはならない。必要資料の公開取得を継続し、公開元への照会が必要になれば外部送信の前に利用者の判断を求める。照会していない現状を、取得済み・回答見込み確定と表現しない。

別版Pazy、S20、論文中の図を無条件にS10の不足値へ移植しない。代替実験を採用する場合も、同じ入力再構成・不確かさ・流れの適用・校正／検証分離を審査する。該当領域の実験参照が認定できなければ、その実務用途は未受入のまま残す。

## A5-04 実務統合入力P-CURVED-BOX

S10だけで積層・初期曲面・複数面・交差接続を済ませない。以下のSI入力を、解析結果を見る前の統合モデル正本とする。これは実験供試体ではなく、同じ3D材料・梁・境界条件の独立FEMと比較する数値検証問題である。

主翼は片持ち半翼、y∈[0,1.2]m、c(y)=0.4−y/6 m、x_LE=y tan20°、a=x−x_LE、0≤a≤c。基準外板は

$$
z_0(x,y)=0.04y+0.20\,a(c(y)-a).
$$

この係数0.20は1/mであり、zはm。参照グラフは二次多項式で、曲率を持つ後退・テーパ領域を厳密に指定する。根元の全厚さ位置トレースを固定する。

外板は8層各0.125mm、[0,+45,−45,90,90,−45,+45,0]°、下から上の順。各層E1=135GPa、E2=E3=10GPa、G12=G13=5GPa、G23=3.8GPa、ν12=ν13=.30、ν23=.35、ρ=1600kg/m³。相反Poisson比はコンプライアンス対称性から作り、全3D正定値性を検査する。材料軸1は局所span方向への参照接平面内投影で定義する。

前桁はa=.25cの外板からz方向へ12.5mm下の参照線を持つ梁、矩形断面4×25mm、E=70GPa、ν=.33、ρ=2700kg/m³。ねじり定数は矩形Saint-Venant解の収束した級数／独立断面Poissonから求める。後桁はa=.65c、深さ25mm、厚さ0.8mmのshell web。ribはy=.3/.6/.9m、a∈[.1c,.9c]、外板から下へ25mm、厚さ0.8mm。同じ等方材を使用し、全厚さの位置接続を採用接続仕様で定義する。梁接合の偏心・frameを保持する。

尾翼は別の片持ち揚力面、y∈[0,.5]m、c_t=.16−.096y m、x_LE,t=1.2+y tan5°、z_t=.15+.10 a_t(c_t−a_t)、同じ積層順・各層0.075mm。根元は同じ剛な試験基準へ固定。厚い胴体の空力を追加せず、両揚力面と両後流を相互に評価する。

構造検証荷重は重力(0,0,−9.80665)m/s²、外板上の基準面積当たり(0,0,10)Paのdead面力と、別ケースの現在面法線圧力10Pa。空力密度1.225kg/m³、root迎角0/3/5°、速度5〜70m/sを5m/s刻みとする。変形平衡が成立しない点も結果として記録し、静的不安定をflutterに数えない。

小規模／中規模の独立3D Hex27比較では同じ層、梁断面、物理接続・面力・支持を使う。格子と求積を3段階以上に変え、外板変位、3断面のwrench、層両側stress、対象modeを照合する。接続特異点から厚さの2倍以内は局所点stressの比較対象から事前除外し、その領域のstress精度を主張しない。合力と断面力の検査からは除外しない。

## A5-05 時間・周期・安定性の評価手順

過渡は同じ確定平衡から初期変位／速度、指定運動、PCHIP既知gustを個別・組合せで実行する。最大変形、ピーク時刻、固定窓の尺度付きL2時系列、各断面wrench、層stressを比較する。過渡の追加目標は振幅・L2誤差5%、ピーク時刻は対象第一周期の1%。対応する全数値予算は1%／周期の0.2%。event前後・初期Taylorと実履歴を区別する。

周期入力は独立参照のk=.05/.2/.8の上下・pitch・gustを用い、共通のexp(+iωt)基準へ変換する。4周期窓の複素係数が2窓連続でその数値予算の10%以内に安定してから8周期を評価する。振幅が記録ノイズ床の10倍を超えない観測の位相は判定対象としない。ノイズ床・対象除外は計算結果を見る前に固定する。

flutterは同じ離散平衡と全履歴Pについて、成長率不確かさ込みで安定／不安定端点を作り、MACと部分空間でbranchを追う。全数値1%の初期配分は構造空間.20%、空力空間.20%、時間.15%、後流保持.15%、コア.10%、求解／固有抽出.05%、速度区間／mode被覆.15%。これらは最大配分で、実測した相関・交差項を含め合計認定する。個別各1%を許す配分ではない。

hを半減するときは保持時間を保ってRを倍増し、同じ物理問題を再構築する。後流保持は1/1.5/2倍、空力・構造は独立3段階、コアは独立感度で評価する。成長率が小さい領域の時間照合は少なくとも10周期・振幅線形性2水準を確認し、成長率の推定不確かさで窓を延長する。資源上限で認定できなければ未解決。周波数順の最初の数modeだけで終了しない。

## A5-06 運用・性能

PERF-FL-01の幾何・材料・支持はP-CURVED-BOXとする。空力は主翼40×100、尾翼20×50の計5,000束縛panel。構造はn=2有界patchを基準に、span／chord分割を固定列(16,4),(16,6),(16,8),(24,4),(24,6),(24,8),(32,4),(32,6),(32,8),(48,8),(64,8)の順でprepareし、全部材の拘束後自由度が95,000〜105,000となる最初の分割を採る。尾翼とweb／ribの分割は主翼の物理セル幅以下に切上げる。解を見て自由度を調整しない。該当分割がなければモデル生成失敗として列を改版し、物理計算前に再固定する。

代表性能プロファイルは**64GB搭載・プロセス48GiB上限**とする。要求された32〜64GB範囲内の主測定機を具体化した技術判断であり、現48GiB機を64GB機と呼ばない。32GB／24GiBプロファイルは別途容量成立範囲を報告する。CPU、物理／logicalコア、clock／電源、OS、storage、thread上限、依存・最適化・全解決設定を測定前に固定する。

速度は5〜70m/sを5m/s刻みの14点を全て評価し、不確かさ込みの符号変化区間を全て探索する。初期走査14点＋追加最大34点、合計48点を資源上限とする。上限到達を合格とせず未解決にする。初期走査の不安定島を飛ばしていないか中点／mode被覆で検査する。無交差の場合は「区間内の臨界なし」とし、臨界速度の性能受入ケースとしては不適格。受入にはこのモデルか、同じ構造・複数面条件の改版で少なくとも一つの認定した臨界区間が必要で、改版は性能結果の計測前に独立参照で固定する。

h、後流時間、求積、固有抽出範囲はA5-05で精度認定したresolved-caseを先に固定してから時間計測する。保持時間やmode数の上限を性能都合だけで短縮しない。checkpointは各速度完了時と探索queue更新時、過渡照合は10周期ごと。最後の状態と選定平衡の可視化・層stress、全診断のI/Oを含める。

総計測は入力からprepare、全速度平衡・線形化・探索、必要な微小過渡照合、再試行、出力完了まで12時間以内。事前準備cacheは同じ実行で作った費用を含める。ビルド・外部参照全V&Vのみを除外する。現時点でこの性能ケースはSpecifiedであり、メッシュ生成・臨界区間・誤差予算・製品性能の認定は未実施。

最初の性能gateは実装初期の実PG結合graph、5,000密空力、自由後流とP作用の個別計測。演算／メモリ内訳から12時間／48GiBが成立しないなら高速化・前処理・近似の誤差認定を設計し直す。未検証FMMを既に採用したことにせず、直接基準経路を残す。品質・ANCF・実務範囲を下げて性能合格にしない。

## A5-07 終了条件と未受入の扱い

作業5で閉じるのは、実装後に実行するケース入力／参照の構成方法、評価、閾値、取得できた資料、不足時の停止規約である。全製品の実験合格や性能合格を実装前に要求しない。表の凍結fixtureには既に独立参照があり、最初の実装パッケージの試験入力を供給できる。

S10の実験参照認定、実PG全規模・全自由後流の性能、fresh環境の配布、全必須用途のAcceptedは未完了。これらをリスク台帳と実装チケットへ引き継ぐ。判定保留・未達を公開仕様の品質達成へ昇格させず、実務適用の初回リリースは全必須ゲートが成立するまで行わない。
