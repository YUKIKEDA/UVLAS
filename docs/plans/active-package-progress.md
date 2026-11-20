# IM-01・IM-03・IM-06の進捗と残作業

2026-11-20時点で、実装中の3パッケージは**全43段階のうち14段階が完了（約33%）**。残り29段階のうち26段階が機能の実装・検証、3段階が費用検証などを含む各パッケージの終了審査。対象のmainは [`5340e6c`](https://github.com/YUKIKEDA/UVLAS/commit/5340e6c)。

| パッケージ | 完了／全体 | 残り | 現在の到達点 |
|---|---:|---:|---|
| [IM-01 構造](#im-01-構造) | 4／12（33%） | 8 | 原辞書と積分分割まで実装済み。n4で全候補の厳密再構成を確認し、四次の2領域で行尺度の共有条件を検証 |
| [IM-03 空力](#im-03-空力) | 5／13（38%） | 8 | 局所P1核、大域幾何、後縁の放出所有まで実装済み。既知流れではPCHIPの係数、区間評価、節点検索、時刻の処理、境界の連続性検査、片側微分の選択を統合 |
| [IM-06 入出力](#im-06-入出力) | 5／18（28%） | 13 | 型付きschemaとSI変換まで実装済み。File／Memoryの取得元との結合と、参照宣言の事前検査を統合 |

割合は[既存の作業分解](next-package-work-breakdown.md)の段階を各1件として数えたもの。段階の成果がmainに入り、その段階の範囲を満たした時点で完了とする。完了欄の分子に、試作・部分実装・レビュー中の変更は加えない。段階ごとに規模が異なり、全43段階には現在の3パッケージの作業だけが含まれる。工数と初回実務適用の進捗を出すには、後続パッケージも含む個別の見積もりが必要になる。統合で到達点が変わったときは、該当行と冒頭の集計をこのページで更新する。

## IM-01 構造

[#17](https://github.com/YUKIKEDA/UVLAS/issues/17) の終了審査まで8段階。局所rankの資源条件と全列再構成が、条件付け以降の先行条件。参照幾何は、グラフと領域分割が揃った範囲から進められる。

| 段階 | 状態 | 到達点と残作業 |
|---|---|---|
| 厳密代数 | 完了 | 有理数・疎多項式の準備。 [#91](https://github.com/YUKIKEDA/UVLAS/pull/91) |
| グラフとセル | 完了 | 有理グラフ、half-cell、矩形の活動map。 [#103](https://github.com/YUKIKEDA/UVLAS/pull/103) |
| 全原辞書と full V | 完了 | アンカー補完、必須場と並進の厳密な写像。 [#113](https://github.com/YUKIKEDA/UVLAS/pull/113) |
| 多角形と全積分分割 | 完了 | 被覆を検査した分割と資源計画。 [分割](https://github.com/YUKIKEDA/UVLAS/pull/128)・[資源計画](https://github.com/YUKIKEDA/UVLAS/pull/149) |
| 局所独立化 | 検証中 | 固定領域の証明器と難条件2列の復元は統合済み。全候補の有限体探索の固定領域実証も [#186](https://github.com/YUKIKEDA/UVLAS/pull/186) で統合済み。n1四次・全背景の全383列の有理復元と全1252行の厳密証明は [#196](https://github.com/YUKIKEDA/UVLAS/pull/196) で統合済みの限定実証。n2四次・全背景でも [#204](https://github.com/YUKIKEDA/UVLAS/pull/204) で全730列の復元と全5008活動行の厳密証明を確認した。[#210](https://github.com/YUKIKEDA/UVLAS/pull/210) では、n4四次U形切断の全1,496候補を復元し、元式8,428,464行と独立証明書の745,008座標を厳密に照合した。[#217](https://github.com/YUKIKEDA/UVLAS/pull/217) でn4平坦1,237候補と四次1,496候補の全背景を検証し、元式の全活動行と独立証明書の全座標を厳密に照合した。[#227](https://github.com/YUKIKEDA/UVLAS/pull/227) では、四次の全背景とU形切断で、共通行尺度が補正桁と残差更新を保つ条件を検証した。共有尺度を使った固定バッチへの接続、残る活動条件の検証、辞書・分割を含む全準備の資源審査を経て、製品実装の開始を判定する。 [数学契約 #133](https://github.com/YUKIKEDA/UVLAS/issues/133)・[資源方式 #140](https://github.com/YUKIKEDA/UVLAS/issues/140) |
| 参照幾何と微分 | 残 | グラフと領域分割を使い、参照Jacobian、可変法線、実法線厚さ、片側を含む高階微分を実装する |
| 保証付き条件付け | 残 | 局所独立化と参照幾何を結合し、幾何内積と基底変換を誤差上界付きで検証する |
| PreparedPatchと参照作用 | 残 | 条件付け済み基底を、層・面・片側の参照作用と反復用Workspaceへつなぐ |
| 全位置traceと拘束 | 残 | 辞書・分割・基底を使い、非整合界面の全厚位置、元拘束、双対量を保持する |
| 固定pivot疎chart | 残 | 元拘束と疎解法から写像と微分・転置作用を実装し、拘束と仕事を照合する |
| 早期R-02 | 残 | 実拘束グラフができた時点で、疎分解のfillと準備・再利用の容量を測る |
| 終了審査 | 残 | V-PGで曲面・切断・非整合を通し、IM-02へ引き渡す |

## IM-03 空力

[#19](https://github.com/YUKIKEDA/UVLAS/issues/19) の終了審査まで8段階。局所核を共有辺・四辺形へ組み立てる作業と既知流れの実装は並行できる。循環とその微分が、全荷重・後流へ進む共通の先行条件。

| 段階 | 状態 | 到達点と残作業 |
|---|---|---|
| 厳密な局所述語 | 完了 | 退化、内外、距離とfar判定。 [#90](https://github.com/YUKIKEDA/UVLAS/pull/90) |
| 所有幾何・ID・query | 完了 | 幾何snapshot、整数接続、厳密な物質位置式。 [#102](https://github.com/YUKIKEDA/UVLAS/pull/102) |
| 包絡と局所packet | 完了 | 位置・幾何包絡と局所適格性。 [#112](https://github.com/YUKIKEDA/UVLAS/pull/112) |
| P1核と全13変数の微分 | 完了 | 3基底の局所誘導核と全Jacobian。 [#129](https://github.com/YUKIKEDA/UVLAS/pull/129) |
| 共有辺・四辺形・全作用 | 残 | 局所核と接続を使い、共有出現の鎖則・相殺、JVPとVJP、双対内積を検証する |
| 大域幾何と後縁所有 | 完了 | 大域接触・不適格の検査、alias運動学、後縁の放出所有。 [#153](https://github.com/YUKIKEDA/UVLAS/pull/153) |
| 既知流れと指定運動 | 実装中 | [PCHIP係数](https://github.com/YUKIKEDA/UVLAS/pull/179)は統合済み。指定区間の値・時間微分も [#187](https://github.com/YUKIKEDA/UVLAS/pull/187) で統合済み。閉定義域内の節点検索と左右区間候補の取得は [#197](https://github.com/YUKIKEDA/UVLAS/pull/197) で統合済み。[#205](https://github.com/YUKIKEDA/UVLAS/pull/205) で、範囲外のRejectとHold、およびHold端点の内外候補の取得を統合した。[#211](https://github.com/YUKIKEDA/UVLAS/pull/211) で、宣言した周期と節点幅の一致を検査し、負時刻を含む周期時刻の厳密な正規化を統合した。[#218](https://github.com/YUKIKEDA/UVLAS/pull/218) で、要求区間、Hold端点のゼロ傾き、周期端の値と時間微分の一致を厳密に検査する部品を統合した。[#228](https://github.com/YUKIKEDA/UVLAS/pull/228) で、微分を評価する側の選択と、Automaticで必要になる導関数の厳密な一致判定を統合した。準備入力の所有と世代照合、評価と再開への接続、既知流れ・運動への組込みが残る |
| 循環と全微分 | 残 | 核、大域幾何、流れと共通求解器を結合し、非貫通の原残差、全束縛・全後流作用、循環微分・随伴を検証する |
| 全ベクトル荷重 | 残 | 循環・共有核・幾何から全荷重と微分を実装する。構造への双対写像にはIM-01／02の成果が必要 |
| 後流値・移流・試行 | 残 | 同じ循環・荷重状態から自由／指定輸送、全履歴、試行の親不変性と全微分を検証する |
| J1状態への接続 | 残 | IM-04側とstep契約を共有し、BDF2と新生台形、輸送別の初期微分、受理済み循環の時間微分、全履歴、一回限りの状態確定を検証する。強連成Newtonの実装はIM-04が担当 |
| 参照と数値予算 | 一部検証済み | [局所核の独立参照](https://github.com/YUKIKEDA/UVLAS/pull/129)を確認済み。全経路を組み立て、V-AERO・V-MULTIの観測量と数値予算を照合する。参照準備は並行可能 |
| 費用と終了審査 | 残 | 全履歴の相互作用数、資源、並列時の結果一致とR-03を確認し、パッケージを審査する |

## IM-06 入出力

[#22](https://github.com/YUKIKEDA/UVLAS/issues/22) の終了審査まで13段階。安全に取得した同じバイト列を、外部入力の形式ごとの読取りに渡す。NastranとGmsh、CLIとPythonは、それぞれの共通基盤が揃い次第並行できる。

| 段階 | 状態 | 到達点と残作業 |
|---|---|---|
| 厳格なscalar変換 | 完了 | 強いIDとIndexへの検査付き変換。 [#89](https://github.com/YUKIKEDA/UVLAS/pull/89) |
| parser資源計画 | 完了 | 固定JSON依存の例外経路を補修し、見積もりと受入条件を実装。 [依存補修](https://github.com/YUKIKEDA/UVLAS/pull/101)・[資源計画](https://github.com/YUKIKEDA/UVLAS/pull/104) |
| 所有文書 | 完了 | 元バイト列、digest、型、十進字句を保持。 [#104](https://github.com/YUKIKEDA/UVLAS/pull/104) |
| 全schema | 完了 | 全caseの型、入れ子、参照と相互条件の検証。 [#130](https://github.com/YUKIKEDA/UVLAS/pull/130) |
| SI・精度変換 | 完了 | 元の厳密値、SI倍率、成分別の精度管理。 [#154](https://github.com/YUKIKEDA/UVLAS/pull/154) |
| 参照snapshot | 実装中 | 宣言収集、パスの正規化、[単一ハンドルの照合](https://github.com/YUKIKEDA/UVLAS/pull/177)は統合済み。Windowsの根相対adapterは [#184](https://github.com/YUKIKEDA/UVLAS/pull/184)、Linuxの根相対adapterは [#195](https://github.com/YUKIKEDA/UVLAS/pull/195) で統合済み。[#202](https://github.com/YUKIKEDA/UVLAS/pull/202) で、マニフェストの元バイト列とdigestを単一のファイルハンドルから取得する部品を統合した。[#219](https://github.com/YUKIKEDA/UVLAS/pull/219) で、Fileは元マニフェストのハンドルと実体ID、Memoryは呼出し側が開いた根を検査済みケースとともに保持するようにした。[#229](https://github.com/YUKIKEDA/UVLAS/pull/229) で、参照ファイルを開く前に、宣言順と位置を保持して参照先を正規化・索引化する部品を統合した。実体とバイト列の照合、依存抽出・循環検査、所有参照表への結合が残る |
| 励振と既知入力 | 残 | schema・SI・参照snapshotと流れ／J1の契約を、時系列・イベント・必要微分へつなぐ |
| parametric入力とprepare | 残 | 型付きcaseを共通のRawModelへ渡す。物理準備には対応するIM-01／02／03の成果が必要 |
| Nastran | 残 | schema・SI・snapshotを入力に、採用カード、座標系、sidecarを解析する |
| Gmsh | 残 | 同じ入力基盤を使い、採用したMSH 4.1 ASCII、物理名、sidecarを解析する |
| 共通参照への変換 | 残 | 外部面を共通の曲面・厚さ・接続へ変換し、質量・仕事などの誤差を評価する。対応するIM-01／02の成果が必要 |
| 結果と診断 | 残 | 物理回復の出力をstatus、来歴、CSV／VTU／PVDと対応付ける |
| checkpoint容器 | 残 | snapshotと状態基盤を使い、形式、整合検査、確定checkpointを守る書込みを実装する |
| 実状態とresume | 残 | 全BDF2ステンシル、受理済みの微分・荷重・後流履歴、探索状態などを保存する。IM-04／05の実状態を使い、同版・同環境での再開結果のビット一致と非互換の拒否を確認する |
| CLI | 残 | 共通処理にvalidate／prepare／run／resume／inspectと終了状態を接続する |
| Python | 残 | 同じC++処理にPythonの辞書・ファイル・実行・再開と所有権を接続する |
| ランタイムと新規環境 | 基盤あり | CLIとbindingの基盤は [IM-00](im00-package-acceptance.md) にある。実際に配布するCLI／wheelの依存・ABIと、開発環境に依存しない導入・実行を確認する |
| 終了審査 | 残 | V-USE・V-STATEで3入力経路とCLI／Pythonの物理結果、両OSの導入・中断・再開を照合する |

schemaのLinux sanitizer実行に関する [#142](https://github.com/YUKIKEDA/UVLAS/issues/142) は継続。数値表の段階数には、同じ段階に属する個々の子Issueを重ねて加えていない。

## 初回実務適用までの位置

IM-00の共通基盤は終了審査を通過した。IM-01の成果をIM-02の材料・梁・接続・構造全体へ渡し、IM-02とIM-03の成果をIM-04の強連成へつなぐ。その後にIM-05の安定性・探索、IM-07の実務規模の数値・性能検証、IM-08の実験・実務受入がある。IM-06は並行しながら、実行と再開の段階で揃った物理状態を接続する。 [全体の依存とゲート](implementation-start-review.md#3-後続パッケージと依存)
