# R&Dの記録

2026-10-29から2026-11-06までの実装前R&Dの設計メモ、検証記録、実験スクリプト。製品の仕様は [docs/specs/](../docs/specs/) にあり、ここにあるのは、その仕様に至った経緯と根拠、採用しなかった候補である。**ここにある式やコードを、現在の採用仕様と取り違えない。**

各系統のディレクトリは次の形をしている。

| 場所 | 内容 |
|---|---|
| `<family>/*.md` | 定式化、契約、候補の設計メモ |
| `<family>/records/` | 検証記録（実験仕様 `.md`、結果・監査 `.json`） |
| `<family>/experiments/` | 記録を作ったスクリプト |
| `<family>/fixtures/`、`figures/` | 入力deckと図 |

## 系統と結論

結論の根拠は各系統の記録。数値は [実装前の経緯](planning/preimplementation-history.md) の該当行から引いた。右端はメモ／記録／スクリプトの件数。

### 構造：主シェル（調べた順）

| 系統 | 結論 | 件数 |
|---|---|---|
| [ancf-3423](structure/main-shell/ancf-3423/) | 不採用。層別・曲面で反例 | 2／3／0 |
| [q9-plate](structure/main-shell/q9-plate/) | 不採用。3833型Q9は歪んだ薄肉で未達 | 10／48／23 |
| [q9-curved](structure/main-shell/q9-curved/) | 不採用。曲面のひずみ拡張、EAS21 | 10／63／37 |
| [q9-edge-trace](structure/main-shell/q9-edge-trace/) | 不採用。辺モーメント、トレース | 6／39／27 |
| [q9-dynamics](structure/main-shell/q9-dynamics/) | Q9上で行った予荷重・過渡・再開の試作。時間法は製品仕様ではない（製品はJ1） | 6／27／9 |
| [q3](structure/main-shell/q3/) | 不採用。四次場の誤差5.53%／32.75%、薄肉増幅11.46倍 | 5／24／17 |
| [mq3](structure/main-shell/mq3/) | 不採用。歪み柱面で59.63%が残る（M-Q3、M-Q3I1） | 9／71／40 |
| [gsp-q3](structure/main-shell/gsp-q3/) | 不採用。曲面四次誤差33.0%／59.6% | 5／38／17 |
| [mixed-order-150](structure/main-shell/mixed-order-150/) | 不採用。薄肉増幅15.51倍 | 1／17／9 |
| [q4-projection](structure/main-shell/q4-projection/) | 不採用。6構成とも歪みで14.86%が残る | 7／41／33 |
| [mortar](structure/main-shell/mortar/) | 不採用。仕事適合mortarと共同tyingでも14.86% | 2／22／10 |
| [hybrid-gradient](structure/main-shell/hybrid-gradient/) | 不採用。S3+に偽零方向3、S4は増幅21.36倍 | 1／18／11 |
| [polynomial-trace](structure/main-shell/polynomial-trace/) | 不採用。増幅26.3倍／83.6倍 | 1／28／17 |
| [reference-work](structure/main-shell/reference-work/) | 不採用。RW-1〜3、残差補正後も増幅76.04倍 | 6／74／50 |
| [moment-enriched](structure/main-shell/moment-enriched/) | 不採用。ME-1は指数場の誤差が細分化で減らない | 4／35／31 |
| [mixed-mean-stress](structure/main-shell/mixed-mean-stress/) | 不採用。MMS-1は増幅27.00〜38.99倍 | 4／63／44 |
| [geometry-compatible](structure/main-shell/geometry-compatible/) | PGの前段の候補。採否は系統内の記録を読む | 3／26／25 |
| [paired-geometry](structure/main-shell/paired-geometry/) | **採用**。参照係数PG-ANCF（R1〜R5） | 6／78／127 |
| [review](structure/main-shell/review/) | 原因診断、定式化の審査、進め方と見通しの履歴 | 14／2／0 |

### 構造：主シェル以外

| 系統 | 結論 | 件数 |
|---|---|---|
| [material](structure/material/) | **採用**。3D材料、積層、初期ひずみ | 3／11／6 |
| [beam](structure/beam/) | **採用**。SE(3)梁 | 0／14／18 |
| [connection](structure/connection/) | **採用**。接続、支持、偏心慣性 | 0／11／14 |
| [integration](structure/integration/) | **採用**。荷重、減衰、回復、構造統合 | 0／10／11 |
| [conventions](structure/conventions/) | 厚さ座標と参照幾何の初期契約 | 2／0／0 |

### 空力

| 系統 | 結論 | 件数 |
|---|---|---|
| [kernel](aero/kernel/) | **採用**。有限3D誘導と局所微分（2.1） | 3／16／14 |
| [wake-loads](aero/wake-loads/) | **採用**。後流、組立、状態、荷重（2.2） | 8／30／20 |
| [finite-wing](aero/finite-wing/) | 定常の有限翼と抗力の参照（2.3の一部） | 1／16／24 |
| [unsteady-reference](aero/unsteady-reference/) | 非定常の独立参照と誤差予算（2.3） | 0／41／51 |
| [multi-surface](aero/multi-surface/) | 面間相互作用と大域写像 | 0／8／10 |
| [adoption](aero/adoption/) | UVLM統合採用（2.4）の終了判断 | 0／4／2 |
| [avl-reference](aero/avl-reference/) | AVLによる参照診断。処理と監査は成立、**精度受入の参照には不採用** | 5／159／120、deck 97 |

### 連成、基盤、受入

| 系統 | 結論 | 件数 |
|---|---|---|
| [coupled/j1](coupled/j1/) | **採用**。物理速度BDF2、新生点台形、初期整合、全履歴線形化、再開 | 1／24／24 |
| [foundation/build](foundation/build/) | ビルド基盤の検証記録（MSVC 14.44、GCC 14.4） | 0／2／0 |
| [foundation/numerics](foundation/numerics/) | **採用**。数値依存の限定試作（C++） | 0／4／13 |
| [acceptance/plan](acceptance/plan/) | 受入計画と開始審査の監査 | 0／5／3 |
| [acceptance/s10](acceptance/s10/) | 公開実験S10の資料診断。SourceAvailableで、ExperimentQualifiedではない | 0／2／3 |
| [planning](planning/) | 実装前の計画と台帳：当時の現在地、実装前設計計画、設計決定台帳、要求と証拠の対応、初期のバックログと検討資料、知識パッケージの検証記録。**現在の計画ではない**（現在は [ロードマップ](../docs/roadmap.md) と [開始審査](../docs/plans/implementation-start-review.md)） | 12／0／0 |

「採用」は方式と契約の採用であり、製品の実装・数値認定・実務受入を意味しない。

## 整理の方法と確かめ方

R&Dの時点では、全てが三つの平坦なディレクトリにあった。2026-11-11に系統別へ並べ直した（[ADR 0003](../docs/adr/0003-rd-material-sorted-by-family.md)）。

| 基準点 | 意味 |
|---|---|
| タグ `rd-baseline-2026-11-06` | 並べ直す前の状態。元のパスとバイトのまま。旧マニフェストはここで検証できる |
| [tools/rd/relocation_plan.py](../tools/rd/relocation_plan.py) | 各ファイルの行き先を決める規則。手で維持するページ（`MAINTAINED`）と、後から削除したファイル（`REMOVED`）もここに書く |
| [registry.json](registry.json) | ファイルごとの旧パス、新パス、両時点のSHA-256 |

配置を変えるときは、`relocation_plan.py` を直して `reorganize.py reissue` を実行する。手で `git mv` しない。

```text
python tools/rd/reorganize.py verify --rev HEAD                        # ベースライン＋ツール＝現在のツリー、をバイトで再現
python tools/rd/check_tree.py                                          # リンク、コンパイル、ハッシュ対（Python 3.11以上）
python tools/rd/evidence.py check                                      # 取得済みの証拠と索引の照合
```

### 2026-11-11に確かめたこと

- 再配置ツールが書く2,238ファイル（`rd/`、`docs/`、Skillなど）が、ベースラインに再配置ツールを適用した結果とバイト一致。手で維持する4ページ（`README.md`、`AGENTS.md`、`docs/build.md`、`docs/engineering/README.md`）と、動かしていない製品ファイルは比較の対象外。
- Markdownの相対リンク切れ0。実験スクリプト全てがコンパイルできる。
- 記録の「パス：SHA-256」の対のうち、リポジトリ内のファイルを指す6,085組、ファイル名だけで指す611組、`build/` のファイルを指す7,436組（開発機の実ファイルと照合）が一致。
- 既存の証拠整合ツールで、終了判断13件から辿れる1,576ファイル・174マニフェストが整合。
- `audit_*`、`check_*` など477本を、ベースラインと整理後の両方で実行して終了状態を比べた。466本は同じ結果、11本が違った。違いのうち、負荷による時間切れのずれを除くと、原因は兄弟スクリプトとルートの探し方の書き換え漏れ（`Path(__file__).with_name()`、`.parent/`、`source.parents[2]`、同じディレクトリのモジュールという判定）で、再配置ツールを直した。

### 確かめていないこと

- **上の修正の後、スクリプトを再実行していない。** 修正が効いていることは、該当の記法が残っていないことの検索と、`rd_paths` の単体の動作でしか確かめていない。再実行しなかった理由は次節。
- 数値実験の本体（`run_*`、`prepare_*` など365本）は、どちらのツリーでも実行していない。
- ハッシュの値がキーとパスの対になっていない箇所（`source_sha256` など12,739組）は、対応するファイルを機械的に特定できないため照合していない。再配置ツールは、ベースラインのファイルと同じハッシュの値を、書き換え後のハッシュへ置き換えている。

## スクリプトを実行するときの注意

**R&Dのスクリプトを、証拠の実物に対して実行しない。** 多くのスクリプトは `build/analysis/` へ結果を書き、既存のファイルを上書きする。証拠はgitの外にあり、上書きは戻せない。この件の後、証拠の置き場所を定めた（下の「大きい証拠ファイル」）。

2026-11-11の比較実行で、これを実際に起こした。`build/analysis/` の既存18ファイルを上書きし、33ファイルを新しく作った。18件のうち元の内容を失ったのは5件である。

| 対象 | 対処と結果 |
|---|---|
| ハッシュで引用されていた4件のうち2件（`finite-wing-2-3/drag-quadrature-m4-n24.json`、`-m16-n96.json`） | 計測時間の欄が変わった。記録に全文が埋め込まれていたため、そこから書き戻した。記録済みのSHA-256と一致 |
| ハッシュで引用されていた残り2件（`drag-common-limit/far-reference-series.json`、`finite-wing-2-3/drag-elliptic-operator.json`） | 失われていない。出力に時刻を含まず、上書き後の内容が記録済みのSHA-256と一致 |
| 同名の記録がある9件（`q9-*-result.json` 6件、`3833-bending-distortion-result.json` ほか2件） | ベースラインの同名の記録の内容で書き戻した。元の `build/` の内容と同じだったかは確かめられない。再生成された内容との違いは、計測時間・メモリの欄と手順書のパスだけだった |
| ハッシュの引用も同名の記録もない5件（`3833-rank-result.json`、`q9-load-release-progress.json`、`q9-multielement-prestress-progress.json`、`q9-prestress-equilibrium-result.json`、`uvlm-trailing-edge-stopped.json`） | **元の内容は失われた。** いまあるのは2026-11-11に再生成された内容。4件は不採用になった3833型／Q9の試作出力、1件はUVLMの途中停止時の診断出力で、採用仕様の根拠や受入計画のfixtureには入っていない。元と同じ数値だったかは確かめられない |
| 新しく作られた33件 | 削除した |

実行が必要なときは次のようにする。

- 必要なファイルを `out/rd-evidence/build/` から作業ツリーの `build/` へ複製し、それに対して実行する。ジャンクションやシンボリックリンクを張らない。
- `close_*` と `sync_*` は実行しない。README等の現状表示を、当時の文面へ書き戻す。
- `audit_*` は記録を上書きすることがある。実行した後は `git status` で確かめて戻す。
- Python 3.11以上が要る。当時の環境は `build/design-fem-env`（Python 3.12、NumPy 2.3.5、SciPy 1.16.3）。
- スクリプトは記録や他のスクリプトを [rd_paths.py](rd_paths.py) を通して名前で探す。系統のディレクトリを動かしても、スクリプトを直す必要はない。

## 大きい証拠ファイル

記録がハッシュで引用するファイルのうち、大きいもの（大半は `.npz` の数値配列）は、リポジトリに入っていない。リポジトリにあるのは索引 [evidence-index.json](evidence-index.json) で、引用される名前（`build/analysis/…`）ごとに、SHA-256、バイト数、置き場所を持つ。

| 種類 | 件数 | 置き場所 |
|---|---:|---|
| 自前の解析データ、ログ、AVL由来のファイル | 2,355 | GitHubのリリース `rd-evidence-2026-11` のアセット2本（1.73GiBと1.65GiB）。索引の `bundle` |
| 第三者の物（依存の配布アーカイブ、公開実験S10のデータ、AeroSandboxとPteraのインストール物） | 424 | 保存しない。取得元を索引の `source` に書く |
| ビルド生成物（数値試作の実行ファイル） | 2 | 保存しない。ビルドの経路を含む。`rd/foundation/numerics/experiments/cpp-numerics` から作り直す |

```text
python tools/rd/evidence.py summary   # 索引の内訳
python tools/rd/evidence.py fetch     # 束を取得し、全ファイルをSHA-256で照合する（out/rd-evidence/ へ）
python tools/rd/evidence.py check     # 取得済みのファイルを照合する
```

- `out/rd-evidence/` は取得したものの置き場で、消してよい。もう一度取得できる。
- 記録の引用ハッシュと索引の照合は `tools/rd/check_tree.py` が行う。実体を取得していなくても照合できる。
- スクリプトを走らせるときは、必要なファイルを `out/rd-evidence/build/` から作業ツリーの `build/` へ複製する。`build/` は一時領域で、消してよい。

置き場所は三度変えた。無視される `build/` に置いたまま、次に作業ツリーの外、次にGit LFS（3.8GiB）。どれも所有者に指摘された：無視される場所にあるものは消してよいものに見え、リポジトリの外にあるものは何のものか分からなくなり、数GiBのデータをLFSに入れるとリポジトリが重くなる。大きい物の置き場所は、作る前に所有者に確認する。

## 伏せ字と、リポジトリに入れなかったもの

R&Dの記録は、作られた機械のディレクトリと利用者名を、ログ、記録、一部のスクリプトに書き込んでいた。リポジトリに入れる前に、`<repo>` と `<home>` に置き換えた（[tools/rd/redact.py](../tools/rd/redact.py)）。置き換えたのはテキストだけで、数値配列（`.npz`）に文字列は無い（全16,119配列が数値型であることを確かめた）。

[redaction-record.json](redaction-record.json) が、変えた全てのファイルについて、変える前と後のSHA-256を持つ。

| 対象 | 件数 |
|---|---:|
| リポジトリ内のテキストの伏せ字 | 22ファイル、30か所 |
| 証拠のログの伏せ字 | 21ファイル |
| AVL由来でGPL-2.0-or-laterのファイル（Fortran 8、パッチ17）。リポジトリに入れず、証拠の束に、ライセンスの表記と一緒に入れた。製品（Apache-2.0）の一部ではない | 25ファイル |

記録が、置き換える前のハッシュで引用していたものは、再配置の道具が新しいハッシュへ再発行している。置き換える前の内容は、リポジトリにも、その履歴にも無い。

## 削除したファイル

| ファイル | 理由 |
|---|---|
| `tools/toolchains/`（3件） | 2026-10-29のLinux検証で使った、GCC 14の隔離環境の再現仕様。製品のビルドとCIは使わない |

ベースラインのタグには残っている。[foundation/build](foundation/build/) の記録 `cpp23-2026-10-29.md` は、本文でこのパスを挙げている。記録は書き換えないので、その記述は当時のままである。
