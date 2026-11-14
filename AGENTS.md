# UVLAS：エージェントの入口

実務で利用できる品質のUVLM–ANCFシェル連成ソルバー（C++23、CLIとPython）。実務範囲・品質を、学習目的や小規模デモへ置き換えない。このファイルが入口。人向けの手順は [CONTRIBUTING.md](CONTRIBUTING.md)、守るべき規則は `.cursor/rules/` にある。

## 正本

| 知りたいこと | 場所 |
|---|---|
| 作業の順序と現在地 | [docs/roadmap.md](docs/roadmap.md) |
| 製品要求、受入精度、受入計画 | [docs/requirements/](docs/requirements/) |
| 採用仕様（構造・空力・連成・数値基盤・入出力） | [docs/specs/](docs/specs/)、優先関係は [レビュー・プロファイル](docs/engineering/uvlas-review-profile.md) |
| 実装パッケージIM-00〜08と残リスクR-01〜07 | [docs/plans/implementation-start-review.md](docs/plans/implementation-start-review.md) |
| ディレクトリ構成とモジュール | [docs/conventions.md](docs/conventions.md) |
| 決定の理由 | [docs/adr/](docs/adr/) |
| 品質保証：追跡、基準結果、結果の変更管理、版の意味 | [docs/quality/quality-plan.md](docs/quality/quality-plan.md)。あるIDの検証の所在は `just trace <ID>` |
| R&Dの経緯と証拠（採用・不採用） | [rd/README.md](rd/README.md) |
| ビルド、試験、ゲートのコマンド | [justfile](justfile)、[docs/build.md](docs/build.md) |

現在の行のIDは `docs/roadmap.md` だけが持つ。`.dev/` は一時置き場で、そこにあるものは決定ではない。

## 役割

本体の実装は実装担当が行い、レビューを別に行う。各依頼の対象と役割に従う。レビュー依頼では仕様自身も審査し、修正実装の依頼がなければ指摘・証拠・採否判断を返す。

## 必ず読む規則

作業に合う規則を読む。`alwaysApply` の規則は全ての変更に適用する。

| 規則 | 一行 |
|---|---|
| `workflow` | Grill、Issue、ブランチ、PR、人のmerge、振り返り。ゲート。commit前にブランチを確かめる |
| `consent` | 質問は決定ではない。何を利用者が決め、何を自分で決めるか |
| `git` | 破壊的な操作は名指しの指示で。commitとpushは頼まれたときだけ。ステージを読む |
| `conventional-commits` | `type(scope): subject`。`perf` は使わない |
| `pull-requests` | テンプレートの見出し。Verificationは完了条件の各行から書く |
| `github-closing` | 番号の前の閉じる動詞は、否定文でもIssueを閉じる |
| `review-comments` | 指摘は、それが書かれたheadについてのもの |
| `scope` | 範囲は要求が決める。半端な能力を出さない。状態を分けて言う |
| `defer` | 先送りは、名前の付いたロードマップの行だけ |
| `pace` | 行の外のものは切り出す。調査の限り。次の作業を決める。長い実行は複製の上で |
| `open-ends` | 未了の作業や、説明できない数値を残したままPRを出さない |
| `spec-checks` | 試験は、仕様の文の全体を実装する |
| `similar-findings` | 同じ種類の問題を、全体を検索して直す |
| `record-corrections` | 訂正は、同じ変更の中で規則・Skill・仕様・ADRへ書く |
| `scratch` | 一時ファイルを残さない。未追跡のファイルをまとめて消さない |
| `no-dead-code` | 取り除いた経路だけに仕えていたものを全て取り除く |
| `layout` | 置き場所。必要になる前に作らない。役目を終えたものを残さない |
| `design-docs` | 仕様の置き場所と変え方。他のコードについての主張。現在の最良を判断する |
| `dev-docs` | 図はMermaid、数式はLaTeX |
| `solver-engineering` | SE-01〜14とUVLAS固有の回帰観点への入口 |
| `quality` | 試験にIDのタグ。結果を変える変更には結果変更記録。来歴と環境の認定 |
| `publish` | 公開できる状態、順序、誰が行うか |
| `evidence` | R&Dの証拠を書き換えない。再配置は道具を通す。消してはいけないものを無視される場所に置かない |
| `reference-data` | 製品を検証する参照の独立性、R&Dの記録の使い方、不確かさ |
| `api-contract` | 所有する結果、読取り専用の親、値としての失敗、例外を出さない、決定的 |
| `types` | 不正な状態を表現できなくする |
| `cpp` | 未定義動作、検査付きの寸法、警告とlint、命名 |
| `numerics` | 非有限、比較、安定な形、残差の尺度、精度の保持と投影 |
| `numeric-tests` | 独立な参照、固定の許容差、変異の確認、性質の試験 |
| `testing` | Catch2、層とラベル、構成、系列の全成員を試す |
| `performance` | 費用の数え方、ビット一致を保つ並列、確保、Eigen、高速化の採否 |
| `bench` | 測り方、比べ方、数値の書き方 |
| `python` | `tools/` のPythonはruffとpyright strictを通す。除外は根に固定する |
| `dependencies` | vcpkg manifestで版を固定。足す前の確認 |
| `review` | レビューする側の範囲と指摘の書き方 |
| `language` | 文書は日本語、コードとcommit件名は英語 |
| `powershell-shell` | この開発機のシェル。スクリプトでファイルを直すときの決まり |

## Skills

正本は `.agents/skills/`。`.claude/skills/` は同じ内容の複製で、Skillを変えるcommitが両方を変える（`just skills-sync`）。`skills-lock.json` が、外から来たSkillを固定する。それらは、上流の `LICENSE` を自分のフォルダに持つ。

- `grilling`：設計が分岐するとき、Issueの前に。
- `github-closing`：`gh pr create` の前に、PRが終えるものだけを閉じることを確かめる。
- `natural-japanese`：日本語の文書を書き、確かめる。書いた本人ではない者が読み直す。
- `solver-rd`：方式の選定、未達原因の切り分け、試作の採否。
- `solver-design-contracts`：要求を、式・離散残差・I/F・不変条件・検証計画へ落とす。開始判定。
- `solver-package-implementation`：開始審査を通ったパッケージを、製品コードと独立試験へ移す。
- `solver-performance`：要求から予算を作り、計測と内訳に基づいて改善し、結果の同一性を確かめる。
- `solver-implementation-review`：実装・変更・完了主張のレビュー。証拠ハッシュの整合検査ツールを含む。

参照データを作り直す手順、プロファイルを採って読む手順は、まだSkillにしていない。その流れが製品に無いためで、最初のfixtureと最初のkernelを足すIssue（IM00-11、IM00-12）で書く。

## 環境

Windows 11（PowerShell 7、Visual Studio 2026のMSVC）とLinux（GCC 14）。CMake 3.25以上、Ninja、vcpkg（`VCPKG_ROOT`）。Windowsのコマンドは `pwsh tools/dev/msvc.ps1 <command>` を通すか、`just` のレシピを使う。
