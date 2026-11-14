# 開発の手順

エージェント向けの入口は [AGENTS.md](AGENTS.md)、規則は `.cursor/rules/` にある。ここには人が行う手順をまとめる。

## 順序

Grill（設計が分岐するとき）→ Issue → ブランチ → 作業 → PR → レビュー → squash merge → 振り返り。CIはPR上で非同期に実行し、mergeに必要な成功はmerge前に確認する。

- Issueはテンプレートから作る。白紙のIssueは作れない。Issueには [docs/roadmap.md](docs/roadmap.md) の行IDを書く。
- ブランチ名は `type/<Issue番号>-<slug>`。`main` へ直接commitしない。
- 1 Issueは1 PR。mergeは人が行い、mergeしたブランチは消す。
- エージェントへ実行を依頼したタスクは、検査・commit・通常のpush・PR作成まで進める。commit・pushごとの確認は不要で、変更内容は原則PRレビューで確認する。設計や範囲の判断が分かれる場合は、PR前に確認する。
- エージェントはCIの完了を原則待たず、PRのURLと現在の状態を報告して返す。完了までのポーリングは行わない。利用者が待機を指定した場合や、結果が次の処理の必須入力となる場合だけ待つ。
- commitの件名は `type(scope): subject`（英語）。

| テンプレート | ラベル | commitのtype |
|---|---|---|
| feat | `type:feat` | `feat` |
| bug | `type:bug` | `fix` |
| refactor | `type:refactor` | `refactor` |
| test | `type:test` | `test` |
| design | `type:docs` | `docs` |
| task | `type:chore` | `chore` `build` `ci` |
| spike | `type:spike` | `chore` |

## 振り返り

merge後、PRに3点を書く：得た教訓、それを書き込むファイル、追加のリファクタリングの要否。教訓があれば `chore` のIssueを立て、規則・Skill・仕様・ADRへ書く。

## ゲート

差分と依存先に影響する検査を選ぶ。下表を毎回全て実行する必要はない。作業中は局所試験を使い、全体検証が必要な場合は、文書・規約を含む編集を終えた最後に1回へまとめる。通過済みの対象は、影響する変更か失敗の修正がある場合だけ再実行する。手元とCIで同じ環境・対象を重複検証しない。

| コマンド | 内容 | ゲート |
|---|---|---|
| `just gate` | 構成、ビルド、試験（Windows：Debug／Release。Linux：Debug／Release／sanitizer） | はい |
| `just format-check` | clang-formatの差分がないこと | はい |
| `just lint` | clang-tidyの指摘がないこと（構成済みであること） | はい |
| `just lint-python` | 道具のPython：ruff（検査と整形）、pyright strict | はい |
| `just trace-check` | 試験と記録が名指しするIDが全て定義されていること。追跡表が最新であること | はい |
| `just skills-check` | `.claude/skills` が `.agents/skills` の複製であること | はい |
| `just rd-check` | R&Dツリーのリンク、コンパイル、ハッシュ対 | はい |
| `just rd-verify` | 整理後のR&Dツリーがベースラインの変換と一致すること（遅い） | いいえ |
| `just format`、`just skills-sync` | 整形、Skillの複製 | いいえ |

CIは変更範囲からジョブを選び、集約ジョブ `ci` が必要な検査の成功を確認する。影響しないジョブの予定された省略は許容する。共有部分や影響を限定できない変更は対象を広げる。環境の準備は [docs/build.md](docs/build.md)。

自動CIはPRで実行する。merge後に同じ変更を再度全体検証するpush起動は行わない。全体を明示的に確認する場合は `workflow_dispatch` を使う。独立したtoolchain試験だけの変更は対応するビルド対象・CTest・静的解析に絞り、CMake・共有ソース・CI自体の変更や未知のパスは広く検証する。

## 並列に進める

着手前にタスクの先行条件・成果物・編集ファイル・共有状態を整理する。依存する作業は順に行い、結果待ちや同時編集の競合がない作業は並列に進める。CI・レビュー待ちには独立した許可済みの作業を進める。編集中の同じファイルをビルド・解析させず、検証対象は固定する。

## 実装パッケージの進め方

- 実装は [開始審査](docs/plans/implementation-start-review.md) のIM-00〜08の順と依存に従う。先行パッケージの完了gateが、次の入口条件になる。
- 各PRは、実装した決定ID、独立参照、許容差、負対照、走らせた検証を書く。主張する完了状態（実装適合／数値認定／実験妥当性／性能／実務受入）を分けて示す。
- 製品の検証記録は `docs/verification/<IM-ID>/` に新しく作る。`rd/` の記録と既存の証拠は書き換えない。
- 仕様の欠陥を見つけたら、designのIssueにする。最小の反例、変えるべき式・時刻・状態・精度、影響する既存結果を書く。

## 言語

文書、Issue、PRの本文は日本語（常体）。コード、コメント、commitの件名、PRの題名は英語。

## 手元だけのファイル

`out/`（ビルド生成物）、`build/`（R&D時代の作業領域。証拠を含むので消さない）、`.dev/`（一時置き場）はgitの対象外。
