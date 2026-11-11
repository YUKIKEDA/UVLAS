# 開発の手順

エージェント向けの入口は [AGENTS.md](AGENTS.md)、規則は `.cursor/rules/` にある。ここには人が行う手順をまとめる。

## 順序

Grill（設計が分岐するとき）→ Issue → ブランチ → 作業 → PR → CI成功 → レビュー → squash merge → 振り返り。

- Issueはテンプレートから作る。白紙のIssueは作れない。Issueには [docs/roadmap.md](docs/roadmap.md) の行IDを書く。
- ブランチ名は `type/<Issue番号>-<slug>`。`main` へ直接commitしない。
- 1 Issueは1 PR。mergeは人が行い、mergeしたブランチは消す。
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

CIは同じ内容を実行し、集約ジョブ `ci` が唯一の必須チェックになる。環境の準備は [docs/build.md](docs/build.md)。

## 実装パッケージの進め方

- 実装は [開始審査](docs/plans/implementation-start-review.md) のIM-00〜08の順と依存に従う。先行パッケージの完了gateが、次の入口条件になる。
- 各PRは、実装した決定ID、独立参照、許容差、負対照、走らせた検証を書く。主張する完了状態（実装適合／数値認定／実験妥当性／性能／実務受入）を分けて示す。
- 製品の検証記録は `docs/verification/<IM-ID>/` に新しく作る。`rd/` の記録と既存の証拠は書き換えない。
- 仕様の欠陥を見つけたら、designのIssueにする。最小の反例、変えるべき式・時刻・状態・精度、影響する既存結果を書く。

## 言語

文書、Issue、PRの本文は日本語（常体）。コード、コメント、commitの件名、PRの題名は英語。

## 手元だけのファイル

`out/`（ビルド生成物）、`build/`（R&D時代の作業領域。証拠を含むので消さない）、`.dev/`（一時置き場）はgitの対象外。
