# ソルバー開発知識の入口

R&Dと設計試作で得た知識を、実装者、レビュー担当、別のプロジェクトで再利用するためのもの。製品本体は実装担当が実装し、レビューを別に行う。

## Skill

正本は `.agents/skills/`。`.claude/skills/` は同じ内容の複製で、Skillを変えるcommitが両方を変える（`just skills-sync`、`just skills-check`）。

| 使う場面 | Skill | 残した手順 |
|---|---|---|
| 候補の採否が決まらない、検証が長引く | [solver-rd](../../.agents/skills/solver-rd/SKILL.md) | 仮説と反例、保存データの再利用、識別実験、未達の分類、研究の終了条件 |
| 要求は合意したが実装仕様が足りない | [solver-design-contracts](../../.agents/skills/solver-design-contracts/SKILL.md) | 観測量からの逆算、値の順方向追跡、離散式と時刻・I/F・不変条件・ケース、開始判定と引継ぎ |
| 開始審査を通ったパッケージを実装する | [solver-package-implementation](../../.agents/skills/solver-package-implementation/SKILL.md) | 決定の列挙、試作から製品へ移すものと移さないもの、独立参照と負対照、状態を分けた報告 |
| 性能目標に届かせる、高速化の採否を決める | [solver-performance](../../.agents/skills/solver-performance/SKILL.md) | 要求からの予算、条件を固定した計測、内訳、一つずつの変更、結果の同一性 |
| 製品実装・変更・完了主張をレビューする | [solver-implementation-review](../../.agents/skills/solver-implementation-review/SKILL.md) | 元系と実消費値、全必要微分・履歴、独立参照、再開・性能、仕様欠陥と実装欠陥の区別 |

各Skillは自己完結したフォルダで、必要な詳細だけreferencesから読む。特定のrepoパス、UVLM／ANCF、言語、許容差、数値依存の版を、汎用Skillの前提にしない。

## Rules

常に守る境界は `.cursor/rules/` にあり、[AGENTS.md](../../AGENTS.md) が一覧を持つ。

| 規則 | 元になったもの |
|---|---|
| `solver-engineering`、`scope`、`defer` | [SE-01〜14](solver-engineering-rules.md)、P-01、P-03、P-05、P-07 |
| `pace`、`open-ends`、`similar-findings` | P-03、P-04。2026-11-11の作業：計測基盤の行の欠落、数えずに書いた件数、形ごとに直さなかった書き換え漏れ |
| `spec-checks`、`numeric-tests`、`testing` | SE-05、SE-06、SE-10、N-01、N-02、N-05 |
| `reference-data`、`evidence` | SE-04、SE-10、SE-11、N-08、N-20、N-22。2026-11-11の再配置と、`build/` を上書きした件 |
| `api-contract`、`types`、`cpp`、`numerics` | N4の型・所有権・失敗・精度の契約、N-06、N-07、N-10〜N-14、N-18 |
| `performance`、`bench` | N4-09、N4-03、N4-04、N4-07、IO4-06、A5-06、SE-12、N-19 |
| `design-docs`、`dependencies` | N4-01、ADR 0001、P-02。仕様の文言に合わせて独自のロック形式を作った件 |
| `workflow`、`consent`、`git`、`scratch`、`no-dead-code`、`layout`、`powershell-shell` | 2026-11-11の作業での訂正：`main` へのcommit、独断での方式選択、未追跡ファイルの一括削除、部分一致の置換、役目を終えた文書の残置 |
| `review`、`review-comments` | SE-14、[UVLASレビュー・プロファイル](uvlas-review-profile.md) |

規則の形（1ファイル1関心事、BADとGOODの例、起きたことの記録）と、進め方の規則の多くは、同じ持ち主の他のリポジトリ（interpx、AuraORM、convx）の規則から移した。数値試験、計測、型の規則は、それらをこの製品とC++へ合わせて書き直したものである。

[知見と根拠](lessons-and-evidence.md)は、進め方の改善点7件（P-01〜07）と技術的知見22件（N-01〜22）を、既存の記録へ対応付けている。成功した手順だけでなく、過早な開始宣言、核心の後回し、現在地の不整合、事前仕様不足、証拠の過大解釈も、再発防止として残している。

[UVLASレビュー・プロファイル](uvlas-review-profile.md)が、採用仕様、追加規約、利用者決定への入口である。

## 維持

新しい失敗、仕様欠陥、レビュー指摘が得られたら、再現条件、根拠、適用範囲を先に記録する。そのうち将来の判断を変えるものだけをSkillとRulesへ取り込む。記述量を増やすことを改善としない。技術手順はSkill、常時の境界はRules、製品固有の数式と閾値は採用仕様、決定の理由はADRに置き、正本を重複させない。

2026-11-06時点の検証の記録は [知識パッケージの検証記録](../../rd/planning/knowledge-validation.md) にある。当時の配置は `tools/skills/` で、2026-11-11に `.agents/skills/` へ移した。構文検査や自己点検を、独立した別エージェントによる行動評価と呼ばない。`solver-package-implementation`、`solver-performance`、`.cursor/rules/` は2026-11-11に書いたもので、実際の実装作業での評価はまだない。とくに `performance` 規則は、製品のkernelがまだ無い時点で採用仕様から導いたものである。最初のkernelの実測（IM00-12、IM00-6）で、当たらない条項を直す。
