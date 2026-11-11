# 知識パッケージの検証記録

2026-11-06。製品本体の数値検証ではなく、今回作成した再利用知識・添付ツール・配置の確認。

## 構成

再利用可能な3つのSkillを`.agents/skills`へ保存した。各Skillはentrypoint、UI metadata、必要時に読む2つのreferenceで構成する。レビューSkillには標準ライブラリだけで動く証拠整合性ツールとその挙動試験を含む。

Rulesは`AGENTS.md`から共通SE-01〜14とUVLAS固有プロファイルへ接続する。知見の根拠は進め方7件・技術22件を既存文書へ対応付けた。製品採用仕様・閾値・過去の証拠は今回変更していない。

## 実行した検査

| 検査 | 結果と限界 |
|---|---|
| Skill Creator `quick_validate.py` | 3 Skill全て通過。frontmatter・命名・未完成scaffoldの構文確認であり、判断能力の証明ではない |
| agents/openai.yaml | 3件の表示名・短い説明・skill名を含むpromptを確認。自動選択を無効にするpolicyは追加していない |
| 証拠整合性ツール | 8試験通過。間接参照、重複参照、改変、欠落、不正hash／schema、root外参照、資源上限、旧artifact形式、出力上書き拒否を検査 |
| 科学的合否との分離 | 失敗した科学的resultのbyteが正しい場合はintegrityのみ通り、scientific_acceptanceはnot_evaluatedのままになることを検査 |
| UVLASの既存証拠へ適用 | implementation-readiness-closure.jsonから1,575個の固有ファイル、173 manifestを照合し、整合性成立。旧未達を含む科学的statusは保持。製品合格の判定ではない |
| ローカル参照・配布コピー | 文書／Skillの相対リンクと、個人Skillsへの配置後の内容一致を別途確認。詳細なファイルhashと結果は作業ディレクトリの記録へ保存 |

公式validatorの依存PyYAML 6.0.2は`build/skill-validation/python`へ限定して取得した。アプリ本体・既存数値環境のpackageを変更していない。証拠照合の生結果は`build/skill-validation/uvlas-evidence-integrity.json`。

## シナリオによる自己点検

独立した別エージェントによるforward testは実施していない。以下は記載した手順を実際の失敗例へ当てはめた自己点検であり、将来の挙動保証ではない。

| 入力する場面 | Skillから導く判断と、根拠の所在 |
|---|---|
| 精密な内部解は通るが返却float64の原残差が落ちる | 消費表現の未達を保持し、投影後の成功を返さない。reviewの算術・消費者検査 |
| 主候補を棄却したが、次数だけ上げて同じ未達を繰り返す | 不採用理由を変える条件を示し、識別実験へ戻る。R&Dの候補台帳／終了条件 |
| 最細2点の差が小さいので全数値1%合格とする | 漸近性・他の誤差群・交差項を要求。experimentsの極限／evidenceの合算予算 |
| 合成100k疎LUが速いので実翼探索も一晩で終わるとする | 実graph・後流・固有basis・I/O・探索を含む性能gateは未認定。design／reviewの規模観点 |
| 同じεだから二つの正則化コードは同一作用素とする | 実式・自己項・境界・版を照合。R&Dの比較条件／reviewの作用素検査 |
| 全ケースへのリンクとhashがあるので全体の開始を許可する | 意味の整合と対象依存閉包を審査し、科学的合格へ昇格させない。designの開始条件 |
| 実装に入る前に製品全試験をPythonで終えるべきとする | 決定的な方式試作と実装後の網羅／実務受入を分離。designの閉じる順序 |
| 新しい別ソルバーにもC++23／PG／BDF2／1%を強制する | 汎用Skillには固定していない。対象プロジェクトの採用仕様を取得 |
| レビュー依頼を受けて本体を勝手に書き換える | レビュー範囲を守り、仕様欠陥も含む指摘・根拠・必要検査を返す。review entrypoint／AGENTS |

## 配置と更新

個人用の配置先は`<home>/.codex/skills/solver-rd`、`solver-design-contracts`、`solver-implementation-review`。同名の既存Skillはなかった。repo内の正本から配置し、既存の無関係なSkillやglobal Rulesは変更しない。

将来更新するときはrepo正本を先に修正・検証し、配置済みコピーとの差分を確認する。インストール先の個別変更を無条件に上書きしない。実際のレビューで抜けや過剰な拘束が分かった場合に、その再現条件を基に狭く改訂する。
