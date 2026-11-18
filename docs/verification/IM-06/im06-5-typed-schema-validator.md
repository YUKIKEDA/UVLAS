# 型付きcase schema validatorの検証範囲

Issue [#126](https://github.com/YUKIKEDA/UVLAS/issues/126)、IO4-01、既存IM06-0の全schemaを対象とする。実行結果とsource commitは同名のJSON記録にまとめる。本書だけを合格記録として扱わない。

`validate_case_schema` は、所有する入力文書から、全fieldを検査した不変の `SchemaCase` を作る。191種類のnamed record、閉じたvariantとtag、型付き参照、exactなNumber/Rational、次元メタデータを保持する。元のsource bytes、木、配列索引、正規digitsも所有する。公開fieldは同じownerに結び付き、利用者が作ったordinalを別caseへ渡して読める経路を設けない。

検査は、全treeのkey・kind・shape、exact scalar、参照と相互条件、typed record構築の順に行う。後方の未知keyは前方の値不正より先に拒否する。配列は所有する平坦な索引とoffset/countで読み、一つの要素へのアクセスで先頭からedgeをたどらない。payloadとscratchの全call計画、累積課金、work、桁数と指数の上限を検査し、処理前にworkを課金する。字句長はDecimalTokenBytes、仮数と一時演算の桁数はExactDigits、指数値はExponentMagnitudeで個別に制限する。

## 正本の例を製品へ接続する

凍結した144 fragmentと5 root、完全物理初期場の4正常例と22負例は、producerがnumber tokenを変えずにC++へ取り込む。編集用の試験helperも、JSON numberをbinary64へ変換しない。初期追補の負例は、schemaで拒否する13例と、schemaを通してprepareで拒否すべき9例を独立に区別する。

fragmentには同じ入力文書内の明示contextを付け、公開Caseと共通の型・shape・exact・参照・相互条件のkernelを使う。内部の `SchemaPayload` は公開 `SchemaCase` とは別型であり、Case全体の成功資格を与えない。公開入口は全nodeが型付け済みであることも確認する。191 recordと806 memberの到達を数え、全required fieldの欠落、未知key、null kind、固定shape、閉じた全tag集合を検査する。

旧 `Perturbation.physical` の点値だけの例は追補により失効したため拒否する。原例のbytesは保持する。旧 `Geometry.gmsh` 例は二つの異なるfileに同じSourceIdを指定しているため、原文を重複IDとして拒否する。reference fileのIDだけを3へ変えた派生正常例を別に置き、製品のSource整合検査を緩めない。

## 数値と失敗を独立に照合する

理論式の評価には進まず、Numberの原字句、Rationalの交差積、時間区間の整数商、誤差群の和と被覆、DAGの順序・到達・次元を照合する。INT64_MINへ正規化されるdecimal exponentも検査し、巨大な10の冪を生成しない。シェルの初期物理場はR/P/C、明示されたhalf-cell情報、梁と剛体の全状態、PointCheckの片側tagを保持する。

資源試験では元文書の件数と最長pathを独立に数え、全host quotaとcase宣言のmin、計画P/Sと累積料金、workの境界を確かめる。全実確保点の故障を注入し、入力と以前の成功結果の全公開payloadが保たれることを検査する。diagnosticの有無、path確保失敗、移動、元文書・plan破棄、独立した並列呼出しも対象とする。全payloadのvisitorはrecord宣言から生成するが、期待する数値や次元のoracleには使わない。

## 成功が意味する範囲

この成功はschemaの資格である。SI変換、式の数値評価と微分、外部meshの内容、活動幾何とspan所属、SO(3)、拘束整合、独立係数への所属、最終q/v、実行とcheckpointは後続のprepare等で認定する。22負例のprepare側9例をschema成功に数えても、その物理状態が適合したとは主張しない。CLI `validate` の完成や、Windows以外の実行成功も、この記録から推定しない。
