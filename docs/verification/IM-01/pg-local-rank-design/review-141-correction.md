# PR #141：anchor核fixtureの訂正

2026-11-18。レビューr4613633044、対象head c00c882e2ea39db3010277391fb7135c32861697。

指摘を再現した。元fixtureのJ[:,3]はC1の単位成分+1を−1へ上書きし、Csを加えていたため、C成分がs−1となり背景anchor s=0で−1だった。AJ=BZは任意の原列結合でも成立するので、この誤りを検出できなかった。正しい持上げはC1−Csであり、J[:,3]=C1−(C1−Cs)=Csへ修正した。

同型のJ全8列を点検した。R側は元anchor値1/0/0/0/0/0/2/1、P側は0/0/1/0/0/0/0/0、C側は0/0/0/1/0/0/0/0。全t半cellでこの手書きliteralを照合し、C以外の減算も正しいことを確認した。fixtureはtに関して定数なのでt係数は0次だけであり、全field・全t半cell・全J列の全係数を調べる。一般の製品辞書の全t次数を実行した主張ではない。

修正版は全背景係数から局所座標u=−1を直接代入し、各J列のtraceを別集計する。n=1/2/4×活動領域3種で合計1,008件をexact零判定した。活動領域がanchor外でも、活動行ではなく全背景から検査する。

負対照は従来6種に2種を追加した。旧C誤代入を復元すると、全t半cellのC/J3が−1となり検査が失敗する。持上げの減算を加算へ変えたJbad=2I−Jでは、各traceが手書き元anchor値の2倍となり検査が失敗する。両方について失敗を捕捉し、失敗しない実装ならproducer全体を失敗させる。9行列の独立Fraction Gaussとの選択・全再構成・尺度復元は許容差0で一致し、全8負対照を検出した。

旧producer・旧出力・旧manifestのbytesをreview-141-originalへ保存した。旧design-result.jsonとverification.mdも保持し、その6変異だけではJのanchor核所属を確認できなかったことを本記録で訂正する。現manifestは修正producerとdesign-result-v2.jsonを指す。旧manifestの.json.txtは保存資料であり、現在のファイル配置への照合用ではない。

数学上のJ定義と資源不採用判断は不変。製品rank実装・実務性能・全PG辞書の認定を主張しない。製品コード変更がないため、C++の無関係suiteは実行しない。

対象検査：trace 208 IDs・36 with tests・42 with records・0 problems。日本語lintは文長の均一性についてwarn1件。定義・再現・正常・負対照・保存範囲を段落ごとに分けて独立通読し、数値や範囲を曖昧にする書換えは行わない。

独立レビュー担当が修正版producerを別worktreeで再実行し、全9例・全8負対照の完走を確認した。出力SHA-256は `7e63ae0872fc220a6284fb0cc3d2416f371ded18b82742d49de9f22601c9cd7e` でdesign-result-v2.jsonとbyte一致した。
