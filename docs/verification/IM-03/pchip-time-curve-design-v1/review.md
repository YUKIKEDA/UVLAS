# PCHIP時間場の開始審査

Issue #158、固定witness source `5b700f3`。対象はAE-07のPCHIP係数準備とcanonical有理数の時刻jetである。この部品の設計を実装へ渡せる状態と判定する。製品実装・数値認定・実験妥当性・実務受入は未実施で、全KnownFlowの開始条件が閉じたとは判定しない。

## 依存と境界

IM-00の値failure/所有/Plan規約、AE-07、FD-02/06とSI元exact値の契約を使う。FileRef snapshot #157や共有quad核には依存しない。U0とaは別の節点列を持てる。PCHIPのexact入力はここで完結するが、SiCase recipeからのadapter、constant/polynomial/expressionの他variant、sine空間評価・全微分、有限表現への投影、C2支持運動、循環・後流は既存IM03-0/IM06-0の対応段に残る。これらの下流が未実装のまま、解析の成功を返す経路は作らない。

| Issueの完了条件 | この審査の根拠 |
|---|---|
| slope・実係数・微分の一意性 | 契約2節。手計算12曲線、独立Hermite基底/積分に165 exact等値 |
| 外側・時刻・C1・side | 契約3節。最小shape/重複/逆順/範囲/hold/周期/二階sideの10拒否。負時刻・巨大周期数・C1 hold正常例 |
| 精度・所有・全call | 契約4〜6節。admission前は不変borrowの語表現走査のみ（数学的canonical性のGCDはadmission後）、checkedワード丸め、bit/node/slot/cost/保持recordを具体化。C++での実slot・allocator検査は実装gate |
| 独立参照・負対照 | 手計算slope、係数式と別のHermite基底、独立積分恒等式。正常11例121値をSciPy1.16.2へ別対照。7意味変異を検出 |
| 状態・記録・独立審査 | 本書とmanifest。旧本文/証拠を保持し、現行profileのaero行へ新契約を案内。親担当の独立レビューで指摘された丸め/検査順を反映 |

端値保持は端傾きが非零ならC1ではない。この穴をイベント指定で隠さず、要求区間で拒否する。periodicも端値だけでなく一階微分までexact一致を要求する。有限語へ先に丸めてPCHIPの分岐を選ぶ経路は認めない。

## 実行した範囲

CPython3.12.8、Windows11。主witnessは標準ライブラリFractionだけを使用した。12曲線には2点、非等間隔、零/符号、左右clip/等号、平坦端、負値、 $$ 1 $$ と $$ 1+2^{-200} $$ の異なる節点を含む。大きい共通原点 $$ 2^{900} $$ 、時刻尺度 $$ 2^{\pm200} $$ 、負周期と巨大周期数も確認した。値・一階・二階と独立積分の等値は許容差0である。

165はコードのassert文数ではなく、result.jsonのexact_equalities counterで数えた独立等値の件数である。GCDの $$ 2B+1 $$ 上界は別に1〜127の全16,129整数対で照合した。有限の全探索を一般上界の証明とはせず、契約5節の二回剰余によるbit減少の導出と区別する。

主witnessの7変異は、元時刻のdouble化、harmonicを算術平均へ変更、端clip省略、時間微分のh除算省略、負周期のfloorを切捨てへ変更、holdのC1検査省略、周期の傾き照合省略である。正常controlが完了した後に各分岐を切り替え、AssertionError/ValueErrorを実際に確認した。C++のcompile/link隔離変異やallocation故障試験としては数えない。

二次対照はSciPy1.16.2/NumPy2.5.5を検証専用uv環境へ固定した。11正常曲線のslopeとu正規化係数121値で、 $$ |a-r|/\max(1,|r|) $$ の最大は5.828670879282072e-16だった。事前の二次照合基準1e-12を満たす。この基準は製品の入力/QoI誤差予算ではない。極近時刻の1例は、binary64へ変換すると節点が潰れるためSciPyの正常対照へ含めず、exact主参照とdouble化変異で検査した。SciPyは実際に実行し、_cubic.pyのSHA256をscipy-result.jsonへ保存した。

両producerのRuff check/formatと専用strict Pyrightは成功した。rootの通常設定はdocsを除外するため、明示include・空excludeの一時configで両ファイルが実検査されたことを確認した。日本語quick lintは技術用語の反復に対するlexical-diversity情報1件で、意味を一定にするため残した。数式と固定名を言い換える修正は行っていない。

## 再現と実装後gate

新しい出力先を指定し、`python witness.py --output <new-path>`を実行する。二次対照は`uv run --python 3.12.8 --with scipy==1.16.2 --with numpy==2.5.5 python scipy_check.py --reference <primary-result> --output <new-secondary-path>`とする。主・二次とも既存出力への上書きを拒否する。別出力で再生成し、主e09f1fe0ba8e322e0c52527d90abe1cf675990dcba183a01a702ebef72792a5a、二次74831fb7e55b18e0efa1a7f20cbb4cb963f5388e5ce96b4f7490cf3450dac5d2のbyte一致を確認した。

今回の成果は設計本文と限定試作である。C++/CMake/CI・製品API・既存R&Dを変更していないため、無関係な製品suiteや全体CIは回さない。実装PRでは契約6節の全系列、全確保fault、cancel/reentry、実slot/workとcountingの一致、整数primitiveの独立対照、1/2/4thread、対象D/R・静的検査と工程/資源実測を閉じる。bit幅に対する保守的上界の実規模受入と、UVLM consumerへの有限精度投影を、この設計witnessの成功から認定しない。

採用result.jsonは過去の採用を記録した成果物として原本byteを照合した。新試作の数値入力としては使っていない。初回manifestでこれをsourcesに分類した際は、旧後流記録の未展開NPZへ再帰してMissing fileとなった。成果物をartifactsへ正しく分け、今回の限定監査が旧後流NPZ全体の再監査を含まないことを明示する。旧JSON・索引・hash・参照先の変更や、工具の検査除外は行っていない。
