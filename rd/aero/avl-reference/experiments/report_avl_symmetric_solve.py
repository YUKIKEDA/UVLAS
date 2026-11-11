"""Report completed audited solve variations without a continuum-error claim."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    r=json.loads((F/'avl-symmetric-solve-sensitivity-result.json').read_text());audit=json.loads((F/'avl-symmetric-solve-artifact-audit.json').read_text());small=json.loads((F/'avl-symmetric-solve-small-check.json').read_text())
    assert r['all_criteria_pass'] and audit['all_criteria_pass'] and small['all_criteria_pass']
    main=json.loads((F/'avl-symmetric-solve-main-runs.json').read_text());weight=json.loads((F/'avl-symmetric-solve-weight-result.json').read_text())
    response_max=max(v['response_change_ratio'] for c in r['comparisons'] for s in c['states'] for v in s['channels'].values())
    envelope_max=max(v['individual_sum_ratio'] for c in r['comparisons'] for s in c['states'] for v in s['channels'].values())
    names={'tau0':'τ帯平均','tau_p1':'τ帯P1','k0':'K帯平均','j_line':'有限辺モーメント密度','j_point':'元点モーメント密度'}
    lines=['# 対称参照の循環・荷重求解感度：実施結果','','2026-10-31。[設計](../symmetric-reference-solve-sensitivity-design.md) と [事前仕様](avl-symmetric-solve-sensitivity-experiment.md) に従い、小系2件の検証後に一様／tip-rich／tip-poorの両翼6系を診断した。元の行列生成・LU・循環・速度・力を変更せず、別の計測版で局所的な求解変化を評価した。','',
        '## 1. 荷重写像の検証','',
        '[重み検査](avl-symmetric-solve-weight-result.json) は100行の全翼観測を半翼へ写し、保存q₀・製造循環・直接積分・負側τ一次係数の符号・転置を確認した。100行はτ平均20、τ一次20、K平均20、有限辺／元点モーメント密度各20である。5チャンネルの重みを保持し、異なる単位の100成分を一つの物理ノルムにはしていない。','',
        '転置の検査は循環と保存観測の双対写像に対する代数的な検査であり、ANCFの構造荷重伝達における仮想仕事を検証したものではない。','',
        '初回の製造例は作用点をx方向だけに動かして支持線外となり、既存の支持検査で拒否された。[初回記録](avl-symmetric-solve-weight-initial-result.json) とそのソースを保持し、支持線上の鏡映点へ修正した。実流れの入力・許容差は変更していない。','',
        f'修正後の8系は{weight["seconds"]:.2f}秒・保存{weight["output_bytes"]:,} bytes。独立のDecimal80列積分は56本の半翼列を照合し、最大尺度付き差は{audit["independent_decimal_column_error"]:.6e}だった。','',
        '## 2. 小系行列と元状態の不変','',
        '[小系検査](avl-symmetric-solve-small-check.json) は両翼120未知数で、局所元行列とUH binaryがbyte一致した。NumPyの通常／転置解、反復改良の漸化式、Decimal80残差と100観測も条件内だった。','',
        '| 小系 | 通常解の尺度付き差 | 転置解の尺度付き差 | Decimal80対REAL16残差の最大絶対差 |','|---|---:|---:|---:|']
    for s in small['records']:lines.append(f'| {s["case"]} | {s["direct_error"]:.6e} | {s["transpose_error"]:.6e} | {s["decimal_quad_absolute_difference"]:.6e} |')
    lines+=['',
        '活性LU・ピボット・GAMのスナップショットを保持し、診断後にbitwise比較した。全8件で元状態が不変だった。UQ前後と親実行体のFT／FS／FE・UHもbyte一致し、after ZIP内部の元ファイルも監査した。REAL16は残差の別計算に用い、LUや解そのものをREAL16へ交換したわけではない。','',
        '## 3. 主対象の資源と残差','',
        f'[主実行記録](avl-symmetric-solve-main-runs.json) は6件・各7680未知数。合計{main["seconds"]:.2f}秒、ピーク{max(s["peak_resident_bytes"] for s in main["runs"]):,} bytes、保存{main["total_bytes"]:,} bytesで、各300秒・計2400秒・5GiB・128MiBの条件内だった。元行列とLUスナップショットの2局所配列を追加した状態で測定した。','',
        '| 翼・配分 | 元q₀の残差∞ | 1回改良後 | 2回改良後 |','|---|---:|---:|---:|']
    for s in r['systems']:lines.append(f'| {s["case"]}・{s["allocation"]} | {s["residual_infinity"][0]:.6e} | {s["residual_infinity"][2]:.6e} | {s["residual_infinity"][3]:.6e} |')
    lines+=['',
        '表は一次右辺と局所A_sによる補償和残差で、迎角1度のUH速度残差とは正規化・集計順が異なる。直接解は同じLUで右辺を解き直した別状態、反復改良は元q₀から開始した。大系のREAL16検査は固定8行×4状態だけであり、全行の高精度照合とは扱わない。','',
        f'後処理は{r["seconds"]:.2f}秒・ピーク{r["peak_bytes"]:,} bytes。100行の元観測再現、差への直接写像と大きな観測値の引き算の一致、選択行の独立残差、配分応答の写像整合性など{len(r["checks"])}処理条件が成立した。','',
        '## 4. 元の配分応答に対する求解変更の感度','',
        '各状態の循環差へW_sを直接適用してから、一様・配分変更の両系の数値変化を差し引いた。表はその差のノルム／元の配分応答ノルムの、直接解・1回・2回改良の3状態での最大値である。割合は無次元の比で、%表示ではない。','',
        '| 翼・配分 | τ帯平均 | τ帯P1 | K帯平均 | 有限辺モーメント | 元点モーメント |','|---|---:|---:|---:|---:|---:|']
    for c in r['comparisons']:lines.append('| '+c['case']+'・'+c['allocation']+' | '+' | '.join(f'{max(s["channels"][k]["response_change_ratio"] for s in c["states"]):.6e}' for k in names)+' |')
    lines+=['',
        '相殺で数値感度を小さく見せないよう、各系の変化ノルムの和／元配分応答ノルムも保存した。次表は同じ3状態での最大値。これは三角不等式による今回の感度の包絡であり、未知の真の解への上界ではない。','',
        '| 翼・配分 | τ帯平均 | τ帯P1 | K帯平均 | 有限辺モーメント | 元点モーメント |','|---|---:|---:|---:|---:|---:|']
    for c in r['comparisons']:lines.append('| '+c['case']+'・'+c['allocation']+' | '+' | '.join(f'{max(s["channels"][k]["individual_sum_ratio"] for s in c["states"]):.6e}' for k in names)+' |')
    lines+=['',
        '循環差への写像は大きな観測値同士の引き算を避けるが、浮動小数点評価の誤差を0にはしない。既存配分応答と今回のq₀写像の微差、大きな観測の引き算との差、各状態の絶対変化もJSONへ併記した。','',
        '## 5. 判定と残る作業','',
        f'今回のM12配分応答を1とすると、試した求解変更による応答変化の最大比は{response_max:.6e}、各系の変化ノルム和による最大比も{envelope_max:.6e}だった。いずれも元配分応答の100万分の1未満で、今回試した求解方法の変更では、その配分応答を説明できなかった。係数生成誤差や未検証の数値要因まで否定した結論ではない。','',
        '[成果物監査](avl-symmetric-solve-artifact-audit.json) は重み・幾何・ソース・再利用オブジェクト・ビルド・原本回帰・状態不変・ZIP・資源・数値感度を検証した。大系の格納作用素を独立に全係数検証したわけではなく、2回の改良の変化を前進誤差や係数生成誤差の上限へ変換しない。','',
        '内側と翼端を同時に変更する交絡は残っている。次は [2因子の独立対照](avl-allocation-factorial-experiment.md) として内側だけ／翼端だけを変更する両翼4入力を [準備した](avl-allocation-factorial-preparation.json)。新規6144／9216未知数は別資源仕様に従って次の作業で実行する。今回その4計算は未実施。N1280のM24尾部・合算予算は未評価、参照全体の採用・UVLAS比較・J1・製品数値コアの開始は未完了のままである。','']
    path=F/'avl-symmetric-solve-sensitivity-results.md';path.write_text('\n'.join(lines),encoding='utf-8')
    (F/'avl-symmetric-solve-report.json').write_text(json.dumps(dict(source_sha256=digest(__file__),parent_sha256=digest(F/'avl-symmetric-solve-sensitivity-result.json'),audit_sha256=digest(F/'avl-symmetric-solve-artifact-audit.json'),report_sha256=digest(path)),indent=2)+'\n')
    print('Wrote audited symmetric solve sensitivity report')
if __name__=='__main__':main()
