"""Audited sensitivity tables and a scientific SVG chart; no precision acceptance."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,math
from reportlab.graphics.shapes import Drawing,String,Rect
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics import renderSVG
from reportlab.lib import colors
import avl_reference_resources as provenance
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;digest=provenance.digest
NAMES={'tau0':'τ帯平均','tau_p1':'τ帯P1','k0':'K帯平均','j_line':'有限辺モーメント密度','j_point':'元点モーメント密度'}
def main():
    raw=json.loads((F/'avl-tip-redistribution-runs.json').read_text());r=json.loads((F/'avl-tip-redistribution-result.json').read_text());audit=json.loads((F/'avl-tip-redistribution-artifact-audit.json').read_text())
    assert raw['completed'] and r['methods_criteria_pass'] and audit['all_criteria_pass']
    components=json.loads((F/'avl-tip-couple-component-result.json').read_text());assert components['all_criteria_pass']
    pairs=[p for p in r['pairs'] if p['from_allocation']=='uniform'];mapping={(p['case'],p['to_allocation']):p for p in pairs}
    values=[[100*mapping[(case,allocation)]['channels'][c]['relative_to_uniform'] for case in ['FW-V01','FW-V02'] for c in NAMES] for allocation in ['tip-rich','tip-poor']]
    draw=Drawing(650,365);draw.add(String(45,341,'Fixed-count allocation sensitivity (M12, N1280)',fontName='Helvetica-Bold',fontSize=14))
    draw.add(String(60,323,'Norm of change / uniform-field norm (%); this is not an error estimate',fontSize=9))
    chart=VerticalBarChart();chart.x=65;chart.y=100;chart.width=530;chart.height=208;chart.data=values
    chart.categoryAxis.categoryNames=[f'{case} {c}' for case in ['01','02'] for c in ['tau0','tauP1','K0','Jline','Jpoint']]
    chart.categoryAxis.labels.angle=35;chart.categoryAxis.labels.boxAnchor='ne';chart.categoryAxis.labels.fontSize=8
    maximum=max(max(v) for v in values);step=10**math.floor(math.log10(maximum));chart.valueAxis.valueMin=0;chart.valueAxis.valueMax=math.ceil(maximum/step)*step;chart.valueAxis.valueStep=step
    palette=[colors.HexColor('#D55E00'),colors.HexColor('#0072B2')]
    for i,col in enumerate(palette):chart.bars[i].fillColor=col;chart.bars[i].strokeColor=colors.white
    draw.add(chart)
    for j,name in enumerate(['tip-rich','tip-poor']):draw.add(Rect(140+j*220,35,12,9,fillColor=palette[j],strokeColor=None));draw.add(String(160+j*220,35,name,fontSize=11))
    draw.add(String(65,15,'Tip and inner allocations change together; no accuracy ranking is implied.',fontSize=9))
    folder=F/'figures';folder.mkdir(exist_ok=True);figure=folder/'avl-tip-redistribution-response.svg';renderSVG.drawToFile(draw,str(figure))
    source=F/'avl-tip-redistribution-result.json';manifest=dict(source_sha256=digest(__file__),parent_sha256=digest(source),figure_sha256=digest(figure),values_percent=values,
        couple_component_sha256=digest(F/'avl-tip-couple-component-result.json'))
    (F/'avl-tip-redistribution-report.json').write_text(json.dumps(manifest,indent=2)+'\n')
    lines=['# 固定総パネル数での配分感度：実施結果','','2026-10-31。[事前仕様](avl-tip-redistribution-experiment.md) と [実行追加仕様](avl-tip-redistribution-execution-addendum.md) に従い、両翼の4ケースを実行した。M12・物理N1280・半翼640帯・7680未知数を固定し、一様配分の既存結果と比較した。','',
        '## 1. 実行と保存','',
        '[準備済み入力](avl-tip-redistribution-preparation.json) のSECTION分割だけを変更した。物理断面、弦分割、cos分点、核、基準量、同じKSMAX修正版の実行体は保持した。正側10区間の一様64に対し、tip-richは内側8区間48・外側2区間128、tip-poorは内側72・外側32とした。','',
        '| 翼 | 配分 | 時間 | ピークRSS bytes |','|---|---|---:|---:|']
    for run in raw['runs']:lines.append(f'| {run["case"]} | {run["allocation"]} | {run["seconds"]:.2f}秒 | {run["peak_resident_bytes"]:,} |')
    lines+=['',f'[実行記録](avl-tip-redistribution-runs.json) の合計は{raw["seconds"]:.2f}秒。各300秒・計1200秒・群RSS5GiB・空きとswapの条件内だった。元出力の圧縮済み保存量は{raw["total_bytes"]:,} bytes、派生場を含む合計は{r["combined_output_bytes"]:,} bytesで128MiB以内。','',
        '計測前後の8出力を各ケースで比較した。after一組はZIPに無損失保存し、全32内部ファイルを独立監査でoutとバイト一致させた。ZIP検証後に新規afterの重複実体だけを除去し、元出力と既存実験は保持した。[監視試験](avl-tip-watchdog-result.json) は境界値・拒否・実プロセス群停止を含む。','',
        '## 2. 6場の処理検査','',
        f'[後処理](avl-tip-redistribution-result.json) は基準2場と新規4場を区別した。実際の区間別帯数、640帯・7680未知数、物理全翼15360パネルを確認した。元残差は最大相対{max(m["residual_relative"] for m in r["meshes"]):.6e}。荷重／元点モーメント保存、帯の独立集計、6場対の差の上限と地域の二乗加算を含む{len(r["checks"])}処理条件が全て成立した。後処理は{r["seconds"]:.2f}秒、ピーク{r["peak_bytes"]:,} bytes。','',
        '上限は二つの保存場の差に対するもので、連続解への誤差上界ではない。以下の荷重値は迎角に関する静的一次係数で、有限迎角のFT係数の保存検査とは区別する。','',
        '## 3. 一様配分に対する応答','',
        '表は差の重み付きノルム／同じ翼の一様M12N1280のノルム。値が小さいことは真の解に近いことを意味しない。','',
        '| 翼・配分 | τ帯平均 | τ帯P1 | K帯平均 | 有限辺モーメント | 元点モーメント |','|---|---:|---:|---:|---:|---:|']
    for p in pairs:lines.append('| '+p['case']+'・'+p['to_allocation']+' | '+' | '.join(f'{100*p["channels"][c]["relative_to_uniform"]:.6g}%' for c in NAMES)+' |')
    lines+=['','![同じ総未知数での分布応答](figures/avl-tip-redistribution-response.svg)','',
        '| 翼・配分 | 揚力傾斜の符号付き差 | 元点モーメント傾斜の符号付き差 |','|---|---:|---:|']
    for p in pairs:lines.append(f'| {p["case"]}・{p["to_allocation"]} | {p["scalar_changes"]["lift"]["signed"]["total"]:+.8e} | {p["scalar_changes"]["moment_point"]["signed"]["total"]:+.8e} |')
    lines+=['','## 4. 高密度／低密度配分の方向と寄与','',
        '次の余弦はtip-rich−一様とtip-poor−一様の差の方向。一様格子を基準に反対へ応答するかを表し、真の解に向かう方向は示さない。','',
        '| 翼 | τ帯平均 | τ帯P1 | K帯平均 | 有限辺モーメント | 元点モーメント |','|---|---:|---:|---:|---:|---:|']
    for response in r['responses']:lines.append('| '+response['case']+' | '+' | '.join('未定義' if response['channels'][c]['direction_cosine'] is None else f'{response["channels"][c]["direction_cosine"]:.6f}' for c in NAMES)+' |')
    lines+=['','地域の二乗寄与は、基準と同じ根元4帯・中間12帯・翼端4帯へ分けた。次表は各配分差の二乗ノルムに占める翼端側の割合。','',
        '| 翼・配分 | τ帯平均 | K帯平均 | 元点モーメント |','|---|---:|---:|---:|']
    for p in pairs:lines.append('| '+p['case']+'・'+p['to_allocation']+' | '+' | '.join(f'{100*p["channels"][c]["regional_fractions"]["tip"]:.3f}%' for c in ['tau0','k0','j_point'])+' |')
    lines+=['','## 5. 元点モーメントの補正偶力の項','',
        '[事後の追加仕様](avl-tip-couple-diagnostic-addendum.md) で、保存済みの元点モーメント差を有限辺の項と補正偶力の項へ分解した。[成分診断](avl-tip-couple-component-result.json) は6組で再構成を確認し、最大絶対差は1e−12以内。次の値は各成分ノルム／元点差ノルムであり、交差項があるため足して1になる割合ではない。','',
        '| 翼・配分（一様との差） | 有限辺の項／元点差 | 補正偶力の項／元点差 |','|---|---:|---:|']
    for c in components['pairs']:
        if c['from_allocation']=='uniform':lines.append(f'| {c["case"]}・{c["to_allocation"]} | {c["line_to_point_norm"]:.6f} | {c["couple_to_point_norm"]:.6f} |')
    lines+=['',
        'tip-richの元点差では補正偶力の項のノルムが元点差ノルムの約0.94倍に相当し、有限辺の項は約0.064〜0.066倍だった。これは保存観測上の分解であり、補正偶力の変更には循環と作用点の両方が関係する。幾何だけの原因や真の誤差と断定しない。','',
        'tip-richでは、有限辺モーメントの帯間相殺比はFW-V01で約0.00367、FW-V02で約0.000760となった。総量差が極小でも分布差は非ゼロであり、総量だけで判断しない。','',
        '## 6. 判定と残る切分け','',
        '[独立成果物監査](avl-tip-redistribution-artifact-audit.json) は入力・実行体・原本・ZIP・派生場のハッシュ、時間・メモリ・保存量、実区間配分、Decimalでの差ノルム、方向余弦と係数差を再検証した。結果図の数値とハッシュは [報告記録](avl-tip-redistribution-report.json) に保持した。','',
        '今回の比較では翼端と内側の帯数を同時に変えている。観測された応答を翼端だけの因果効果、精度改善率、真の解への距離と呼ばない。また非相似格子のため、Richardson次数や尾部推定を返さない。M12の診断でM24の合算誤差予算を補完しない。参照全体の受入、UVLAS対参照の一致、J1と製品数値コアの開始は未完了。','',
        'tip-richのτ応答は基準ノルム比で約3.2〜3.4e−8と小さい。残差が条件内でも、それだけではこの微小応答の前進誤差上限を与えない。次は [対称縮約下の循環・荷重求解感度](../symmetric-reference-solve-sensitivity-design.md) の写像・処理・回帰と資源仕様を具体化して検査し、その後に内側と翼端を独立に変える対照を設計する。','']
    (F/'avl-tip-redistribution-results.md').write_text('\n'.join(lines),encoding='utf-8');print('Wrote audited sensitivity tables and SVG')
if __name__=='__main__':main()
