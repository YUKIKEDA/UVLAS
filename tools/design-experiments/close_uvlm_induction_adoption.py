"""2.1 closure with immutable structural parent and individual-column gates."""
from pathlib import Path
import hashlib,json,re,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(c,m):
    if not c:raise AssertionError(m)
def main():
    start=time.monotonic();dest=ROOT/'docs/verification/uvlm-induction-adoption-disposition.json';require(not dest.exists(),'immutable destination');sources={};artifacts={};parentpath=ROOT/'docs/verification/structural-adoption-closure.json';parent=json.loads(parentpath.read_text());a=json.loads((ROOT/'docs/verification/structural-adoption-audit.json').read_text());require(parent['passed'] and parent['structural_review_complete'] and a['passed'],'structural prerequisite')
    for group in ['sources','artifacts']:
        for name,h in a[group].items():require(sha(ROOT/name)==h,'changed structural evidence '+name)
        (sources if group=='sources' else artifacts).update(a[group])
    for group in ['immutable','status_sync_source']:
        for name,h in parent[group].items():require(sha(ROOT/name)==h,'changed structural closure '+name);sources[name]=h
    sources[parentpath.relative_to(ROOT).as_posix()]=sha(parentpath)
    resultpath=ROOT/'docs/verification/uvlm-induction-adoption-v2-result.json';columns=ROOT/'docs/verification/uvlm-induction-column-result.json';result=json.loads(resultpath.read_text());col=json.loads(columns.read_text());require(result['passed'] and col['passed'] and not result['failures'] and not col['failures'],'selected induction evidence');require(len(result['records'])==8 and len(col['records'])==2,'coverage');records={r['case']:r for r in result['records']};individual=col['records'][0];require(len(individual['velocity_columns'])==81 and len(individual['geometric_derivative_columns'])==117,'individual coverage')
    for r in individual['velocity_columns']+individual['geometric_derivative_columns']:require(r['absolute']<=1e-10 or r['relative']<=1e-8,'individual column accuracy')
    require(col['records'][1]['previous_returned_unqualified_result'] and col['records'][1]['current_rejection']=='far_induction_precision_unresolved','nonfinite failure regression')
    near=records['general_triangle_near_and_declared_self_each_P1_basis']['cases'];require(len(near)==15 and len([r for r in near if r['distance']==0])==6,'near and self');require(all(r['relative_error']<=1e-8 or r['absolute_error']<=1e-10 for r in near),'near accuracy');require(len(records['classification_and_failure_guards']['rejected'])==10 and records['classification_and_failure_guards']['input_unchanged'],'failure guards')
    shared=records['shared_edge_spanwise_traces_and_incident_mesh_node'];require(max(shared[k] for k in ['aggregated_velocity_error','shared_vertex_chain_rule_error','incident_node_velocity_error'])<1e-8 and shared['trace_averaging_negative_control']>.01,'shared trace and geometry');require(records['constant_ring_and_nonzero_whole_self']['whole_ring_self_norm']>1e-3,'self retained')
    for r in records['rigid_scale_orientation_covariance']['cases']:require(max(r['velocity_error'],r['Jacobian_error'])<1e-8,'covariance')
    queue=[resultpath,columns];seen=set()
    while queue:
        p=queue.pop()
        if p in seen:continue
        seen.add(p);r=json.loads(p.read_text());sources[p.relative_to(ROOT).as_posix()]=sha(p)
        for name,h in r.get('sources',{}).items():
            q=ROOT/name;require(sha(q)==h,'changed induction source '+name);sources[name]=h
            if q.suffix=='.json':queue.append(q)
        if isinstance(r.get('artifact'),dict):
            ar=r['artifact'];require(sha(ROOT/ar['artifact'])==ar['artifact_sha256'],'changed operator snapshot');artifacts[ar['artifact']]=ar['artifact_sha256']
    # Older experiments use a different provenance schema. Preserve them and
    # their recorded status; do not turn an old failed experiment into a pass.
    legacy={}
    for stem,runner in [('uvlm-column-update','check_uvlm_column_update.py'),('uvlm-adaptive-quadrature','check_uvlm_adaptive_quadrature.py'),('uvlm-wp1-three-dimensional','check_uvlm_wp1_three_dimensional.py')]:
        p=ROOT/'docs/verification'/f'{stem}-result.json';r=json.loads(p.read_text());require(sha(ROOT/'tools/design-experiments'/runner)==r['source_sha256'],'legacy source '+runner);sources[p.relative_to(ROOT).as_posix()]=sha(p);sources['tools/design-experiments/'+runner]=r['source_sha256'];legacy[stem]=r['all_criteria_pass']
        for name,h in r['dependency_sha256'].items():q=ROOT/'tools/design-experiments'/name;require(sha(q)==h,'legacy dependency '+name);sources[q.relative_to(ROOT).as_posix()]=h
    require(legacy['uvlm-column-update'] and legacy['uvlm-adaptive-quadrature'],'qualified limited legacy results')
    for name in ['docs/verification/avl-reference-diagnostic-disposition.md','docs/verification/avl-reference-diagnostic-disposition.json']:
        sources[name]=sha(ROOT/name)
    for name in ['uvlm-induction-exploration.log','uvlm-induction-qualified.log','uvlm-induction-qualified-v2.log','uvlm-induction-qualified-v3.log']:
        p=ROOT/'build/analysis'/name;require(p.exists(),'diagnostic log');artifacts[p.relative_to(ROOT).as_posix()]=sha(p)
    stored=np.load(ROOT/result['artifact']['artifact']);require(stored['basis'].shape==(3,3) and stored['Jacobian'].shape==(3,3,13),'operator dimensions');require(np.isfinite(stored['basis']).all() and np.isfinite(stored['Jacobian']).all(),'operator finite')
    # Snapshot eps=.02 belongs to the local triangle; the shared fixture uses
    # eps=.08, recorded explicitly here rather than inferred from that scalar.
    context=dict(local_triangle_eps=.02,shared_mesh_eps=.08,mesh_Jacobian_order='point3, shared_vertices27, eps1')
    spec=ROOT/'docs/uvlm-induction-adopted-specification.md';text=spec.read_text(encoding='utf-8');decisions=[f'UI-D{i:02d}' for i in range(1,7)];require(all(x in text for x in decisions),'adopted decisions');metrics=dict(verification_groups=10,individually_qualified_velocity_columns=81,individually_qualified_geometric_columns=117,max_individual_velocity_relative=max(r['relative'] for r in individual['velocity_columns']),max_individual_geometry_relative=max(r['relative'] for r in individual['geometric_derivative_columns']),max_near_relative=max(r['relative_error'] for r in near),shared_vertex_chain_rule_error=shared['shared_vertex_chain_rule_error'],trace_averaging_defect=shared['trace_averaging_negative_control'],incident_mesh_node_error=shared['incident_node_velocity_error'],negative_controls=10)
    decision=dict(passed=True,uvlm_2_1_complete=True,uvlm_review_complete=False,structural_review_complete=True,product_implementation_start=False,decisions=decisions,metrics=metrics,legacy_recorded_status=legacy,next_work='2.2 wake assembly state topology and loads',specification=spec.relative_to(ROOT).as_posix(),default_core_radius_selected=False,finite_wing_reference_qualified=False,free_wake_time_update_adopted=False)
    dest.write_text(json.dumps(decision,indent=2),encoding='utf-8');report=ROOT/'docs/verification/uvlm-induction-adoption-results.md';report.write_text(f'''# UVLM 2.1：有限3D誘導核の採用

2026-11-05。**指定した核・半径のW-P1誘導作用素と局所幾何微分を採用し、2.1を完了した。** 次は2.2の後流・組立・状態・荷重。UVLM全体・既定半径・自由後流更新・有限翼荷重の採用は未完了。[仕様v1.0](../uvlm-induction-adopted-specification.md)、[判定JSON](uvlm-induction-adoption-disposition.json)、[監査](uvlm-induction-adoption-audit.json)。

## 成立した数値比較

| 比較 | 結果 |
|---|---|
| 平面・skew・非平面、3コア条件 | 18配置のW-P1両列を、別の正則化面積分へ照合 |
| 一般3D三角形の面内／辺／頂点への接近 | 理想核の面外距離1e−2/1e−3/1e−5と、ε=.02/.08の面上自己を含む15配置。各P1基底を検証、最大相対差 {metrics['max_near_relative']:.6g} |
| 点・頂点・半径のJacobian | 全13変数を別の面積分のcomplex-step微分へ照合。incident点の物質写像、共有頂点の全出現箇所の鎖則も一致 |
| 弱い基底の独立判定 | [追加監査](uvlm-induction-column-result.json)で速度81列・微分117列を個別に判定。相対1e−8または絶対1e−10を全列が満たした |
| 自己・共有・帯トレース | 一定μで独立な解析渦輪へ戻り、渦輪全体の自己速度は非零。共有頂点鎖則の差 {metrics['shared_vertex_chain_rule_error']:.6g}、帯トレースを平均する負対照の差 {metrics['trace_averaging_defect']:.2%} |
| 共変性と評価分岐 | 速度／Jacobianの回転・並進・尺度1e−3/1e3、向き反転、同一平面外側の極限、距離8/1e3/1e6の遠方で成立 |
| 失敗処理 | 非incident接触、偽incidence、理想核の面上自己、負／極小コア、退化、非有限点、遠方の非有限微分、共有辺の半径不整合の10負対照を拒否。入力不変 |

[主8群](uvlm-induction-adoption-v2-result.json)と追加2群の計10群。近接の境界線と一定面内渦層を解析積分し、遠方は同じ面積分のGauss12/8で速度と微分を照合した。格子・核・正則化項・精度条件を変えて近接誤差を受け入れたものではない。

## 不合格と修正の扱い

最初の近接面積分参照は32次→48次の差が約6.46e−7で、参照自身の条件を満たさなかった。投影点を頂点にしたfanと対数状の半径分割、局所座標でのR評価・補償和を使い、48次／64次の双方を照合して参照を適格化した。物理式や許容差は変えていない。向き反転検査の差を二重に取る誤記も、検査式だけを訂正した。

[v1](uvlm-induction-adoption-result.json)はその検査範囲で合格したが、極端に遠い有限入力では非有限微分を返す経路が残った。[v2](uvlm-induction-adoption-v2-result.json)は遠方の全値・全微分の有限性確認を追加し、その入力を明示的に拒否する。両ソース・結果と[回帰検査](uvlm-induction-column-result.json)を保存した。非有限値を誤差推定のNaN比較で成功扱いしない。

## 維持する既存判断と次作業

幅一列・連続指定後流の局所平均圧力は、既存の4096分割で成立した範囲を再利用する。既存3D固定求積の初期結果は未達を含むため、今回の採用で旧結果を合格へ書き換えていない。AVLの参照診断終了・精度受入への不採用も維持し、その不足は2.3に残す。

2.2ではこの局所核を、安定IDによる自己／共有・全域交差診断、後流の非貫通／放出／移流、履歴・試行／確定、有限幅の荷重と接線へ統合する。2.3で有限翼・後退テーパ・複数面・非定常の独立参照と誤差予算、2.4でAE系の全体採用を閉じる。今回の正則化自己評価を実渦核の物理妥当性や自由後流の時間安定性の合格とはしない。
''',encoding='utf-8')
    for p in [spec,dest,report,Path(__file__)]:sources[p.relative_to(ROOT).as_posix()]=sha(p)
    pending=[ROOT/p for p in sources if p.endswith('.py')];visited=set()
    while pending:
        p=pending.pop()
        if p in visited:continue
        visited.add(p)
        for module in re.findall(r'^(?:from|import)\s+([A-Za-z_]\w*)',p.read_text(encoding='utf-8'),re.M):
            q=p.parent/(module+'.py')
            if q.exists():
                key=q.relative_to(ROOT).as_posix();h=sha(q)
                if key in sources:require(sources[key]==h,'dependency changed '+key)
                else:sources[key]=h
                pending.append(q)
    auditpath=ROOT/'docs/verification/uvlm-induction-adoption-audit.json';links=0
    for p in [spec,report]:
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            q=(p.parent/target.split('#')[0]).resolve();require(q.exists() or q==auditpath,'missing link '+target);links+=1
    for name,h in sources.items():require(sha(ROOT/name)==h,'changed final source '+name)
    auditpath.write_text(json.dumps(dict(passed=True,uvlm_2_1_complete=True,uvlm_review_complete=False,structural_review_complete=True,product_implementation_start=False,metrics=metrics,sources=sources,artifacts=artifacts,snapshot_context=context,legacy_recorded_status=legacy,checked_links=links,seconds=time.monotonic()-start),indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,metrics=metrics,sources=len(sources),artifacts=len(artifacts)),indent=2))
if __name__=='__main__':main()
