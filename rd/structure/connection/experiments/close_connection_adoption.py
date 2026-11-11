"""Close connection 1.4 only after physical, coupled and precision gates."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,re,time
import numpy as np
ROOT=rd_paths.ROOT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(c,m):
    if not c:raise AssertionError(m)
def main():
    start=time.monotonic();dest=ROOT/'rd/structure/connection/records/connection-adoption-disposition.json';require(not dest.exists(),'immutable destination');sources={};artifacts={}
    parentpath=ROOT/'rd/structure/beam/records/beam-adoption-closure.json';parent=json.loads(parentpath.read_text());pa=json.loads((ROOT/'rd/structure/beam/records/beam-adoption-audit.json').read_text());require(parent['passed'] and parent['beam_1_3_complete'] and pa['passed'],'beam prerequisite')
    for group in ['sources','artifacts']:
        for name,h in pa[group].items():require(sha(ROOT/name)==h,'changed parent '+name)
        (sources if group=='sources' else artifacts).update(pa[group])
    for group in ['immutable','status_sync_source']:
        for name,h in parent[group].items():require(sha(ROOT/name)==h,'changed beam closure '+name);sources[name]=h
    sources[parentpath.relative_to(ROOT).as_posix()]=sha(parentpath)
    names=['connection-core','connection-beam-shell','connection-affine','connection-fold-thickness','connection-snapshot','connection-partition-mass'];paths=[rd_paths.RECORDS/f'{n}-result.json' for n in names];results=[json.loads(p.read_text()) for p in paths];require(all(r['passed'] and not r.get('failures') for r in results),'connection result gates');core,bs,affine,fold,snapshot,mass=results
    require(len(core['records'])==8 and len(bs['records'])==3 and len(affine['records'])==4 and len(fold['records'])==18 and len(snapshot['records'])==2,'coverage')
    require(len(mass['records'])==6 and all(r['rigid_energy_error']<1e-9 and r['volume_error']<1e-10 and r['min_child_mass_eigenvalue']>0 for r in mass['records']),'partition volume and physical inertia')
    for r in bs['records'][:2]:require(max(r[k] for k in ['kkt_residual','elimination_residual','physical_state_difference','reaction_difference'])<1e-8,'nonlinear KKT');require(max(s['gap'] for s in r['kkt_history'][-1:])<1e-9,'nonlinear trace')
    require(bs['records'][0]['series_error']<1e-8,'independent series')
    modes=bs['records'][2]['modal'];require([r['n'] for r in modes]==[4,8,16] and max(modes[-1]['axial_frequency_errors'])<.004 and min(modes[-1]['mac'])>.99,'independent axial modes')
    for r in modes:require(max(r[k] for k in ['spectrum_difference','mass_orthogonality','original_reduced_eigen_residual'])<1e-8,'modal algebra')
    require([r['kind'] for r in affine['records']]==['L','T','X','step'],'joint configurations')
    for r in affine['records']:require(max(r[k] for k in ['affine_field_error','free_residual','holdout_gap','support_work_difference','global_internal_wrench'])<1e-8 and r['external_surfaces_only'],'affine physics')
    finest=[];coarse=[]
    for h in [.1,.05]:
        rr=[r for r in fold['records'] if r['h']==h];require(len(rr)==9,'joint refinement grid');refs=next(x['records'] for x in fold['references'] if x['h']==h);require(len(refs)==3 and all(r['free_residual']<1e-8 for r in refs),'independent reference residual');require(abs(refs[-1]['tip']/refs[-2]['tip']-1)<.002,'reference convergence')
        fine=next(r for r in rr if r['n']==4 and r['layers']==4);base=next(r for r in rr if r['n']==4 and r['layers']==1);require(fine['tip_error']<.01 and fine['tip_error']<base['tip_error'],'folded joint quality');finest.append(fine['tip_error']);coarse.append(base['tip_error'])
        for layers in [1,2,4]:
            seq=[r['tip'] for r in rr if r['layers']==layers];require(all(b>=a*(1-1e-8) for a,b in zip(seq[:-1],seq[1:])),'surface refinement monotonicity')
        for r in rr:require(r['free_residual']<1e-8 and r['holdout_gap']<1e-9 and r['exact_full_space_average_closure'] and r['trace']['verified_kernel_bound']<1e-20,'full PG invariant subspace and trace')
    # Every selected result and its named sources is immutable. Preserve the
    # entire prior evidence chain; never rewrite failed alternatives.
    queue=paths[:];seen=set()
    while queue:
        path=queue.pop()
        if path in seen:continue
        seen.add(path);result=json.loads(path.read_text());sources[path.relative_to(ROOT).as_posix()]=sha(path)
        for name,h in result.get('sources',{}).items():
            p=ROOT/name;require(sha(p)==h,'changed connection source '+name);sources[name]=h
            if p.suffix=='.json':queue.append(p)
        if isinstance(result.get('artifact'),dict):
            a=result['artifact'];require(sha(ROOT/a['artifact'])==a['artifact_sha256'],'changed operator snapshot');artifacts[a['artifact']]=a['artifact_sha256']
    a=snapshot['artifact'];stored=np.load(ROOT/a['artifact']);require(stored['K_internal'].shape==(stored['J'].shape[1],)*2,'operator dimensions');require(np.linalg.norm(stored['R'].T@stored['R']-np.eye(3))<1e-12,'retained SO3');require(np.linalg.eigvalsh(stored['shell_mass']).min()>0,'retained actual shell mass');require(snapshot['records'][0]['relative_error']<1e-8 and snapshot['records'][0]['omitted_multiplier_hessian_work']>1e-8,'nonlinear shoe Hessian')
    # Retain exploratory failures, including the rejected ill-conditioned
    # independent reference. These logs do not substitute for frozen-source
    # acceptance results; the selected reference is re-generated and qualified.
    for name in ['connection-fold-exploration.log','connection-fold-refinement.log','connection-fold-aligned-full.log','connection-reference-corrected.log','connection-modal-exploration.log','connection-modal-factored.log']:
        p=ROOT/'build/analysis'/name;require(p.exists(),'diagnostic log');artifacts[p.relative_to(ROOT).as_posix()]=sha(p)
    # Record all directly imported local numerical utilities as reproducibility
    # inputs, including those reached through earlier immutable source chains.
    pending=[ROOT/p for p in sources if p.endswith('.py') and (ROOT/p).exists()];visited=set()
    while pending:
        p=pending.pop()
        if p in visited:continue
        visited.add(p)
        for module in re.findall(r'^(?:from|import)\s+([A-Za-z_]\w*)',p.read_text(encoding='utf-8'),re.M):
            q=p.parent/(module+'.py')
            if q.exists():
                key=q.relative_to(ROOT).as_posix();h=sha(q)
                if key in sources:require(sources[key]==h,'source mismatch '+key)
                else:sources[key]=h
                pending.append(q)
    metrics=dict(core_groups=8,physical_affine_joints=4,nonaffine_cases=18,nonlinear_beam_shell_cases=2,modal_refinements=3,finest_fold_errors=finest,unresolved_single_thickness_errors=coarse,max_nonlinear_original_residual=max(r['kkt_residual'] for r in bs['records'][:2]),independent_axial_series_error=bs['records'][0]['series_error'],finest_axial_modal_error=max(modes[-1]['axial_frequency_errors']),max_modal_formulation_difference=max(r['spectrum_difference'] for r in modes),max_nonaffine_holdout_gap=max(r['holdout_gap'] for r in fold['records']))
    spec=ROOT/'docs/specs/structure/connection-adopted-specification.md';decisions=[f'CN-D{i:02d}' for i in range(1,10)];require(all(x in spec.read_text(encoding='utf-8') for x in decisions),'specification decisions')
    decision=dict(connection_1_4_complete=True,main_shell_adopted=True,material_1_2_complete=True,beam_1_3_complete=True,structural_review_complete=False,product_implementation_start=False,decisions=decisions,metrics=metrics,next_work='1.5 loads damping recovery and structural integration',specification=spec.relative_to(ROOT).as_posix(),scope='fixed-reference full-position bonded interfaces; explicit rigid shoes; prescribed supports; particles and eccentric rigid bodies',arbitrary_joint_stress_validation=False)
    dest.write_text(json.dumps(decision,indent=2),encoding='utf-8')
    report=ROOT/'rd/structure/connection/records/connection-adoption-results.md';report.write_text(f'''# 接続1.4の採用と終了判断

2026-11-05。**全位置接着・明示的な剛体取付・指定支持・集中／偏心慣性の数値仕様を採用し、1.4を完了した。** [仕様v1.0](../connection-adopted-specification.md)、[判定JSON](connection-adoption-disposition.json)、[監査](connection-adoption-audit.json)。次は1.5荷重・減衰・回復・構造統合。構造全体・製品実装開始は未完了。

| 検証 | 成立した内容 |
|---|---|
| [基本8群](connection-core-result.json) | 取付接線フレーム、偏心剛体のM／慣性、重力・接線、粒子、指定運動の混合時間微分、従属反力、全SE(3)閉路、質量所有と退化拒否 |
| [梁―シェル実平衡](connection-beam-shell-result.json) | 軸・混合荷重をKKTと独立な非線形消去で比較。最大原平衡残差 {metrics['max_nonlinear_original_residual']:.6g}、軸の独立直列解との差 {metrics['independent_axial_series_error']:.6g} |
| 組合せモード | 同じ拘束M/K、n4/8/16。独立な連続直列棒の対象3モードへ収束。最細誤差最大 {metrics['finest_axial_modal_error']:.6g}、0.4%以内。KKT零空間／独立消去の全スペクトル差最大 {metrics['max_modal_formulation_difference']:.6g} |
| [実接合4配置](connection-affine-result.json) | 折れ・T・交差・厚さ段差の有限厚さ物理面。独立な外面応力からの荷重、全位置、支持仕事、原平衡と合力／偶力が成立 |
| [非アフィン18ケース](connection-fold-thickness-result.json) | 2厚さ×面内n1/2/4×厚さ1/2/4分割。独立P2三角形連続体の別細分化と比較。最細誤差はh=.1で {finest[0]:.4%}、h=.05で {finest[1]:.4%}。基準1%を維持 |
| [保存演算子・二階微分](connection-snapshot-result.json) | 実平衡の状態・元のK/J/H・乗数・荷重・質量を保存。非線形取付の二階変分と乗数H項を負対照で確認 |
| [分割の質量・慣性6配置](connection-partition-mass-result.json) | 厚さ分割後の体積・質量・偏心剛体速度のエネルギーを、元の2物体の独立体積モーメントへ照合。体積の重複・欠落なし |

## 不合格を変えた内容

全位置連続と代数残差だけでは、接合の変形精度を保証できなかった。接合の足元をセル途中で拘束した最初のモデルは、h=.1・n4でも約3.43%の変位誤差。足元をパッチ境界に合わせると改善したが、厚さ1分割の最細でもh=.1で {coarse[0]:.4%}、h=.05で {coarse[1]:.4%}。少なくとも厚い側は1%基準に未達であり、未解像のまま採用しない。

採用版は同じPG基底・同じ3D弾性エネルギーと全位置接着を保ち、厚さ方向にもPG子パッチを分けた。面内・厚さの二つの収束を分けて比較し、4分割の最細が両厚さで1%に入った。材料層の積分点を増やしただけではない。任意の実務接合で4分割が常に十分だとはしない。

非アフィン比較は微小変形に線形化した3D弾性問題であり、規格化した単位荷重の変位の大きさを有限変形解とは扱わない。一様幅・ν=0では、y平均が元の3D PG空間へ厳密に含まれることを有理数係数で検査した。yに依存しない変位・uy=0の問題は、同じ接着・支持・荷重に対する3D最小化と一致する。単に平均トレースへ拘束を減らしたものではない。粗いケースは全3D計算とも照合した。一般曲面・非一様幅・ν≠0へこの削減を流用しない。

## 数値参照と失敗の保存

薄い側の最初の独立三角形格子では、浮動小数の重複座標が細いスリバーを作り、原平衡残差が1を超えた。その値は参照として拒否し、座標重複を除去した上で、全参照の原残差1e−8と細分化差0.2%以内を必須にして再生成した。元ログを監査へ保持する。

組合せ固有値の通常の二次形式eighは、同じモデルの二表現の差が最細で約1.1e−8となり基準未達。採用計算では、同じ物理ひずみ作用素の平方根を質量で白色化しSVDを用いた。人工剛性や固有値の切捨ては加えず、元K/Mの固有残差も検査した。旧結果と比較ログを保持し、基準は変えていない。

## 適用境界と引継ぎ

剛体取付は物理的な金具／断面の宣言範囲に限定する。今回の円形梁を用いた実平衡を、任意の開断面の拘束反りや局所溶接応力へ転用しない。接触・破損・未解像の接合応力は未採用。一般曲面の全トレースは既存R4/R5の根拠を継承するが、今回の折れ接合を一般曲面接合の実務受入完了とはしない。

1.5で実際の荷重・減衰・回復と構造全体の契約を統合する。3／4には動的拘束の初期化・高階微分・精度型・疎求解・状態寿命、5には用途ごとの接合収束と独立受入を渡す。未達の課題を1.4の採用済みとして隠さない。
''',encoding='utf-8')
    for p in [spec,dest,report,Path(__file__)]:sources[p.relative_to(ROOT).as_posix()]=sha(p)
    auditpath=ROOT/'rd/structure/connection/records/connection-adoption-audit.json';links=0
    for p in [spec,report]:
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            q=(p.parent/target.split('#')[0]).resolve();require(q.exists() or q==auditpath,'missing link '+target);links+=1
    for name,h in sources.items():require(sha(ROOT/name)==h,'changed final source '+name)
    auditpath.write_text(json.dumps(dict(passed=True,connection_1_4_complete=True,structural_review_complete=False,product_implementation_start=False,metrics=metrics,sources=sources,artifacts=artifacts,checked_new_links=links,seconds=time.monotonic()-start),indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,metrics=metrics,sources=len(sources),artifacts=len(artifacts)),indent=2))
if __name__=='__main__':main()
