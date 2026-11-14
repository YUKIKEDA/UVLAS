"""Global trace tests against a same-rule strongly conforming baseline."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
import hybrid_gradient_fit as fit
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def case(k,n,amp,h,kind,baseline,out):
    tick=time.monotonic();a=fit.assemble(k,n,amp,h,kind,baseline);solution=fit.solve(a);observed=fit.observe(a,solution)
    ref=[fit.m.reference_h1(name,k,40) for name in fit.KIND_TARGETS];checks=dict(local_condition=solution['fit']['local_rcond']>=1e-10,
        global_condition=solution['fit']['global_rcond']>=1e-10,local_ranks=all(c['local_ranks']==[c['n']-6]*3 for c in a['local']),
        linear=max(solution['fit']['linear_residual'])<=1e-8,kkt=max(solution['fit']['kkt'])<=1e-8,
        quadrature=max(abs(x-y) for x,y in zip(observed['error'],solution['fit']['error']))<=1e-8,
        H1=all(x<=2*y*(1+1e-8) for x,y in zip(observed['h1'],ref)),resources=time.monotonic()-tick<=180 and memory.memory_bytes()<=2*1024**3)
    path=out/f'{"baseline" if baseline else kind}-k{k:g}-a{amp:g}-n{n}-h{h:g}.npz'
    np.savez_compressed(path,local=solution['local'],trace=solution['trace'],ref=a['original']['ref'],cells=a['original']['cells'])
    graph=fit.rigid_graph(a['grid'])
    checks['rigid_graph']=graph['passed']
    return dict(kappa=k,n=n,amplitude=amp,h=h,kind='baseline' if baseline else kind,fit=solution['fit'],observed=observed,h1_reference=ref,
        graph=graph,checks={k:bool(v) for k,v in checks.items()},seconds=time.monotonic()-tick,
        artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))

def main():
    dest=ROOT/'rd/structure/main-shell/hybrid-gradient/records/hybrid-gradient-shared-result.json';assert not dest.exists();start=time.monotonic();records=[];stop=None;eligible=[]
    parent=json.loads((ROOT/'rd/structure/main-shell/hybrid-gradient/records/hybrid-gradient-local-result.json').read_text());audit=json.loads((ROOT/'rd/structure/main-shell/hybrid-gradient/records/hybrid-gradient-local-audit.json').read_text())
    for kind in ['S4','S3+']:
        selected=[r for r in parent['records'] if r['kind']==kind];audited=[r for r in audit['records'] if r['kind']==kind]
        if len(selected)==9 and all(all(r['checks'].values()) for r in selected) and len(audited)==9 and all(r['passed'] for r in audited):eligible.append(kind)
    out=ROOT/'build/analysis'/('hybrid-gradient-shared-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    try:
        for k,amp,h in itertools.product([0.,.4],[0.,.2],[.01,.0001]):
            for kind,baseline in [('baseline',True)]+[(x,False) for x in eligible]:
                if not eligible:continue
                r=case(k,2,amp,h,kind,baseline,out);records.append(r)
                print(kind,k,amp,h,'err',r['observed']['error'],'H1',r['observed']['h1'],'fit',r['fit'],'checks',r['checks'],flush=True)
                if not r['checks']['resources'] or time.monotonic()-start>900:raise RuntimeError('global_resource_limit')
        dispositions=[]
        for kind in ['S4','S3+']:
            selected=[r for r in records if r['kind']==kind];numerical=len(selected)==8 and all(all(r['checks'].values()) for r in selected)
            pairs=[]
            for r in selected:
                b=next(x for x in records if x['kind']=='baseline' and (x['kappa'],x['amplitude'],x['h'])==(r['kappa'],r['amplitude'],r['h']))
                pairs.append(dict(kappa=r['kappa'],amplitude=r['amplitude'],h=r['h'],ratios=[x/y for x,y in zip(r['observed']['error'],b['observed']['error'])]))
            improved=numerical and all(all(x<.99 for x in p['ratios']) for p in pairs if p['kappa'])
            dispositions.append(dict(kind=kind,local_eligible=kind in eligible,all_global_checks=bool(numerical),all_curved_improved=bool(improved),comparisons=pairs,
                status='finite_thickness_refinement_required' if improved else 'rejected_necessary_gate',main_shell_adopted=False))
    except Exception as exc:stop=repr(exc);dispositions=[]
    sources=[Path(__file__),Path(fit.__file__),Path(fit.f.__file__),ROOT/'rd/structure/main-shell/hybrid-gradient/records/hybrid-gradient-ancf-protocol.md',ROOT/'rd/structure/main-shell/hybrid-gradient/records/hybrid-gradient-shared-preparation.json',ROOT/'rd/structure/main-shell/hybrid-gradient/records/hybrid-gradient-local-result.json',ROOT/'rd/structure/main-shell/hybrid-gradient/records/hybrid-gradient-local-audit.json']
    dest.write_text(json.dumps(dict(records=records,eligible_kinds=eligible,dispositions=dispositions,stop_reason=stop,physical_adoption=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
