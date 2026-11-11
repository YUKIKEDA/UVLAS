"""R4.1 closure: independently recompute gates, state equivalence and provenance."""
from pathlib import Path
import ast,hashlib,json,re
import numpy as np
import pg_physical_pde as pde
import pg_physical_pde_plane as equivalent_route
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    dest=ROOT/'docs/verification/paired-geometry-ancf-physical-pde-audit.json';assert not dest.exists();path=ROOT/'docs/verification/paired-geometry-ancf-physical-pde-batch-v3-result.json';r=json.loads(path.read_text());errors=[]
    if not r['passed'] or len(r['records'])!=16:errors.append('incomplete_physical_cases')
    expected={(mode,amp,n,h) for mode in ['rational','standard'] for amp in [0.,.2] for n in [2,4] for h in [.01,.0001]};lookup={(a['mode'],a['amplitude'],a['n'],a['h']):a for a in r['records']}
    if set(lookup)!=expected:errors.append('case_coverage')
    for a in r['records']:
        if not all(a['checks'].values()) or a['rcond']<1e-10 or a['linear']>1e-8 or a['configuration_elapsed']>300 or a['load_difference']>1e-8:errors.append('scalar_gate')
        for key in ['stationarity','projection_field_energy_difference','projection_field_position_difference','relative_work_defect','relative_force_moment_balance','relative_reaction_resultant_difference']:
            if max(a['metrics'][key])>1e-8:errors.append('metric:'+key)
        if a['metrics']['fixed_boundary']>1e-9:errors.append('fixed_boundary')
        if max(a.get('commutation_error',0.),a.get('load_value_error',0.),a.get('basis_cache_error',0.),a.get('plane_kernel_error',0.))>1e-12:errors.append('equivalent_operator')
        with np.load(ROOT/a['artifact']) as d:
            if not all(np.isfinite(d[k]).all() for k in d.files):errors.append('nonfinite_state')
    spatial=[]
    for mode in ['rational','standard']:
        for amp in [0.,.2]:
            for field in ['error','position']:
                for target in range(3):
                    for h in [.01,.0001]:
                        a=lookup[mode,amp,2,h]['metrics'][field][target];b=lookup[mode,amp,4,h]['metrics'][field][target];floor=max(a,b)<1e-12;spatial.append(dict(mode=mode,amplitude=amp,field=field,target=target,h=h,numerical_floor=floor,decreased=b<a))
                        if not floor and b>=a:errors.append('spatial_decrease')
                    for n in [2,4]:
                        a=lookup[mode,amp,n,.01]['metrics'][field][target];b=lookup[mode,amp,n,.0001]['metrics'][field][target]
                        if max(a,b)>=1e-12 and b>2*a:errors.append('thin_ratio')
                        if n==4 and max(a,b)>.01:errors.append('one_percent')
    # Same-configuration physical fields before and after the resource optimization.
    oldpath=ROOT/'docs/verification/paired-geometry-ancf-physical-pde-result.json';old=json.loads(oldpath.read_text());oa=old['records'][-1];na=lookup['rational',.2,4,.01]
    with np.load(ROOT/oa['artifact']) as d:oc=d['coef'].copy();ob=d['beta_global'].copy()
    with np.load(ROOT/na['artifact']) as d:nc=d['coef'].copy();nb=d['beta_global'].copy()
    prior=json.loads((ROOT/'docs/verification/paired-geometry-ancf-candidate-result.json').read_text());pr=next(a for a in prior['records'] if a['mode']=='rational' and a['amplitude']==.2 and a['n']==4 and a['h']==.01)
    with np.load(ROOT/pr['artifact']) as d:
        geo=dict(refs=list(d['refs']),coeff=d['geo_coeff'].copy(),V=d['geo_V'].copy(),s=d['geo_s'].copy(),rank=len(d['geo_s']),metric=d['metric'].tolist(),n=4,enrich=[dict(coeff=d['enrich0_coeff'].copy(),V=d['enrich0_V'].copy(),s=d['enrich0_s'].copy())])
    pe=pn=se=sn=0.;p=np.array([[x,y,z] for x in np.linspace(-.93,.93,7) for y in np.linspace(-.91,.91,5) for z in [-1.,-.27,.39,1.]])
    for cell,q in enumerate(geo['refs']):
        N,G,B,det,X=pde.f.c.kinematics(4,cell,q,.01,p,geo);E,U=pde.loads.global_fields(q,p,.01,.4,pde.m);ds=B@(nc-oc)+E@(nb-ob);du=(N@(nc-oc).reshape(len(nc)//3,9)).reshape(len(p),3,3)+U@(nb-ob);es=B@oc+E@ob;us=(N@oc.reshape(len(oc)//3,9)).reshape(len(p),3,3)+U@ob;pe+=np.sum(du**2);pn+=np.sum(us**2);se+=np.sum(ds**2);sn+=np.sum(es**2)
    equivalence=dict(position=float(np.sqrt(pe/pn)),strain=float(np.sqrt(se/sn)))
    if max(equivalence.values())>1e-8:errors.append('resource_route_field_equivalence')
    # Isolate series truncation from the combined quadrature/series comparison.
    series=[]
    parents=[]
    for name in ['paired-geometry-ancf-candidate-result.json','paired-geometry-ancf-spatial-result.json']:parents+=json.loads((ROOT/'docs/verification'/name).read_text())['records']
    for a in r['records']:
        pr=next(x for x in parents if x.get('mode','rational')==a['mode'] and x['amplitude']==a['amplitude'] and x['n']==a['n'] and x['h']==a['h'])
        with np.load(ROOT/pr['artifact']) as d:
            geo=dict(refs=list(d['refs']),coeff=d['geo_coeff'].copy(),V=d['geo_V'].copy(),s=d['geo_s'].copy(),rank=len(d['geo_s']),metric=d['metric'].tolist(),n=a['n'],enrich=[dict(coeff=d['enrich0_coeff'].copy(),V=d['enrich0_V'].copy(),s=d['enrich0_s'].copy())]);nr,nd,ncc=map(int,d['sizes']);model=dict(geo=geo,k=.4,n=a['n'],amp=a['amplitude'],h=a['h'],nr=nr,nd=nd,nc=ncc,ns=d['coef'].shape[0]//3)
        eighteen=pde.external(model,a['mode'],15,18)[0]
        with np.load(ROOT/a['artifact']) as d:fourteen=d['F'].copy()
        difference=float(np.linalg.norm(eighteen-fourteen)/np.linalg.norm(fourteen));series.append(difference)
        if difference>1e-10:errors.append('isolated_series_difference')
    disposition=ROOT/'docs/verification/paired-geometry-ancf-physical-pde-disposition.json';queue=[disposition];seen=set();sources={};artifacts={}
    while queue:
        p=queue.pop()
        if p in seen:continue
        seen.add(p)
        if not p.exists():errors.append('missing:'+str(p));continue
        if p.suffix!='.json':continue
        data=json.loads(p.read_text(encoding='utf-8-sig'))
        for name,expected in (data.get('sources',{}) if isinstance(data,dict) else {}).items():
            if not isinstance(expected,str) or not re.fullmatch('[0-9a-f]{64}',expected):continue
            q=ROOT/name
            if not q.exists() or sha(q)!=expected:errors.append('source_hash:'+name)
            if name in sources and sources[name]!=expected:errors.append('conflicting_source:'+name)
            sources[name]=expected;queue.append(q)
        def walk(v):
            if isinstance(v,dict):
                if isinstance(v.get('artifact'),str) and 'artifact_sha256' in v:
                    q=ROOT/v['artifact'];artifacts[v['artifact']]=v['artifact_sha256']
                    if not q.exists() or sha(q)!=v['artifact_sha256']:errors.append('artifact_hash:'+v['artifact'])
                for x in v.values():walk(x)
            elif isinstance(v,list):
                for x in v:walk(x)
        walk(data)
    for name in sources:
        if name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf-8-sig'),filename=name)
    links=0
    for name in ['README.md','docs/roadmap.md','docs/preimplementation-plan.md','docs/structural-model-adoption-review.md','docs/requirements-evidence-map.md','docs/design-decisions.md','docs/paired-geometry-ancf-physical-load-contract.md','docs/verification/paired-geometry-ancf-physical-pde-results.md']:
        p=ROOT/name
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',p.read_text(encoding='utf-8-sig')):
            target=target.strip('<>').split('#')[0]
            if not target or re.match('[a-zA-Z]+://',target):continue
            q=(p.parent/target).resolve();links+=1
            if q!=dest and not q.exists():errors.append('link:'+name+':'+target)
    result=dict(audit_passed=not errors,physical_PDE_bounded_acceptance=not errors,R4_complete=False,main_shell_adopted=False,checked_sources=len(sources),checked_artifacts=len(artifacts),checked_links=links,resource_route_equivalence=equivalence,max_isolated_series_difference=max(series),spatial=spatial,errors=errors,sources={disposition.relative_to(ROOT).as_posix():sha(disposition),Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))});dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print({k:v for k,v in result.items() if k not in ['sources','spatial']});assert not errors

if __name__=='__main__':main()
