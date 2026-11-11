"""Independent uncondensed QR solve for the two original rectangular cases."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import solve_triangular
import reference_work_refined_fit as fit
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;f=fit.f
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def setup(a):
    total=sum(c['ni'] for c in a['local']);free=a['free'];lookup={int(v):i for i,v in enumerate(free)};n=total+3*len(free);maps=[];offset=0
    for cell in a['local']:
        indices=list(range(offset,offset+cell['ni']));offset+=cell['ni']
        indices.extend(total+3*lookup[int(node)]+cart if int(node) in lookup else -1 for node in cell['trace_ids'] for cart in range(3));maps.append(np.array(indices))
    return maps,n,total
def solve_full(a):
    maps,n,total=setup(a);counts=[]
    for cell in a['local']:
        s=f.prepare(cell['q'],a['h']);counts.append(len(s['root']));del s
    A=np.zeros((sum(counts),n+2),order='F');offset=0
    for cell,indices,count in zip(a['local'],maps,counts):
        s=f.prepare(cell['q'],a['h']);W=s['root']/a['h']**1.5;valid=indices>=0;rows=slice(offset,offset+count);offset+=count
        A[rows,indices[valid]]=W[:,valid]
        for j,kind in enumerate(fit.KIND_TARGETS):A[rows.start:rows.start+s['elastic_rows'],n+j]=(fit.target(cell['q'],s['p'],a['h'],a['grid']['kappa'],kind)@fit.m.CH.T*np.sqrt(s['weight'])[:,None]).reshape(-1)/a['h']**1.5
        del s,W
    scale=1/np.linalg.norm(A[:,:n],axis=0);np.multiply(A[:,:n],scale[None,:],out=A[:,:n]);R=f.qrtools.rfactor(A);del A
    T=R[:n,:n];rhs=R[:n,n:n+2];y=solve_triangular(T,rhs,check_finite=False);value=scale[:,None]*y
    def unpack(x):
        local=[]
        for indices in maps:
            u=np.zeros((len(indices),2));valid=indices>=0;u[valid]=x[indices[valid]];local.append(u)
        trace=np.zeros((len(a['grid']['trace_ref']),3,2));trace[a['free']]=x[total:].reshape(-1,3,2);return dict(local=np.array(local),trace=trace)
    sol=unpack(value);metrics,grad=fit.evaluate(a,sol);history=[dict(iteration=0,**metrics)];reason=None
    for iteration in range(1,5):
        if max(metrics['direct_kkt'])<=1e-8:break
        g=np.zeros((n,2))
        for indices,gi in zip(maps,grad):valid=indices>=0;np.add.at(g,indices[valid],gi[valid])
        z=solve_triangular(T.T,scale[:,None]*g,lower=True,check_finite=False);delta=scale[:,None]*solve_triangular(T,z,check_finite=False);accepted=False
        for backtrack in range(4):
            alpha=2.**(-backtrack);trial_value=value-alpha*delta;trial=unpack(trial_value);tm,tg=fit.evaluate(a,trial)
            if max(tm['direct_kkt'])<max(metrics['direct_kkt']):value=trial_value;sol=trial;metrics,grad=tm,tg;history.append(dict(iteration=iteration,alpha=alpha,**metrics));accepted=True;break
        if not accepted:reason='original_stationarity_not_decreasing';break
    return sol,dict(history=history,stop_reason=reason,full_rcond=f.qrtools.rcond(T),**metrics)
def main():
    source=ROOT/'rd/structure/main-shell/reference-work/records/reference-work-refined-pair-result.json';parent=json.loads(source.read_text());dest=ROOT/'rd/structure/main-shell/reference-work/records/reference-work-refined-independent-result.json';assert not dest.exists();records=[];stop=None
    out=ROOT/'build/analysis'/('reference-work-refined-independent-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();start=time.monotonic()
    try:
        for r in parent['records']:
            tick=time.monotonic();a=fit.assemble(r['kappa'],r['n'],r['amplitude'],r['h']);sol,metrics=solve_full(a);obs=fit.observe(a,sol);original=np.load(ROOT/r['artifact']);field_diff=float(np.linalg.norm(sol['local']-original['local'])/np.linalg.norm(original['local']))
            error_diff=max(abs(x-y) for x,y in zip(obs['error'],r['observed']['error']));checks=dict(stationarity=max(metrics['direct_kkt'])<=1e-8,error=error_diff<=1e-8,fields=field_diff<=1e-8,resources=time.monotonic()-tick<=180 and memory.memory_bytes()<=2*1024**3)
            path=out/f"h{r['h']:g}.npz";np.savez_compressed(path,local=sol['local'],trace=sol['trace']);records.append(dict(h=r['h'],metrics=metrics,observed=obs,field_difference=field_diff,error_difference=error_diff,checks={k:bool(v) for k,v in checks.items()},seconds=time.monotonic()-tick,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path)))
            print(r['h'],metrics,'field',field_diff,'checks',checks,flush=True)
            if not all(checks.values()):raise RuntimeError('uncondensed_reference_gate_failed')
    except Exception as exc:stop=repr(exc)
    ps=[Path(__file__),Path(fit.__file__),source,ROOT/'rd/structure/main-shell/reference-work/records/reference-work-refined-protocol.md',ROOT/'rd/structure/main-shell/reference-work/records/reference-work-refined-preparation.json']
    dest.write_text(json.dumps(dict(records=records,stop_reason=stop,all_checks_pass=len(records)==2 and stop is None,independent_global_resolve=True,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in ps},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
