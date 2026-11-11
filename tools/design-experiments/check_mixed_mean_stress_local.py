"""Recheck all local gates for MMS-1 after the necessary approximation gate."""
from pathlib import Path
from types import SimpleNamespace
import hashlib,json,time
import numpy as np
import mixed_mean_stress_element as f
import check_moment_enriched_local_v2 as prior
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]
def relative(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))
def evaluate(q,name,h,k,out):
    adapter=SimpleNamespace(**vars(f));adapter.prepare=lambda q,h,**kwargs:f.prepare(q,h,kappa=k,**kwargs);prior.f=adapter;r=prior.evaluate(q,name,h,out);d=np.load(ROOT/r['artifact']);s=f.prepare(q,h,k);si=f.prepare(q,h,k,order=15,zorder=9,independent=True);T=np.eye(s['nscalar'])
    for a,b in zip(s['faces'],si['faces']):
        sl=slice(b['offset'],b['offset']+b['n']);T[sl,sl]=b['Mcoeff'].T@a['Mcoeff']
    U,force,K,_=f.nonlinear(si,T@d['current'],(T@d['v1'].reshape(-1,3)).reshape(-1),(T@d['v2'].reshape(-1,3)).reshape(-1));ind=dict(energy=relative(U-d['energy'],d['energy']),force=relative((T.T@force.reshape(-1,3)).reshape(-1)-d['force'],d['force']),tangent=relative(K-d['tangent'],d['tangent']))
    Fstar=np.array([[1.003,.002,0.],[0.,1.003,.001],[0.,0.,1.003]]);E=np.einsum('kij,ij->k',f.m.MB,(Fstar.T@Fstar-np.eye(3))/2);S=E@f.m.C.T;P=Fstar@np.einsum('kij,k->ij',f.m.MB,S);z=np.zeros(3*s['nscalar']);Up,fp,_,_=f.nonlinear(s,s['X0']@Fstar.T,z,z);expected=(f.boundary_map(s)@P.T).reshape(-1)
    rng=np.random.default_rng(1104);a=rng.normal(size=(len(s['p']),6));b=rng.normal(size=a.shape);Sa,_,_=f.project(s,a);Sb,_,_=f.project(s,b);aa=np.einsum('p,pk,pk->',s['weight'],a,Sb);bb=np.einsum('p,pk,pk->',s['weight'],b,Sa);idem,_,_=f.project(s,Sa@np.linalg.inv(f.m.C).T)
    metrics=dict(constant_F_energy=relative(Up-s['volume']*E@S/2,s['volume']*E@S/2),constant_P_force=relative(fp-expected,expected),self_adjoint=float(abs(aa-bb)/max(1.,abs(aa),abs(bb))),idempotence=relative(idem-Sa,Sa))
    r.update(kappa=k,independent=ind,mixed_metrics=metrics,stress_space=s['stress_space'],stress_gram_condition=s['stress_gram_condition']);r['checks'].update(independent=max(ind.values())<=1e-8,mixed=max(metrics.values())<=1e-9,stress_gram=s['stress_gram_condition']<=1e10);return r
def main():
    dest=ROOT/'docs/verification/mixed-mean-stress-local-result.json';assert not dest.exists();gate=json.loads((ROOT/'docs/verification/mixed-mean-stress-thin-v6-result.json').read_text());assert gate['necessary_approximation_gate']
    start=time.monotonic();out=ROOT/'build/analysis'/('mixed-mean-stress-local-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();records=[];stops=[]
    for k,a in [(0.,0.),(0.,.2),(.4,0.),(.4,.2)]:
        grid=f.m.mesh(k,2,a);q=grid['ref'][grid['cells'][0]].reshape(-1);name=f'k{k:g}-a{a:g}'
        for h in [.01,.001,.0001]:
            try:
                tick=time.monotonic();r=evaluate(q,name,h,k,out);r['checks']['resources']=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3;records.append(r);print(name,h,r['checks'],r['mixed_metrics'],r['independent'],flush=True)
                if not all(r['checks'].values()):raise RuntimeError('MMS1_local_gate')
            except Exception as exc:stops.append(dict(geometry=name,h=h,reason=repr(exc)));print('STOP',repr(exc),flush=True);break
        if stops:break
    paths=[Path(__file__),Path(f.__file__),Path(f.space.__file__),Path(prior.__file__),ROOT/'docs/verification/mixed-mean-stress-local-protocol.md',ROOT/'docs/verification/mixed-mean-stress-thin-v6-result.json']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_local_checks=len(records)==12 and not stops,main_shell_adopted=False,conditional_benchmark_only=True,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
