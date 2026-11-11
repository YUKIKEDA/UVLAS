"""Value action vs frozen dense recovery and independent moment observations."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import moment_position_recovery_action as r
import check_moment_enriched_local_v2 as oldcheck
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];f=r.f
def main():
    dest=ROOT/'docs/verification/moment-position-recovery-action-result.json';assert not dest.exists();records=[];stops=[];rng=np.random.default_rng(20261104)
    for name,k,n,amp in [('flat',0.,2,0.),('curved-distorted',.4,2,.2),('quartic-reference',None,1,0.)]:
        try:
            start=time.monotonic()
            if k is None:qs=[oldcheck.prior.generic()]
            else:
                grid=f.m.mesh(k,n,amp);qs=[grid['ref'][ids].reshape(-1) for ids in grid['cells']]
            states=[r.light(f.prepare(q,.001)) for q in qs];dense=r.old.Recovery(states,n,clamp=True);global_values=rng.normal(size=(dense.size,4))
            for key,item in dense.facet.items():
                if key[0]==0 and key[1]==0:global_values[item['offset']:item['offset']+item['face']['n']]=0
            values=[N@global_values for N in dense.native];action=r.RecoveryAction(states,n,values,clamp=True);difference=moments=shared=support=0.
            for cell,s in enumerate(states):
                p,w=f.volume_points(15,9);v=action.evaluate(cell,p);reference=dense.evaluate(cell,p)@global_values;difference=max(difference,float(np.linalg.norm(v-reference)/max(1.,np.linalg.norm(reference))));V=r.modal(p,f.VID);moments=max(moments,float(np.linalg.norm(V.T@((w/8)[:,None]*v)-values[cell][:76])/max(1.,np.linalg.norm(values[cell][:76]))))
                for face in s['faces']:
                    p,w=f.old.face_points(face['axis'],face['sign'],15,9);v=action.evaluate(cell,p);target=values[cell][face['offset']:face['offset']+face['n']];moments=max(moments,float(np.linalg.norm(r.facevalue(face,p).T@((w/4)[:,None]*v)-target)/max(1.,np.linalg.norm(target))))
                    if face['axis']==0 and face['sign']==-1 and cell%n==0:support=max(support,float(np.linalg.norm(v)))
            pp=np.array([[a,b] for a in np.linspace(-1,1,9) for b in np.linspace(-1,1,9)])
            for cell in range(n*n):
                j,i=divmod(cell,n)
                for axis,other in [(0,cell+1),(1,cell+n)]:
                    if (axis==0 and i==n-1) or (axis==1 and j==n-1):continue
                    a=action.evaluate(cell,r.points(axis,1,pp));b=action.evaluate(other,r.points(axis,-1,pp));shared=max(shared,float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a))))
            probe=np.zeros((dense.size,1));col=int(np.argmax(dense.native[0][75]));probe[col]=1;extra=r.RecoveryAction(states,n,[N@probe for N in dense.native],clamp=True);p,w=f.volume_points(13,7);response=float(np.linalg.norm(extra.evaluate(0,p)));seconds=time.monotonic()-start;peak=memory.memory_bytes();passed=max(difference,moments,shared,support)<=1e-9 and response>1e-6 and seconds<=900 and peak<=2*1024**3
            records.append(dict(name=name,n=n,difference=difference,moments=moments,shared=shared,support=support,extra_response=response,global_scalar=dense.size,action_columns=action.size,global_identity_allocated=hasattr(action,'eye'),seconds=seconds,peak_bytes=peak,passed=bool(passed and not hasattr(action,'eye'))));print(records[-1],flush=True)
        except Exception as exc:stops.append(dict(name=name,reason=repr(exc)));print('STOP',repr(exc),flush=True);break
    paths=[Path(__file__),Path(r.__file__),ROOT/'tools/design-experiments/moment_position_recovery.py',ROOT/'docs/verification/moment-position-recovery-action-preparation.json']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_checks_pass=len(records)==3 and not stops and all(x['passed'] for x in records),efficient_load_adjoint_implemented=False,physical_PDE_acceptance=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
