"""Objective upper bound from a conforming projected field; no new solve."""
from pathlib import Path
import hashlib,json
import numpy as np
import reference_work_refined_fit_v2 as fit
import diagnose_polynomial_trace_commutation as reference
ROOT=Path(__file__).resolve().parents[2];f=fit.f
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def observe(q,h,u,order,zorder,independent):
    s=f.prepare(q,h,order=order,zorder=zorder,independent=independent);star=fit.target(q,s['p'],h,.4,'quartic');e=s['B']@u.reshape(-1)-star
    elastic=s['weight']@np.einsum('pk,kl,pl->p',e,fit.m.C,e);norm=s['weight']@np.einsum('pk,kl,pl->p',star,fit.m.C,star)
    stab=sum(s['tau']*face['measure']@np.sum((face['jump']@u)**2,axis=1) for face in s['faces']);return elastic,stab,norm
def main():
    source=ROOT/'docs/verification/reference-work-refined-finite-result.json';parent=json.loads(source.read_text());dest=ROOT/'docs/verification/reference-work-refined-admissible-result.json';assert not dest.exists();records=[]
    for r in [r for r in parent['records'] if r['n']==4]:
        data=np.load(ROOT/r['artifact']);h=r['h'];values=[];canonical={};mismatch=0.;support=0.
        for cell,ids in enumerate(data['cells']):
            q=data['ref'][ids].reshape(-1);s=f.prepare(q,h);u=f.reduce(s,lambda p:reference.values(q,p,h,.4))
            for face,key in zip(s['faces'][2:],fit.edge_keys(cell,4)):
                sl=slice(face['offset'],face['offset']+face['n'])
                if key[0]==0 and key[1]==0:support=max(support,float(np.linalg.norm(u[sl])));u[sl]=0.
                if key in canonical:mismatch=max(mismatch,float(np.linalg.norm(u[sl]-canonical[key])));u[sl]=canonical[key]
                else:canonical[key]=u[sl].copy()
            values.append(u)
        primary=np.zeros(3);ind=np.zeros(3)
        for ids,u in zip(data['cells'],values):
            q=data['ref'][ids].reshape(-1);primary+=observe(q,h,u,13,7,False);ind+=observe(q,h,u,15,9,True)
        pe=float(np.sqrt((primary[0]+primary[1])/primary[2]));ie=float(np.sqrt((ind[0]+ind[1])/ind[2]));records.append(dict(n=4,h=h,initial_shared_mismatch=mismatch,initial_support_projection_norm=support,
            final_shared_mismatch=0.,final_support_norm=0.,admissible_primary_error=pe,admissible_independent_error=ie,solver_primary_error=r['fit']['error'][0],solver_independent_error=r['observed']['error'][0],
            admissible_to_solver_ratio=ie/r['observed']['error'][0],better_admissible_field=ie<r['observed']['error'][0],is_minimum=False))
    ps=[Path(__file__),Path(fit.__file__),source,ROOT/'docs/verification/reference-work-refined-admissible-protocol.md']
    dest.write_text(json.dumps(dict(records=records,new_solve=False,solution_modified=False,acceptance_thresholds_changed=False,finite_h_lower_bound=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in ps}),indent=2),encoding='utf-8');print(json.dumps(records,indent=2))
if __name__=='__main__':main()
