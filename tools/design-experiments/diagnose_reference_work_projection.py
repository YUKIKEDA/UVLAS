"""RW-2 projected quartic isometry and exact reference moment observability."""
from pathlib import Path
import hashlib,json,itertools
import numpy as np
import reference_work_ancf_v2 as f
import diagnose_polynomial_trace_commutation as prior
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    dest=ROOT/'docs/verification/reference-work-projection-result.json';assert not dest.exists();records=[]
    for k,amp,h in itertools.product([0.,.4],[0.,.2],[.01,.0001]):
        grid=f.m.mesh(k,2,amp);q=grid['ref'][grid['cells'][0]].reshape(-1);s=f.prepare(q,h);u=f.reduce(s,lambda p:prior.values(q,p,h,k))
        plane,w=f.m.quadrature(15,2);p=np.column_stack((plane,np.zeros(len(plane))));G=f.evaluate_gradient(s,p,True);du=np.einsum('na,pjn->paj',u,G);E=(du+du.transpose(0,2,1))/2
        N,g=f.m.shape(p,0,True);x=(N@q.reshape(75,3)[:25])[:,0];J2=np.einsum('na,pnb->pab',q.reshape(75,3)[:25],g);weight=w*np.sqrt(np.linalg.det(J2.transpose(0,2,1)@J2))
        a=np.column_stack((np.ones(len(x)),np.zeros(len(x)),k*x));B=-x[:,None,None]**2*a[:,:,None]*a[:,None,:]/(1+k*k*x*x)[:,None,None]
        mid=float(np.sqrt(weight@prior.energy(E)));bend=float(np.sqrt(weight@prior.energy(B)));jump=sum(s['tau']*face['measure']@np.sum((face['jump']@u)**2,axis=1) for face in s['faces'])
        records.append(dict(kappa=k,amplitude=amp,h=h,mid_strain_defect=mid,normalized_mid_defect=float(np.sqrt(12)*mid/(h*bend)),scaled_stabilization=float(jump/(h**3*bend*bend/12))))
    # On the affine reference cube these three Legendre products are exactly
    # orthogonal to every volume and face moment, by tensor orthogonality.
    grid=f.m.mesh(0.,1,0.);q=grid['ref'][grid['cells'][0]].reshape(-1);s=f.prepare(q,.01)
    M=f.reduce(s,lambda p:f.basis(p,f.PID)[0]);sv=np.linalg.svd(M,compute_uv=False);blind=[f.PID.index((5,5,j)) for j in range(3)]
    obs=dict(fields=108,expected_rank=105,ranks=f.ranks(sv),blind_indices=blind,blind_products=[[5,5,j] for j in range(3)],blind_reduction_norm=float(np.linalg.norm(M[:,blind])),
        analytic_reason='Tensor Legendre orthogonality: volume is Q44P2, vertical faces are Q4P2, and top/bottom are Q44.',full_polynomial_reproduction_possible=False)
    paths=[Path(__file__),Path(f.__file__),ROOT/'docs/verification/reference-work-protocol.md',ROOT/'docs/verification/reference-work-v2-protocol.md']
    dest.write_text(json.dumps(dict(records=records,observability=obs,solved_minimum=False,finite_h_lower_bound=False,unique_cause_proven=False,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths}),indent=2),encoding='utf-8')
    print(json.dumps(dict(observability=obs,records=records),indent=2))
if __name__=='__main__':main()
