"""Curved commuting-map diagnostic; not a minimum-error or impossibility proof."""
from pathlib import Path
import hashlib,json,itertools
import numpy as np
from scipy.linalg import solve
import polynomial_trace_ancf as f
m=f.m;ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def values(q,p,h,k):
    N,_=m.shape(p,0);x=(N@q.reshape(75,3)[:25])[:,0];z=h*p[:,2]/2
    return np.column_stack((-k*x**5/15-z*x**3/3,np.zeros(len(x)),x**4/12-z*k*x**4/3))
def energy(E):return .3/(1.3*.4)*np.trace(E,axis1=1,axis2=2)**2+2/2.6*np.einsum('pij,pij->p',E,E)
def main():
    dest=ROOT/'docs/verification/polynomial-trace-commutation-result.json';assert not dest.exists();records=[]
    for kind,k,amp,h in itertools.product(['weak','potential'],[0.,.4],[0.,.2],[.01,.0001]):
        mesh=m.mesh(k,2,amp);q=mesh['ref'][mesh['cells'][0]].reshape(-1);s=f.prepare(q,h,kind);phi,_=f.potential_basis(q,h,s['p']);low=phi[:,f.LOW];mu=s['mu'];u=np.zeros((135,3))
        coef=solve(low.T@(mu[:,None]*low),low.T@(mu[:,None]*values(q,s['p'],h,k)),assume_a='pos');u[:75]=f.nodal_volume(coef,h)
        for face,(axis,sign) in zip(s['faces'],f.FACES):
            L=f.face_legendre(face['points'],axis);mf=face['measure']/face['measure'].sum();coef=solve(L.T@(mf[:,None]*L),L.T@(mf[:,None]*values(q,face['points'],h,k)),assume_a='pos')
            u[face['offset']:face['offset']+15]=f.nodal_face(coef,h)
        plane,w=m.quadrature(15,2);p=np.column_stack((plane,np.zeros(len(plane))));Gw=f.evaluate_gradient(s,p,True);grad=np.einsum('na,pjn->paj',u,Gw);E=(grad+grad.transpose(0,2,1))/2
        N,g=m.shape(p,0,True);x=(N@q.reshape(75,3)[:25])[:,0];J2=np.einsum('na,pnb->pab',q.reshape(75,3)[:25],g);weight=w*np.sqrt(np.linalg.det(J2.transpose(0,2,1)@J2))
        a=np.column_stack((np.ones(len(x)),np.zeros(len(x)),k*x));B1=-x[:,None,None]**2*a[:,:,None]*a[:,None,:]/(1+k*k*x*x)[:,None,None]
        midnorm=float(np.sqrt(weight@energy(E)));bending=float(np.sqrt(weight@energy(B1)));stabilization=sum(s['tau']*face['measure']@np.sum((face['jump']@u)**2,axis=1) for face in s['faces'])
        records.append(dict(kind=kind,kappa=k,amplitude=amp,h=h,mid_strain_defect=midnorm,normalized_mid_defect=float(np.sqrt(12)*midnorm/(h*bending)),
            scaled_stabilization=float(stabilization/(h**3*bending*bending/12)),exact_reference_mid_strain=0.,scope='L2 reduction of exact quartic isometry, not a solved minimum'))
    sources=[Path(__file__),Path(f.__file__),ROOT/'docs/verification/polynomial-trace-commutation-protocol.md',ROOT/'docs/verification/polynomial-trace-commutation-preparation.json']
    dest.write_text(json.dumps(dict(records=records,finite_h_lower_bound=False,unique_cause_proven=False,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources}),indent=2),encoding='utf-8')
    print([(r['kind'],r['kappa'],r['amplitude'],r['h'],r['mid_strain_defect'],r['normalized_mid_defect']) for r in records])
if __name__=='__main__':main()
