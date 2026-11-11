"""Shared constant-stress patch work against exterior surface integration."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
import q9_eas21 as method
import q9_developable_bending_approximation as reference
import audit_q9_curved_strain_screen as independent
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def boundary(q,h,faces):
    p,w=np.polynomial.legendre.leggauss(9);force=np.zeros((27,3,6));stress=independent.ST.copy();stress[3:]/=np.sqrt(2.)
    for axis,sign in faces:
        points=np.zeros((81,3));points[:,axis]=sign;other=[a for a in range(3) if a!=axis];points[:,other]=np.array(list(itertools.product(p,p)));weights=np.array([a*b for a,b in itertools.product(w,w)]);N,J,_=independent.data(q,h,points);area=sign*np.cross(J[:,:,(axis+1)%3],J[:,:,(axis+2)%3]);traction=np.einsum('kab,pb->pak',stress,area);force+=np.einsum('p,pi,pak->iak',weights,N,traction)
    return force.reshape(81,6)
def main():
    target=ROOT/'rd/structure/main-shell/q9-curved/records/q9-eas21-shared-work-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();records=[];h=.0001;n=2;points,w=method.quadrature(11)
    for k in [0.,.4]:
        xy,cells,coords,u,free=reference.mesh(k,n,.2);volume=np.zeros((len(xy),9,6));external=np.zeros_like(volume);constant=0.
        for i,j,ids,_ in cells:
            q=coords[ids].reshape(81);state=method.prepare(q,h);J,_,R,_=method.evaluate(state,points);weight=w*np.linalg.det(J);affine=np.einsum('kij,nj->nik',method.MB,q.reshape(27,3)).reshape(81,6);constant=max(constant,float(np.max(abs(R@affine-np.eye(6)))))
            faces=[(2,-1),(2,1)]
            for axis,index in [(0,i),(1,j)]:
                if index==0:faces.append((axis,-1))
                if index==n-1:faces.append((axis,1))
            volume[ids]+=np.einsum('p,pki->ik',weight,R).reshape(9,9,6);external[ids]+=boundary(q,h,faces).reshape(9,9,6)
        defect=float(np.linalg.norm(volume-external)/np.linalg.norm(external));checks=dict(constant=constant<=1e-9,shared_work=defect<=1e-9);assert all(checks.values());records.append(dict(kappa=k,n=n,h=h,constant_error=constant,exterior_work_error=defect,checks=checks))
    sources=[Path(__file__),Path(method.__file__),Path(independent.__file__),Path(reference.__file__),ROOT/'rd/structure/main-shell/q9-curved/records/q9-eas21-shared-work-experiment.md']
    target.write_text(json.dumps(dict(records=records,all_checks_pass=True,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8');print(records)
if __name__=='__main__':main()
