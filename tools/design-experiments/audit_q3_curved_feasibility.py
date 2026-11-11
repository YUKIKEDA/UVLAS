"""Saved Q3 witnesses: independent derivatives, manufactured strain and work."""
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
import q3_curved_independent as independent
ROOT=Path(__file__).resolve().parents[2]
MU=1/2.6;LAM=.3/(1.3*.4);ST=np.zeros((6,3,3))
for i in range(3):ST[i,i,i]=1
for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):ST[k,i,j]=ST[k,j,i]=1/np.sqrt(2.)
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def energy(E):return LAM*np.trace(E,axis1=1,axis2=2)**2+2*MU*np.einsum('pij,pij->p',E,E)
def analytic(k,x,z):
    J=np.zeros((len(x),3,3));dJ=np.zeros_like(J);J[:,0,0]=1-k*z;J[:,1,1]=1;J[:,2,0]=k*x;J[:,0,2]=-k*x;J[:,2,2]=1;dJ[:,0,0]=-k*x*x-z;dJ[:,2,0]=x-2*k*x*z;dJ[:,0,2]=-x;dJ[:,2,2]=-k*x*x;G=dJ@np.linalg.inv(J);return (G+G.transpose(0,2,1))/2
def trace(coords,ids,ij,n,h):
    rng=np.random.default_rng(20261102);u=rng.standard_normal((len(coords),9));lookup={tuple(v):cell for v,cell in zip(ij,ids)};p=np.array(list(itertools.product(np.linspace(-1,1,21),[-1.,-.5,0.,.5,1.])));maximum=0.
    for axis in [0,1]:
        for i,j in ij:
            if (i if axis==0 else j)==n-1:continue
            neighbor=(i+1,j) if axis==0 else (i,j+1);outputs=[]
            for sign,cell in [(1.,lookup[i,j]),(-1.,lookup[neighbor])]:
                points=np.zeros((len(p),3));points[:,axis]=sign;points[:,1-axis]=p[:,0];points[:,2]=p[:,1];phi,J,grad=independent.data(coords[cell].reshape(144),h,points);du=np.einsum('na,pnb->pab',u[cell].reshape(48,3),grad)@J;cov=(J.transpose(0,2,1)@du+du.transpose(0,2,1)@J)/2;select=[1-axis,2];outputs.append((phi@u[cell].reshape(48,3),du[:,:,select],cov[:,select,:][:,:,select]))
            maximum=max(maximum,max(float(np.max(abs(a-b))) for a,b in zip(*outputs)))
    return maximum
def main():
    target=ROOT/'docs/verification/q3-curved-feasibility-audit.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();parent=ROOT/'docs/verification/q3-curved-feasibility-result.json';data=json.loads(parent.read_text());assert data['all_checks_pass'];records=[]
    for p,sha in data['sources'].items():assert digest(ROOT/p)==sha
    p,w=np.polynomial.legendre.leggauss(11);points=np.array(list(itertools.product(p,p,p)));weights=np.array([a*b*c for a,b,c in itertools.product(w,w,w)])
    for record in data['records']:
        path=ROOT/record['artifact'];assert digest(path)==record['artifact_sha256'];saved=np.load(path);coords,u,ids,ij=[saved[key] for key in ['coords','u','ids','ij']];k,h,n=record['kappa'],record['h'],record['n'];local=[];work=np.zeros((len(coords),9,6));outer=np.zeros_like(work)
        for ci,(cell,(i,j)) in enumerate(zip(ids,ij)):
            q=coords[cell].reshape(144);phi,J,grad=independent.data(q,h,points);z=h*points[:,2]/2;pos=phi@q.reshape(48,3);x=pos[:,0]/(1-k*z);star=analytic(k,x,z);G=np.einsum('na,pnb->pab',u[cell].reshape(48,3),grad);E=(G+G.transpose(0,2,1))/2;weight=weights*np.linalg.det(J);err=float(np.sqrt(weight@energy(E-star)/(weight@energy(star))));difference=abs(err-record['local'][ci]['energy_error']);value=np.column_stack((-k*x**3/3-x*z,np.zeros(len(x)),x*x/2-k*x*x*z));disp=float(np.max(abs(phi@u[cell].reshape(48,3)-value)));vforce=np.einsum('p,pnb,kab->nak',weight,grad,ST).reshape(144,6);bforce=independent.boundary(q,h);work_difference=float(np.linalg.norm(vforce-bforce)/np.linalg.norm(bforce));checks=dict(strain=err<=1e-8,primary_difference=difference<=1e-8,displacement=disp<=1e-11,work=work_difference<=1e-9);assert all(checks.values());local.append(dict(cell=ci,energy_error=err,primary_difference=difference,displacement=disp,work_difference=work_difference,checks=checks));work[cell]+=vforce.reshape(16,9,6);faces=[(2,-1),(2,1)]
            for axis,index in [(0,i),(1,j)]:
                if index==0:faces.append((axis,-1))
                if index==n-1:faces.append((axis,1))
            outer[cell]+=independent.boundary(q,h,faces).reshape(16,9,6)
        shared_work=float(np.linalg.norm(work-outer)/np.linalg.norm(outer));trace_error=trace(coords,ids,ij,n,h);assert shared_work<=1e-9 and trace_error<=1e-10;records.append(dict(kappa=k,amplitude=record['amplitude'],n=n,h=h,local=local,shared_work=shared_work,actual_thickness_trace=trace_error));print(k,record['amplitude'],n,h,max(v['energy_error'] for v in local),flush=True)
    sources=[Path(__file__),Path(independent.__file__),parent,ROOT/'docs/verification/q3-curved-feasibility-experiment.md']
    target.write_text(json.dumps(dict(records=records,all_checks_pass=len(records)==16,physical_adoption=False,independent_mesh_generation=False,independent_resolve=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
if __name__=='__main__':main()
