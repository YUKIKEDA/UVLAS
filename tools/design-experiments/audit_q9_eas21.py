"""Changed modal basis, independent gradients, boundary work and tensor energy."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import qr,solve_triangular
import audit_q9_curved_strain_screen as independent
import q9_eas21 as primary
import q9_developable_bending_approximation as mesh_source
import check_q9_curved_load as geometry
ROOT=Path(__file__).resolve().parents[2]
MU=1/2.6;LAM=.3/(1.3*.4)
MB=independent.E.copy();MB[3:]*=np.sqrt(2.)
C=2*MU*np.eye(6);C[:3,:3]+=LAM;CH=np.linalg.cholesky(C).T
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def source(q,h,points):
    _,J,g=independent.data(q,h,points);B=np.zeros((len(points),6,81))
    for a in range(3):
        for b in range(3):B[:,:,a::3]+=MB[:,a,b][None,:,None]*g[:,None,:,b]
    s,t=points[:,:2].T;gs=s**3-s;gt=t**3-t;d=np.zeros((len(s),7,3))
    for j in range(3):
        d[:,j,0]=(3*s*s-1)*t**j;d[:,j,1]=gs*j*(t**(j-1) if j else np.ones(len(t)))
        d[:,j+3,0]=j*(s**(j-1) if j else np.ones(len(s)))*gt;d[:,j+3,1]=s**j*(3*t*t-1)
    d[:,6,0]=(3*s*s-1)*gt;d[:,6,1]=gs*(3*t*t-1);grad=d@np.linalg.inv(J);M=np.zeros((len(points),6,21))
    for a in range(3):
        for b in range(3):M[:,:,a::3]+=MB[:,a,b][None,:,None]*grad[:,None,:,b]
    return J,B,M
def prepare(q,h):
    points,w=primary.quadrature(9);J,B,M=source(q,h,points);w=w*np.linalg.det(J);mean=np.einsum('p,pki->ki',w,M)/w.sum();M-=mean
    A=(np.einsum('ij,pjk->pik',CH,M)*np.sqrt(w)[:,None,None]).reshape(-1,21);W=(np.einsum('ij,pjk->pik',CH,B)*np.sqrt(w)[:,None,None]).reshape(-1,81);scale=1/np.linalg.norm(A,axis=0);Q,R=qr(A*scale,mode='economic');co=scale[:,None]*solve_triangular(R,Q.T@W)
    return q,h,mean,co
def evaluate(st,points):
    q,h,mean,co=st;J,B,M=source(q,h,points);return J,B-(M-mean)@co
def target(x,z):
    J=np.zeros((len(x),3,3));dJ=np.zeros_like(J);k=.4
    J[:,0,0]=1-k*z;J[:,1,1]=1;J[:,2,0]=k*x;J[:,0,2]=-k*x;J[:,2,2]=1
    dJ[:,0,0]=-k*x*x-z;dJ[:,2,0]=x-2*k*x*z;dJ[:,0,2]=-x;dJ[:,2,2]=-k*x*x
    G=dJ@np.linalg.inv(J);return (G+G.transpose(0,2,1))/2
def energy(E):return LAM*np.trace(E,axis1=1,axis2=2)**2+2*MU*np.einsum('pij,pij->p',E,E)
def main():
    dest=ROOT/'docs/verification/q9-eas21-audit-result.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();parents=[ROOT/'docs/verification'/f'q9-eas21-{s}-result.json' for s in ['feasibility','shared-approximation']]
    for parent in parents:
        data=json.loads(parent.read_text())
        for path,sha in data['sources'].items():assert digest(ROOT/path)==sha
    shared=json.loads(parents[1].read_text());cases=[(name,geometry.geometry(name),.01) for name in ['RECT','SADDLE','NORMAL-C']];xy,cells,coords,u,free=mesh_source.mesh(.4,2,.2)
    cases.extend((f'parabolic-{i}',coords[ids].reshape(81),.0001) for i,(_,_,ids,_) in enumerate(cells));points,w=primary.quadrature(11);local=[]
    for name,q,h in cases:
        st=prepare(q,h);J,R=evaluate(st,points);_,_,Rp,_=primary.evaluate(primary.prepare(q,h),points);difference=float(np.linalg.norm(R-Rp)/np.linalg.norm(Rp));work=np.einsum('p,pki->ik',w*np.linalg.det(J),R);boundary=independent.boundary(q,h)@np.diag([1,1,1,1/np.sqrt(2.),1/np.sqrt(2.),1/np.sqrt(2.)]);traction=float(np.linalg.norm(work-boundary)/np.linalg.norm(boundary));assert difference<=1e-9 and traction<=1e-9;local.append(dict(name=name,operator_difference=difference,boundary_difference=traction))
    observations=[]
    for record in shared['records']:
        artifact=ROOT/record['artifact'];assert digest(artifact)==record['artifact_sha256'];u=np.load(artifact)['u'];n,h=record['n'],record['h'];xy,cells,coords,_,_=mesh_source.mesh(.4,n,.2);total=error=0.
        for _,_,ids,_ in cells:
            q=coords[ids].reshape(81);J,R=evaluate(prepare(q,h),points);phi,_,_=independent.data(q,h,points);x=(phi@q.reshape(27,3))[:,0]
            # Reference physical x includes -k*x*z; recover the midsurface parameter.
            z=h*points[:,2]/2;x=x/(1-.4*z);star=target(x,z);field=np.einsum('kij,pk->pij',MB,R@u[ids].reshape(81));weight=w*np.linalg.det(J);total+=weight@energy(star);error+=weight@energy(field-star)
        err=float(np.sqrt(error/total));difference=abs(err-record['observation']['error']);assert difference<=1e-8;observations.append(dict(n=n,h=h,error=err,difference=difference));print(n,h,err,difference,flush=True)
    sources=[Path(__file__),Path(independent.__file__),ROOT/'docs/verification/q9-eas21-audit-experiment.md']+parents
    dest.write_text(json.dumps(dict(all_checks_pass=True,independent_resolve=False,local=local,observations=observations,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
if __name__=='__main__':main()
