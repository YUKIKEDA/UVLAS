"""Independent gradient/strain observation and separate condensed SVD solves."""
from pathlib import Path
import hashlib,json,time,itertools
import numpy as np
from scipy.linalg import lstsq,solve_triangular
import q3_curved_independent as independent
import developable_general_reference as reference
ROOT=Path(__file__).resolve().parents[2]
LAM=.3/(1.3*.4);MU=1/2.6
PERM=(9*np.array([5,6,9,10]+[a for a in range(16) if a not in [5,6,9,10]])[:,None]+np.arange(9)).reshape(144)
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def energy(E):return LAM*np.trace(E,axis1=1,axis2=2)**2+2*MU*np.einsum('pij,pij->p',E,E)
def observe(saved,h,u):
    p,w=np.polynomial.legendre.leggauss(11);points=np.array(list(itertools.product(p,p,p)));weight0=np.array([a*b*c for a,b,c in itertools.product(w,w,w)]);total=err=hnorm=area=0.;sp=np.array(list(itertools.product(p,p)));sw=np.array([a*b for a,b in itertools.product(w,w)]);N,g=independent.shape(sp)
    for ids in saved['ids']:
        q=saved['coords'][ids].reshape(144);phi,J,grad=independent.data(q,h,points);z=h*points[:,2]/2;x=(phi@q.reshape(48,3))[:,0]/(1-.4*z);_,star=reference.direct('quartic',.4,np.column_stack((x,np.zeros(len(x)),z)));G=np.einsum('na,pnb->pab',u[ids].reshape(48,3),grad);E=(G+G.transpose(0,2,1))/2;weight=weight0*np.linalg.det(J);total+=weight@energy(star);err+=weight@energy(E-star)
        J0=np.einsum('na,pnb->pab',saved['coords'][ids,0],g);metric=J0.transpose(0,2,1)@J0;aw=sw*np.sqrt(np.linalg.det(metric));v=N@u[ids,3:9];d=np.einsum('na,pnb->pab',u[ids,3:9],g);hnorm+=aw@(np.sum(v*v,axis=1)+np.einsum('pai,pij,paj->p',d,np.linalg.inv(metric),d));area+=aw.sum()
    return float(np.sqrt(err/total)),float(np.sqrt(hnorm/area))
def reconstruct(saved,ub):
    u=np.zeros_like(saved['u']);u[saved['boundary']]=ub.reshape(-1,9)
    for ids,R,scale in zip(saved['ids'],saved['local_R'],saved['local_scale']):
        v=u[ids].reshape(144)[PERM][36:]/scale[36:];inside=solve_triangular(R[:36,:36],R[:36,-1]-R[:36,36:144]@v,check_finite=False);u[ids[[5,6,9,10]]]=(scale[:36]*inside).reshape(4,9)
    return u
def main():
    dest=ROOT/'docs/verification/q3-quartic-shared-approximation-audit.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();parent=ROOT/'docs/verification/q3-quartic-shared-approximation-result.json';data=json.loads(parent.read_text());assert data['all_numerical_checks_pass'];records=[]
    for path,sha in data['sources'].items():assert digest(ROOT/path)==sha
    for record in data['records']:
        path=ROOT/record['artifact'];assert digest(path)==record['artifact_sha256'];saved=np.load(path);error,h1=observe(saved,record['h'],saved['u']);difference=abs(error-record['observation']['error']);h1diff=abs(h1/record['observation']['h1']-1);svd=None
        if record['n']==4:
            A=saved['A']*saved['scale'];b=saved['b'];y,residual,rank,s=lstsq(A,b,cond=1e-11,lapack_driver='gelsd');ub=saved['scale']*y;u=reconstruct(saved,ub);e,_=observe(saved,record['h'],u);qr_objective=float(np.linalg.norm(saved['A']@saved['boundary_solution']-b));svd_objective=float(np.linalg.norm(A@y-b));change=float(np.linalg.norm(saved['A']@(ub-saved['boundary_solution']))/np.linalg.norm(b));svd=dict(rank=int(rank),columns=A.shape[1],qr_objective=qr_objective,svd_objective=svd_objective,objective_difference=abs(qr_objective-svd_objective),energy_change=change,observed_error=e,observation_difference=abs(e-error));assert rank==A.shape[1] and max(svd['objective_difference'],change,svd['observation_difference'])<=1e-8
        checks=dict(observation=difference<=1e-8,h1=h1diff<=1e-8);assert all(checks.values());records.append(dict(n=record['n'],h=record['h'],error=error,error_difference=difference,h1=h1,h1_difference=h1diff,svd=svd,checks=checks));print(record['n'],record['h'],error,difference,'svd',svd,flush=True)
    sources=[Path(__file__),Path(independent.__file__),Path(reference.__file__),parent,ROOT/'docs/verification/q3-quartic-shared-approximation-experiment.md']
    dest.write_text(json.dumps(dict(records=records,all_checks_pass=True,independent_full_resolve=False,condensed_svd_resolves=2,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
if __name__=='__main__':main()
