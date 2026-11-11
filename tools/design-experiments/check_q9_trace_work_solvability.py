"""Range obstruction for a fixed interior lift versus trace-and-work lift."""
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from scipy.linalg import lstsq
from numpy.polynomial.legendre import leggauss
import check_q9_edge_moment_shared_space_v2 as screen
import q9_edge_moment as method
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rank_metrics(D,W):
    U,s,_=np.linalg.svd(D,full_matrices=False);sc=np.linalg.svd(np.column_stack((D,W)),compute_uv=False);scale=sc[0];out={}
    for t in [1e-9,1e-10,1e-11]:
        rd=int(np.count_nonzero(s>t*scale));rc=int(np.count_nonzero(sc>t*scale));Z=W-U[:,:rd]@(U[:,:rd].T@W);rz=int(np.count_nonzero(np.linalg.svd(Z,compute_uv=False)>t*scale))
        out[str(t)]=dict(direction=rd,combined=rc,compatible_w=W.shape[1]-(rc-rd),projected_w=W.shape[1]-rz,range_defect=float(np.linalg.norm(Z)/max(np.linalg.norm(W),1e-30)))
    return out
def boundary_moment(q):
    p,w=leggauss(9);M=np.zeros((2,81));area=0.;ints=np.zeros(9)
    for i,j in itertools.product(range(9),repeat=2):
        N,dx,dy=screen.shape.shape(p[i],p[j]);J=q.reshape(9,3,3)[:,0,:2].T@np.column_stack((dx,dy));a=w[i]*w[j]*np.linalg.det(J);area+=a;ints+=a*N
    for axis,sign in itertools.product(range(2),[-1.,1.]):
        for x,weight in zip(p,w):
            xy=[x,sign] if axis==1 else [sign,x];N,dx,dy=screen.shape.shape(*xy);J=q.reshape(9,3,3)[:,0,:2].T@np.column_stack((dx,dy));normal=sign*np.linalg.det(J)*np.linalg.inv(J).T[:,axis]
            M[:,2::9]+=weight*normal[:,None]*N[None,:]
    M[0,3::9]=ints;M[1,4::9]=ints
    return M/area
def constraints(name,n,h):
    xy,cells,coords,normals,_,_=screen.make_mesh(name,n);side=2*n+1;inside=[j*side+i for j in range(1,side-1) for i in range(1,side-1)];wm={node:i for i,node in enumerate(inside)}
    D=np.zeros((50*n*n,3*len(xy)));W=np.zeros((len(D),len(inside)));ED=np.zeros((8*n*n,D.shape[1]));EW=np.zeros((len(ED),len(inside)));MD=np.zeros((2*n*n,D.shape[1]));MW=np.zeros((len(MD),len(inside)));audit=0.
    for k,(_,_,ids,_) in enumerate(cells):
        q=coords[ids].reshape(81);state=method.prepare(q,h,0.,5);p,w=leggauss(5);points=np.array([[a,b,0.] for a,b in itertools.product(p,p)]);weights=np.array([a*b for a,b in itertools.product(w,w)])
        J,E,B=method.r1(state,points);d=J[:,:,2]*2/h;A=J[:,:,:2];G=A.transpose(0,2,1)@A;L=np.linalg.cholesky(G);row=2*np.einsum('pia,kij,pj->pak',A,method.MB,d);white=np.linalg.solve(L,row)/np.linalg.norm(d,axis=1)[:,None,None];C=np.einsum('pak,pkn->pan',white,method.mandel(B)[:,:,:81]);C=(C*np.sqrt(weights*np.sqrt(np.linalg.det(G)))[:,None,None]).reshape(50,81)
        edge=[]
        for axis,sign in itertools.product(range(2),[-1,1]):
            ep=np.array([[x,sign,0.] if axis==0 else [sign,x,0.] for x in [-1/np.sqrt(3),1/np.sqrt(3)]]);Je,_,Be=method.r1(state,ep);t=Je[:,:,axis];de=Je[:,:,2]*2/h;erow=2*np.einsum('pi,kij,pj->pk',t,method.MB,de)/(np.linalg.norm(t,axis=1)*np.linalg.norm(de,axis=1))[:,None];edge.extend(np.einsum('pk,pkn->pn',erow,method.mandel(Be)[:,:,:81]))
        EC=np.array(edge);vp,vw=screen.raw.trace.base.table(1,.1) if False else (None,None)
        vp,vw=__import__('q9_edge_moment_fields').quadrature(5);Jv,Ev=method.source(q,h,0.,vp);measure=vw*np.linalg.det(Jv);MC=np.einsum('p,pkn->kn',measure,method.mandel(Ev)[:,[4,3],:81])*np.sqrt(2)/measure.sum();BM=boundary_moment(q);columns=np.array([[9*j+2,9*j+3,9*j+4,9*j+5] for j in range(9)]).reshape(-1);audit=max(audit,float(np.linalg.norm((MC-BM)[:,columns])/np.linalg.norm(BM[:,columns])))
        for j,node in enumerate(ids):
            for X,DC,WC,sl in [(C,D,W,slice(50*k,50*k+50)),(EC,ED,EW,slice(8*k,8*k+8)),(MC,MD,MW,slice(2*k,2*k+2))]:
                DC[sl,3*node:3*node+3]+=X[:,9*j+3:9*j+6]
                if node in wm:WC[sl,wm[node]]+=X[:,9*j:9*j+3]@normals[node]/n
    return dict(interior=rank_metrics(D,W),interior_work=rank_metrics(np.vstack((D,MD)),np.vstack((W,MW))),edge_work=rank_metrics(np.vstack((ED,MD)),np.vstack((EW,MW))),w_dofs=len(inside),boundary_audit=audit)
def main():
    target=F/'q9-trace-work-solvability-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();cases=[];stop=None
    parent=F/'q9-edge-moment-shared-ablation-result.json';assert json.loads(parent.read_text())['completed']
    try:
        configs=[('TRAP',n,h) for h,n in itertools.product([.1,.01,.001],[2,4])]+[(name,n,.01) for name,n in itertools.product(['RECT','SKEW'],[2,4])]
        for name,n,h in configs:
            a=constraints(name,n,h);a.update(geometry=name,n=n,h=h);cases.append(a);print(json.dumps(dict(geometry=name,n=n,h=h,w=a['w_dofs'],dimensions={key:a[key]['1e-10']['compatible_w'] for key in ['interior','interior_work','edge_work']},audit=a['boundary_audit'])),flush=True)
            assert time.monotonic()-start<300 and screen.memory.memory_bytes()<1024**3
    except Exception as exc:stop=repr(exc)
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),protocol_sha256=digest(F/'q9-trace-work-solvability-experiment.md'),contract_sha256=digest(ROOT/'docs/q9-trace-work-solvability-contract.md'),dependencies={Path(m.__file__).name:digest(m.__file__) for m in [screen,method,screen.raw,screen.shape]},cases=cases,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=screen.memory.memory_bytes(),diagnostic_checks_pass=bool(stop is None and len(cases)==10 and all(c['boundary_audit']<=1e-10 and all(r['compatible_w']==r['projected_w'] for key in ['interior','interior_work','edge_work'] for r in c[key].values()) for c in cases)),product_shell_adopted=False)
    result['edge_work_all_w']=bool(result['diagnostic_checks_pass'] and all(all(r['compatible_w']==c['w_dofs'] for r in c['edge_work'].values()) for c in cases))
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['cases','dependencies']}))
    if not result['diagnostic_checks_pass']:raise SystemExit(1)
if __name__=='__main__':main()
