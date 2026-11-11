"""Frozen local/kinematic feasibility screen, without a global equilibrium solve."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import svdvals
import q9_eas21 as method
import q9_developable_bending_approximation as reference
import check_q9_curved_load as geometry
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def weighted(B,w):return (np.einsum('ij,pjk->pik',method.CH,B)*np.sqrt(w)[:,None,None]).reshape(-1,B.shape[-1])
def ranks(A):
    s=svdvals(A);return [int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]]
def check(q,h):
    st=method.prepare(q,h);points,w=method.quadrature(11);J,B,R,M=method.evaluate(st,points);w=w*np.linalg.det(J)
    A=weighted(B,w);D=weighted(R,w);scale=1/np.linalg.norm(A,axis=0)
    affine=np.einsum('kij,nj->nik',method.MB,q.reshape(27,3)).reshape(81,6)
    constant=float(np.max(abs(np.einsum('pki,ij->pkj',R,affine)-np.eye(6))))
    wb=np.einsum('p,pki->ki',w,B);wr=np.einsum('p,pki->ki',w,R);work=float(np.linalg.norm(wr-wb)/np.linalg.norm(wb))
    mean11=np.einsum('p,pki->ki',w,M)/w.sum();mode_scale=np.linalg.norm(st['mode_matrix'],axis=0)
    mean_error=float(np.linalg.norm(mean11/mode_scale)/max(1.,np.linalg.norm(st['mean']/mode_scale)))
    mode_rank=ranks(st['mode_matrix']/mode_scale)
    modes=[]
    for axis in np.eye(3):
        v=np.zeros((9,3,3));v[:,0,:]=axis;modes.extend([v.reshape(81),np.cross(np.broadcast_to(axis,(9,3,3)),q.reshape(9,3,3)).reshape(81)])
    As=D*scale
    rigid=max(float(np.linalg.norm(As@(v/scale))/(np.linalg.norm(As)*np.linalg.norm(v/scale))) for v in modes)
    p9,w9=method.quadrature(9);J9,_,R9,_=method.evaluate(st,p9);A9=weighted(R9,w9*np.linalg.det(J9))
    integration=float(np.linalg.norm(A9.T@A9-D.T@D)/np.linalg.norm(D.T@D))
    checks=dict(modes=mode_rank==[21]*3,constant=constant<=1e-9,work=work<=1e-9,mean=mean_error<=1e-8,rigid=rigid<=1e-10,integration=integration<=1e-8)
    return st,dict(raw_ranks=ranks(A*scale),eas_ranks=ranks(As),mode_ranks=mode_rank,constant=constant,work=work,mean=mean_error,rigid=rigid,integration=integration,checks=checks)
def main():
    target=ROOT/'rd/structure/main-shell/q9-curved/records/q9-eas21-feasibility-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();local=[];bending=[]
    for name in ['RECT','SKEW','TRAP','SADDLE','NORMAL-C']:
        for h in [.1,.01]:
            _,log=check(geometry.geometry(name),h);local.append(dict(geometry=name,h=h,**log));print(name,h,log['eas_ranks'],flush=True)
    for n in [1,2]:
        for amplitude in [0.,.2]:
            for h in [.001,.0001]:
                xy,cells,coords,u,free=reference.mesh(.4,n,amplitude);total=0.;raw_error=0.;eas_error=0.;mean_error=0.;logs=[]
                for _,_,ids,_ in cells:
                    q=coords[ids].reshape(81);state,log=check(q,h);logs.append(log);points,w=method.quadrature(11);J,B,R,M=method.evaluate(state,points);weight=w*np.linalg.det(J);N=np.array([method.common.shape(s,t)[0] for s,t,_ in points]);x=(N@xy[ids])[:,0];star=reference.target(.4,x,h*points[:,2]/2);v=u[ids].reshape(81);dr=np.einsum('pki,i->pk',B,v)-star;de=np.einsum('pki,i->pk',R,v)-star
                    total+=float(np.linalg.norm((star@method.CH.T)*np.sqrt(weight)[:,None])**2)
                    raw_error+=float(np.linalg.norm((dr@method.CH.T)*np.sqrt(weight)[:,None])**2)
                    eas_error+=float(np.linalg.norm((de@method.CH.T)*np.sqrt(weight)[:,None])**2)
                    mean=(weight@dr)/weight.sum();mean_error+=float(weight.sum()*np.linalg.norm(method.CH@mean)**2)
                bending.append(dict(n=n,amplitude=amplitude,h=h,raw_interpolation_error=float(np.sqrt(raw_error/total)),eas_interpolation_error=float(np.sqrt(eas_error/total)),mean_obstruction_norm=float(np.sqrt(mean_error/total)),local=logs));print('bending',n,amplitude,h,bending[-1]['eas_interpolation_error'],flush=True)
    sources=[Path(__file__),Path(method.__file__),ROOT/'rd/structure/main-shell/q9-curved/records/q9-eas21-feasibility-experiment.md',Path(reference.__file__),Path(method.raw.__file__),Path(method.common.__file__),Path(geometry.__file__)]
    out=dict(candidate='EAS21-linear-homogeneous',physical_adoption=False,local=local,bending=bending,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes())
    target.write_text(json.dumps(out,indent=2),encoding='utf-8')
if __name__=='__main__':main()
