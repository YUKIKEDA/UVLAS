"""3D physical surface traction tests of finite-face L/T/X/step joints."""
from pathlib import Path
from fractions import Fraction as F
import json,time,hashlib,traceback
import numpy as np
from scipy.linalg import block_diag,null_space,solve
import connection_shell as pg
import connection_exact_trace as exact
from check_connection_fold import nodes,Q
ROOT=Path(__file__).resolve().parents[2]

def case(kind):
    n=1;h=F(1,10);hb=F(1,20) if kind=='step' else h;models=[pg.flat_model(n,float(h),1.,.2)];origins=[np.zeros(3)];rot=[np.eye(3)];pairs=[]
    if kind=='step':
        models.append(pg.flat_model(n,float(hb),1.,.2));origins.append(np.array([1.,0.,0.]));rot.append(np.eye(3));pairs.append((0,1,lambda y,z:(1,y,hb/h*z),lambda y,z:(0,y,z),np.array([float(hb)*.2/2,0,0]),[F(-1),F(1)]))
    else:
        center=1-h/2 if kind=='L' else F(1,2)
        for sign in ([1,-1] if kind=='X' else [1]):
            j=len(models);models.append(pg.flat_model(n,float(h),.5,.2));origins.append(np.array([float(center),0.,sign*float(h)/2]));rot.append(Q if sign==1 else np.array([[0,0,1.],[0,1.,0],[-1.,0,0]]))
            cuts={F(-1),F(1)}
            for k in range(2*n+1):
                z=2*(center-F(k,2*n))/(sign*h)
                if -1<z<1:cuts.add(z)
            pairs.append((0,j,lambda y,z,c=center,s=sign:(c-s*h*z/2,y,s),lambda y,z:(0,y,z),np.array([0.,0.,sign*float(h)*.2/2]),sorted(cuts)))
    ys=nodes(n,9);blocks=[[(0,[(F(0),y,z) for y in ys for z in [F(-1),F(0),F(1)]],1)]]
    for i,j,a,b,area,cuts in pairs:
        zz=[l+(r-l)*F(2*k+1,20) for l,r in zip(cuts[:-1],cuts[1:]) for k in range(10)];blocks.append([(i,[a(y,z) for z in zz for y in ys],1),(j,[b(y,z) for z in zz for y in ys],-1)])
    C,record=exact.assemble(models,blocks);Z=np.kron(null_space(C),np.eye(3));data=[pg.data(m,origin=o,Q=q,order=10) for m,o,q in zip(models,origins,rot)];K=block_diag(*(pg.linear_stiffness(d) for d in data));offset=np.cumsum([0]+[3*m['ns'] for m in models]);A=np.array([[.01,.002,-.003],[.004,-.006,.001],[.002,-.001,.008]]);sigma=1000*(A+A.T)/2;ua=np.concatenate([(d['q0']@A.T).ravel() for d in data]);f=np.zeros(len(ua))
    for i,(m,o,q) in enumerate(zip(models,origins,rot)):
        fi=np.zeros((m['ns'],3))
        for axis in range(3):
            for side in [0.,1.]:
                p,w,area=pg.face(m,axis,side);N,_,_=pg.point(m,p);fi+=N.T@(w[:,None]*((area@q.T)@sigma.T))
        f[offset[i]:offset[i+1]]=fi.ravel()
    # Remove both internal surface tractions. Remaining loads are external only.
    xx,ww=np.polynomial.legendre.leggauss(12);yy=np.r_[(xx+1)/4,(xx+3)/4];wy=np.tile(ww/4,2)
    for i,j,a,b,area,cuts in pairs:
        zz=np.concatenate([float(l)+(float(r-l))*(xx+1)/2 for l,r in zip(cuts[:-1],cuts[1:])]);wz=np.concatenate([float(r-l)*ww/2 for l,r in zip(cuts[:-1],cuts[1:])]);Y,Zz=np.meshgrid(yy,zz,indexing='ij');w=(wy[:,None]*wz).ravel();pa=np.array([a(y,z) for y,z in zip(Y.ravel(),Zz.ravel())],float);pb=np.array([b(y,z) for y,z in zip(Y.ravel(),Zz.ravel())],float)
        for index,p,sign in [(i,pa,1),(j,pb,-1)]:
            N,_,_=pg.point(models[index],p);f[offset[index]:offset[index+1]]-=(N.T@(w[:,None]*((sign*area)@sigma.T)[None,:])).ravel()
    rr=Z.T@(K@ua-f);du=Z@solve(Z.T@K@Z,-rr,assume_a='pos');u=ua+du;free=float(np.linalg.norm(Z.T@(K@u-f))/np.linalg.norm(f));difference=float(np.linalg.norm(du)/np.linalg.norm(ua));gap=0.;rng=np.random.default_rng(812)
    for i,j,a,b,area,cuts in pairs:
        y=rng.uniform(size=39);z=rng.uniform(-1,1,39);pa=np.array([a(t,s) for t,s in zip(y,z)],float);pb=np.array([b(t,s) for t,s in zip(y,z)],float);Na,_,Xa=pg.point(models[i],pa,origins[i],rot[i]);Nb,_,Xb=pg.point(models[j],pb,origins[j],rot[j]);assert np.max(abs(Xa-Xb))<1e-12;gap=max(gap,float(np.max(abs(Na@u[offset[i]:offset[i+1]].reshape(-1,3)-Nb@u[offset[j]:offset[j+1]].reshape(-1,3)))))
    # Recover nodal reaction using the original, uncondensed force residual.
    reaction=K@u-f;work=float(abs(u@reaction-ua@reaction)/max(abs(u@f),1e-30));force=np.zeros(3);moment=np.zeros(3)
    for i,d in enumerate(data):
        fi=(f+reaction)[offset[i]:offset[i+1]].reshape(-1,3);force+=d['translation']@fi;moment+=np.sum(np.cross(d['q0'],fi),axis=0)
    balance=float(max(np.linalg.norm(force),np.linalg.norm(moment))/max(np.linalg.norm(f),1e-30));assert difference<1e-8 and free<1e-8 and gap<1e-9 and work<1e-8 and balance<1e-8
    return dict(kind=kind,trace=record,affine_field_error=difference,free_residual=free,holdout_gap=gap,support_work_difference=work,global_internal_wrench=balance,external_surfaces_only=True,physical_interface_count=len(pairs))

def main():
    start=time.monotonic();records=[];failures=[]
    try:
        for kind in ['L','T','X','step']:
            r=case(kind);records.append(r);print('AFFINE',json.dumps(r),flush=True)
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    import sys
    if '--publish' in sys.argv:
        paths=[Path(__file__),Path(pg.__file__),Path(exact.__file__),ROOT/'tools/design-experiments/check_connection_fold.py',ROOT/'docs/verification/connection-adoption-protocol.md'];dest=ROOT/'docs/verification/connection-affine-result.json';assert not dest.exists();dest.write_text(json.dumps(dict(passed=len(records)==4 and not failures,records=records,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print('CONNECTION_AFFINE_RESULT',len(records)==4 and not failures,flush=True)
if __name__=='__main__':main()
