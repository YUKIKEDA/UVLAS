"""Non-affine finite-thickness folded joint against independent P2 continuum."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,time,traceback
import numpy as np
from scipy.linalg import null_space,solve,block_diag
from scipy.sparse.linalg import spsolve
import connection_shell as pg
import connection_exact_trace as exact
ROOT=Path(__file__).resolve().parents[2]
Q=np.array([[0.,0.,-1.],[0.,1.,0.],[1.,0.,0.]])

def reference(h,division):
    from skfem import MeshTri,Basis,FacetBasis,ElementVector,ElementTriP2,LinearForm,asm,condense
    from skfem.models.elasticity import linear_elasticity
    step=h/division;xs=np.unique(np.round(np.r_[np.linspace(0,1,round(1/step)+1),1-h],12));zs=np.round(np.r_[np.linspace(-h/2,h/2,division+1),h/2+np.linspace(0,.5,round(.5/step)+1)[1:]],12)
    mesh=MeshTri.init_tensor(xs,zs);centers=mesh.p[:,mesh.t].mean(axis=1);outside=np.where((centers[1]>h/2+1e-12)&(centers[0]<1-h-1e-12))[0];mesh=mesh.remove_elements(outside).remove_unused_nodes()
    elem=ElementVector(ElementTriP2());basis=Basis(mesh,elem,intorder=5);K=asm(linear_elasticity(Lambda=0.,Mu=500.),basis)*.2
    top=mesh.facets_satisfying(lambda x:np.isclose(x[1],h/2+.5));fb=FacetBasis(mesh,elem,facets=top,intorder=5)
    @LinearForm
    def traction(v,w):return v[0]/h
    f=asm(traction,fb);D=basis.get_dofs(lambda x:np.isclose(x[0],0.)).all();Kr,fr,x,I=condense(K,f,D=D);x[I]=spsolve(Kr,fr);tip=float(f@x) # f is unit resultant; this is mean tip ux.
    residual=K@x-f;free=np.setdiff1d(np.arange(len(x)),D);err=float(np.linalg.norm(residual[free])/np.linalg.norm(f))
    if err>1e-8:raise ValueError('independent_continuum_reference_residual')
    return dict(h=h,division=division,dofs=len(x),tip=tip,free_residual=err)

def nodes(n,degree):
    return [F(i,2*n)+F(2*j+1,2*(degree+1)*2*n) for i in range(2*n) for j in range(degree+1)]

def restrictions(models,h):
    n=models[0]['n'];dd=models[0]['degree'];hf=F(str(h));ys=nodes(n,dd);zs=[F(-1),F(0),F(1)]
    # Common finite face: skin z=+h/2, x in [1-h,1]; rib s=0.
    cuts={F(-1),F(1)}
    for j in range(2*n+1):
        z=2*(1-hf/2-F(j,2*n))/hf
        if -1<z<1:cuts.add(z)
    cuts=sorted(cuts);zz=[a+(b-a)*F(2*j+1,2*(dd+1)) for a,b in zip(cuts[:-1],cuts[1:]) for j in range(dd+1)]
    pa=[(1-hf/2-hf*z/2,y,F(1)) for z in zz for y in ys];pb=[(F(0),y,z) for z in zz for y in ys]
    root=[(F(0),y,z) for y in ys for z in zs]
    C,record=exact.assemble(models,[[(0,pa,1),(1,pb,-1)],[(0,root,1)]])
    return C,record

def candidate(h,n):
    start=time.monotonic();models=[pg.flat_model(n,h,1.,.2),pg.flat_model(n,h,.5,.2)];origin=[np.zeros(3),np.array([1-h/2,0,h/2])];rot=[np.eye(3),Q]
    C,record=restrictions(models,h);Z0=null_space(C,rcond=1e-11);Z=np.kron(Z0,np.eye(3));print('TRACE',n,h,record,'free',Z.shape,flush=True)
    data=[pg.data(m,origin=o,Q=q,order=m['degree']+1) for m,o,q in zip(models,origin,rot)];K=block_diag(*(pg.linear_stiffness(d) for d in data));sizes=[3*m['ns'] for m in models];f=np.zeros(sum(sizes))
    p,w,area=pg.face(models[1],axis=0,side=1.);N,_,_=pg.point(models[1],p);load=N.T@(w*np.linalg.norm(area,axis=1)/(h*.2));f[sizes[0]::3]=load
    Kr=Z.T@K@Z;fr=Z.T@f;v=solve(Kr,fr,assume_a='pos');u=Z@v;tip=float(f@u);res=float(np.linalg.norm(Z.T@(K@u-f))/np.linalg.norm(f));constraint=float(np.linalg.norm(np.kron(C,np.eye(3))@u))
    # Unused physical points; no projected or averaged trace in this check.
    rng=np.random.default_rng(148);yy=rng.uniform(size=37);zz=rng.uniform(-1,1,37);pa=np.column_stack((1-h/2-h*zz/2,yy,np.ones(37)));pb=np.column_stack((np.zeros(37),yy,zz));Na,_,Xa=pg.point(models[0],pa,origin[0],rot[0]);Nb,_,Xb=pg.point(models[1],pb,origin[1],rot[1]);ua=u[:sizes[0]].reshape(-1,3);ub=u[sizes[0]:].reshape(-1,3)
    gap=float(np.max(np.abs(Na@ua-Nb@ub)));geometry=float(np.max(np.abs(Xa-Xb)))
    result=dict(h=h,n=n,tip=tip,free_residual=res,constraint=constraint,holdout_gap=gap,reference_face_gap=geometry,trace=record,seconds=time.monotonic()-start,dofs=len(u),reduced_dofs=len(v))
    print('PG_FOLD',json.dumps(result),flush=True);return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--n',type=int,nargs='*',default=[1,2,4]);parser.add_argument('--h',type=float,nargs='*',default=[.1,.05]);parser.add_argument('--reference-only',action='store_true');parser.add_argument('--candidate-only',action='store_true');parser.add_argument('--publish',action='store_true');args=parser.parse_args();records=[];failures=[];start=time.monotonic()
    try:
        for h in args.h:
            refs=[];candidates=[]
            if not args.candidate_only:
                for div in [2,4,8]:
                    r=reference(h,div);refs.append(r);print('REFERENCE',json.dumps(r),flush=True)
            if not args.reference_only:
                for n in args.n:candidates.append(candidate(h,n))
            records.append(dict(h=h,reference=refs,candidates=candidates))
        if args.publish:
            for row in records:
                refs=row['reference'];cc=row['candidates'];assert len(refs)==3 and [c['n'] for c in cc]==[1,2,4]
                for c in cc:c['tip_error']=abs(c['tip']/refs[-1]['tip']-1)
                assert abs(refs[-1]['tip']/refs[-2]['tip']-1)<.002
                assert cc[-1]['tip_error']<.01 and cc[-1]['tip_error']<cc[0]['tip_error']
                assert all(c['free_residual']<1e-8 and c['constraint']<1e-9 and c['holdout_gap']<1e-9 and c['reference_face_gap']<1e-12 for c in cc)
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    if args.publish:
        paths=[Path(__file__),Path(pg.__file__),Path(exact.__file__),ROOT/'docs/verification/connection-adoption-protocol.md']
        result=dict(passed=not failures,records=records,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'docs/verification/connection-fold-result.json'
        if dest.exists():raise RuntimeError('immutable result exists')
        dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('CONNECTION_FOLD_RESULT',result['passed'],flush=True)
if __name__=='__main__':main()
