"""Two actual curved solid patches: full trace, face work and traction PDE."""
from pathlib import Path
import hashlib,json,time,traceback,gc
import numpy as np
from scipy.linalg import svd,qr,solve_triangular
import mass_paired_ancf as shell
import exact_paired_finite as engine
from check_mass_paired_body import cofactor
engine.shell=shell
ROOT=Path(__file__).resolve().parents[2]

def trace_matrix(model,s,y):
    p=np.column_stack((np.full(len(y),s),y));n=model['n']
    return np.vstack([shell.polynomial.interpolation(p,n,d,model[key],hermite)[0] for d,key,hermite in
        [(model['r_degree'],'Rcoeff',True),(model['degree'],'Dcoeff',False),(1,'Ccoeff',False)]])

def constraints(models):
    left,right=models;cut=left['ns'];size=cut+right['ns'];y=np.linspace(0,1,4*max(m['n'] for m in models)*(max(m['degree'] for m in models)+2)+1)
    A=trace_matrix(left,0.,y);L=trace_matrix(left,1.,y);R=trace_matrix(right,0.,y)
    mat=np.vstack((np.column_stack((A,np.zeros((len(A),right['ns'])))),np.column_stack((L,-R))))
    norm=np.linalg.norm(mat,axis=1);mat=mat[norm>0]/norm[norm>0,None]
    U,s,V=svd(mat,full_matrices=True,check_finite=False);ranks=[int(sum(s>tol*s[0])) for tol in [1e-9,1e-10,1e-11]]
    assert len(set(ranks))==1
    Z=V[ranks[0]:].T
    hold=np.array([[s,y,z] for y in np.linspace(.001,.999,113) for z in [-1.,-.37,.11,.73,1.] for s in [0.]])
    NL=shell.kinematics(left,hold+np.array([1.,0.,0.]),left['h'])[0];NR=shell.kinematics(right,hold,right['h'])[0]
    defect=float(np.linalg.norm(NL@Z[:cut]-NR@Z[cut:])/max(1.,np.linalg.norm(NL),np.linalg.norm(NR)))
    fixed=shell.kinematics(left,hold,left['h'])[0]@Z[:cut]
    assert max(defect,float(np.max(abs(fixed))))<=1e-9
    return Z,dict(ranks=ranks,scalar_dimension=Z.shape[1],independent_interface_holdout=defect,fixed_holdout=float(np.max(abs(fixed))))

def local_root(model,Zlocal):
    h=model['h'];plane,wp=shell.geometry.quadrature(model['n'],model['degree']+3);z,wz=np.polynomial.legendre.leggauss(7);root=None
    for begin in range(0,len(plane),96):
        xy=plane[begin:begin+96];F=shell.fields(model,xy);Q,R,scales,weight,p=engine.universal(model,xy,z,wz,h,False)
        scalar=[F[1][0],F[2][0],F[0][1],F[1][1],F[2][1],F[0][2],F[1][2],F[2][2]]
        jets=np.zeros((len(xy),24,3*model['ns']))
        for a in range(8):
            for j in range(3):jets[:,3*a+j,j::3]=scalar[a]/scales[3*a+j]
        block=((R@jets)@Zlocal*np.sqrt(wp[begin:begin+96])[:,None,None]).reshape(-1,Zlocal.shape[1])
        root=qr(block if root is None else np.vstack((root,block)),mode='r',check_finite=False)[0][:Zlocal.shape[1]]
    return root

def loads(model,side,stress):
    order=model['degree']+5;n=model['n'];h=model['h'];out=np.zeros((model['ns'],3));internal=np.zeros_like(out)
    plane,wp=shell.geometry.quadrature(n,order);q,w=np.polynomial.legendre.leggauss(order)
    edge=np.concatenate([(j+(q+1)/2)/(2*n) for j in range(2*n)]);ew=np.tile(w/(4*n),2*n);z,wz=np.polynomial.legendre.leggauss(8)
    for axis in range(3):
        for sign in [-1,1]:
            if axis==2:p=np.column_stack((plane,np.full(len(plane),sign)));weight=wp
            else:
                p=np.zeros((len(edge)*len(z),3));p[:,axis]=(sign+1)/2;p[:,1-axis]=np.repeat(edge,len(z));p[:,2]=np.tile(z,len(edge));weight=(ew[:,None]*wz).ravel()
            N,G,B,det,X=shell.kinematics(model,p,h);J=shell.geometry.reference(p[:,:2],p[:,2],h,model['parameters'])[1]
            F=N.T@(weight[:,None]*(sign*cofactor(J,axis)@stress.T))
            if axis==0 and sign==(1 if side==0 else -1):internal+=F
            else:out+=F
    return out,internal

def main():
    dest=ROOT/'docs/verification/mass-paired-interface-result.json';assert not dest.exists();records=[];stops=[];start=time.monotonic()
    out=ROOT/'build/analysis'/('mass-paired-interface-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    for counts in [(1,1),(1,2)]:
        try:
            models=[shell.prepare(n,(.4,.2,.07,.1,-.06,.025,offset,.5),.05,False,'extended') for n,offset in zip(counts,[0.,.5])]
            cut=models[0]['ns'];Zs,trace=constraints(models);Z=np.kron(Zs,np.eye(3));cut3=3*cut
            # Independently prescribed Cartesian strain and physical surface traction.
            A=np.array([[.013,.007,-.004],[.003,-.009,.005],[.006,-.002,.011]])
            strain=np.einsum('kij,ij->k',engine.MB,(A+A.T)/2);stress=np.einsum('kij,k->ij',engine.MB,engine.C@strain)
            ua=np.vstack([model['reference_coordinates']@A.T for model in models]);outer=[];inner=[];internal=[]
            for side,model in enumerate(models):
                ext,inter=loads(model,side,stress);outer.append(ext);inner.append(inter)
                f=np.zeros_like(ext);plane,wp=shell.geometry.quadrature(model['n'],model['degree']+5);z,wz=np.polynomial.legendre.leggauss(7)
                for begin in range(0,len(plane),128):
                    p=engine.points(plane[begin:begin+128],z);N,G,B,det,X=shell.kinematics(model,p,.05)
                    f+=np.einsum('p,pnb,ab->na',(wp[begin:begin+128,None]*wz).ravel()*det,G,stress)
                internal.append(f)
            outer=np.vstack(outer);inner=np.vstack(inner);internal=np.vstack(internal)
            balance=float(np.linalg.norm(Z.T@(internal-outer).ravel())/max(1e-30,np.linalg.norm(Z.T@outer.ravel())))
            cancellation=float(np.linalg.norm(Z.T@inner.ravel())/max(1e-30,np.linalg.norm(Z.T@outer.ravel())))
            leftroot=local_root(models[0],Z[:cut3]);rightroot=local_root(models[1],Z[cut3:]);root=qr(np.vstack((leftroot,rightroot)),mode='r',check_finite=False)[0]
            scale=1/np.linalg.norm(root,axis=0);_,T=qr(root*scale,mode='economic',check_finite=False)
            singular=svd(T,compute_uv=False);ranks=[int(sum(singular>tol*singular[0])) for tol in [1e-9,1e-10,1e-11]];assert ranks==[len(scale)]*3
            residual=Z.T@(outer-internal).ravel()/.05**3
            dy=solve_triangular(T.T,scale*residual,lower=True);dx=scale*solve_triangular(T,dy);correction=(Z@dx).reshape(ua.shape)
            error=0.
            for model,du in zip(models,np.split(correction,[cut])):
                p=np.random.default_rng(193).uniform([.001,.001,-1.],[.999,.999,1.],(257,3));N=shell.kinematics(model,p,.05)[0]
                error=max(error,float(np.max(abs(N@du))))
            checks=dict(physical_equilibrium=balance<=1e-8,interface_work=cancellation<=1e-8,physical_solution=error<=1e-8,
                root_rcond=float(singular[-1]/singular[0])>=1e-10)
            artifact=out/f'n{counts[0]}-{counts[1]}.npz';np.savez_compressed(artifact,Z=Z,ua=ua,correction=correction,root=T,scale=scale,outer=outer,inner=inner)
            record=dict(elements=counts,trace=trace,stiffness_ranks=ranks,rcond=float(singular[-1]/singular[0]),equilibrium=balance,interface_work=cancellation,
                solution_error=error,checks=checks,artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(artifact.read_bytes()).hexdigest())
            records.append(record);print('INTERFACE',record,flush=True);del models;gc.collect()
        except Exception as exc:stops.append(dict(elements=counts,reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',stops[-1],flush=True)
    paths=[Path(__file__),Path(shell.__file__),Path(engine.__file__),ROOT/'tools/design-experiments/check_mass_paired_body.py']
    result=dict(passed=len(records)==2 and not stops and all(all(r['checks'].values()) for r in records),records=records,stops=stops,
        seconds=time.monotonic()-start,R4_complete=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('INTERFACE_RESULT',result['passed'],flush=True)

if __name__=='__main__':main()
