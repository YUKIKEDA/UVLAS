"""Layered nonlinear physical affine patch, exact tangent, full mass and modes."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import solve,eigh,svd,cholesky,solve_triangular
import pg_ritz_mechanics as mechanics
import pg_ritz_trace_v3 as protected_trace
mechanics.trace=protected_trace
ROOT=rd_paths.ROOT

def rel(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))

def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-mechanics-v2-result.json';assert not dest.exists();start=time.monotonic();records=[];stops=[]
    try:
        m=mechanics;model=m.shell.prepare(1,(.4,.2,.05));q0=m.reference_coefficients(model);data=m.prepare(model);Zs,ranks,report=m.trace.coefficients(model,m.shell);Z=np.kron(Zs,np.eye(3));ua=q0@(m.FSTAR-np.eye(3)).T;external,fixed=m.external(model,.05);U,force,K,info=m.evaluate(data,ua,True);reference=rel(force-external-fixed,force);print('AFFINE',reference,flush=True)
        rng=np.random.default_rng(20261105);direction=rng.normal(size=ua.shape);gradient=np.einsum('na,pnb->pab',direction,data['G']);direction*=.003/np.max(abs(gradient));u=ua+direction;U,f,K,state=m.evaluate(data,u,True);eps=1e-28;Uc,fc,_,_=m.evaluate(data,u+1j*eps*direction);metrics=dict(energy_derivative=rel(np.array(Uc.imag/eps-np.sum(f*direction)),np.array(np.sum(f*direction))),force_derivative=rel(fc.imag/eps-(K@direction.ravel()).reshape(u.shape),fc.imag/eps),symmetric=rel(K-K.T,K));angle=.53;R=np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],[-np.sin(angle),0,np.cos(angle)]]);ur=(q0+u)@R.T-q0;Ur,fr,_,_=m.evaluate(data,ur);metrics.update(objective_energy=abs(Ur-U)/max(abs(U),1e-30),objective_force=rel(fr-f@R.T,fr));records.append(dict(case='nonlinear_layered_material_tangent_objectivity',metrics=metrics));assert max(metrics.values())<=1e-9
        split=m.prepare(model,split=True);Us,fs,Ks,_=m.evaluate(split,u,True);metrics=dict(energy=abs(Us-U)/max(abs(U),1e-30),force=rel(fs-f,f),tangent=rel(Ks-K,K),mass=rel(split['M']-data['M'],data['M']));records.append(dict(case='homogeneous_layer_split_same_S',metrics=metrics));print('SPLIT',metrics,flush=True);assert max(metrics.values())<=1e-9
        # Same physical affine traction problem; perturb all free coordinates.
        x=rng.normal(size=Z.shape[1]);du=(Z@x).reshape(ua.shape);du*=.002/np.max(abs(np.einsum('na,pnb->pab',du,data['G'])));current=ua+du;history=[]
        for iteration in range(10):
            U,fi,Ki,st=m.evaluate(data,current,True);g=Z.T@(fi-external).ravel();Kr=Z.T@Ki@Z;scale=1/np.sqrt(np.maximum(abs(np.diag(Kr)),1e-30));res=float(np.linalg.norm(scale*g)/max(1e-30,np.linalg.norm(scale*(Z.T@external.ravel()))));history.append(dict(iteration=iteration,residual=res,min_det=st['min_det']));print('NEWTON',history[-1],flush=True)
            if res<=1e-8:break
            dx=scale*solve(Kr*scale[:,None]*scale[None,:],-scale*g,assume_a='sym');dq=(Z@dx).reshape(current.shape);merit=float(U-np.sum(external*current));accepted=False
            for bt in range(14):
                trial=current+dq/2**bt;Ut,ft,_,stt=m.evaluate(data,trial)
                if stt['min_det']>0 and Ut-np.sum(external*trial)<=merit+1e-14:current=trial;accepted=True;break
            if not accepted:raise RuntimeError('nonlinear_line_search')
        U,fi,K,st=m.evaluate(data,current,True);physical_error=float(np.max(abs(data['N']@(current-ua))));metrics=dict(free_residual=history[-1]['residual'],physical_displacement_error=physical_error,reference_affine_work=reference);records.append(dict(case='nonlinear_physical_traction_equilibrium',metrics=metrics,history=history));assert max(metrics.values())<=1e-8
        scalarM=data['M'];fullM=np.kron(scalarM,np.eye(3));Mr=Z.T@fullM@Z;Kr=Z.T@K@Z;values,vectors=eigh(Kr,Mr,check_finite=False);back=rel(Kr@vectors-(Mr@vectors)*values[None,:],Kr@vectors);orth=rel(vectors.T@Mr@vectors-np.eye(len(values)),np.eye(len(values)));assert values.min()>0 and max(back,orth)<=1e-8
        reference_data=dict(data,e0=np.zeros_like(data['e0']));L=cholesky(fullM,lower=True);G=data['G'];B=np.einsum('kij,pnj->pkni',m.MB,G).reshape(len(G),6,3*model['ns']);CH=np.linalg.cholesky(data['C']).transpose(0,2,1);A=(CH@B*np.sqrt(data['weights'])[:,None,None]).reshape(-1,3*model['ns']);A=solve_triangular(L,A.T,lower=True,check_finite=False).T;sv=svd(A,compute_uv=False,check_finite=False);nullities=[int(sum(sv<=tol*sv[0])) for tol in [1e-9,1e-10,1e-11]];assert nullities==[6]*3
        # Translation and reference rotation velocities recover the physical 6x6 inertia.
        rigid=np.zeros((model['ns'],3,6));rigid[:(model['n']+2)**2,:,:3]=np.eye(3)[None];axes=np.eye(3)
        for j in range(3):rigid[:,:,3+j]=np.cross(axes[j],q0)
        Q=rigid.reshape(3*model['ns'],6);V=np.einsum('pn,nak->pak',data['N'],rigid);direct=np.einsum('p,pai,paj->ij',data['weights']*data['rho'],V,V);inertia=rel(Q.T@fullM@Q-direct,direct);assert inertia<=1e-10
        records.append(dict(case='consistent_mass_modes_and_eccentric_inertia',metrics=dict(eigen_backward=back,mass_orthogonality=orth,spatial_inertia=inertia),free_nullities=nullities,eigen_min=float(values.min()),eigen_max=float(values.max()),mass_min=float(np.linalg.eigvalsh(Mr).min())))
        out=ROOT/'build/analysis'/('pg-ritz-mechanics-'+time.strftime('%Y%m%d-%H%M%S')+'.npz');np.savez_compressed(out,q0=q0,equilibrium=current,Z=Z,M=Mr,K=Kr,omega2=values,modes=vectors,fullM=fullM);artifact=dict(artifact=out.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
    except Exception as exc:stops.append(repr(exc));artifact=None;print('STOP',repr(exc),flush=True)
    paths=[Path(__file__),Path(protected_trace.__file__),Path(mechanics.__file__),ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-mechanics-protocol.md',ROOT/'docs/specs/structure/material-layer-adoption-contract.md',ROOT/'rd/structure/main-shell/q9-dynamics/nonlinear-step-restart-contract.md'];result=dict(passed=len(records)==4 and not stops,records=records,stops=stops,state_artifact=artifact,seconds=time.monotonic()-start,R4_complete=False,main_shell_adopted=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('PASSED',result['passed'],flush=True)

if __name__=='__main__':main()
