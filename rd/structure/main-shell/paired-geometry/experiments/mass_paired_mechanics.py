"""Same verified PG space coupled to existing 3D material kernels."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from types import SimpleNamespace
import numpy as np
from scipy.linalg import svd
import mass_paired_ancf as shell
import pg_ritz_mechanics as kernel

# The existing material routine consumes only these reference operations.
kernel.shell=SimpleNamespace(lift=SimpleNamespace(geometry=shell.geometry),kinematics=shell.kinematics)
MB=kernel.MB
FSTAR=kernel.FSTAR
S0=kernel.S0
prepare=kernel.prepare
evaluate=kernel.evaluate
external=kernel.external


def fixed_trace_kernel(model):
    n=model['n'];dd=model['degree'];q,w=np.polynomial.legendre.leggauss(dd+3)
    y=np.concatenate([(j+(q+1)/2)/(2*n) for j in range(2*n)])
    points=np.column_stack((np.zeros(len(y)),y))
    R=shell.polynomial.interpolation(points,n,model['r_degree'],model['Rcoeff'],True)[0]
    D=shell.polynomial.interpolation(points,n,dd,model['Dcoeff'])[0]
    C=shell.polynomial.interpolation(points,n,1,model['Ccoeff'])[0]
    matrix=np.vstack((R,D,C));norm=np.linalg.norm(matrix,axis=1)
    matrix=matrix[norm>0]/norm[norm>0,None]
    U,s,V=svd(matrix,full_matrices=True,check_finite=False)
    ranks=[int(sum(s>tol*s[0])) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('complete_trace_rank:'+str(ranks))
    Z=V[ranks[0]:].T
    holdout=np.array([[0.,yy,z] for yy in np.linspace(0,1,113) for z in [-1.,-.29,0.,.43,1.]])
    N,G,B,det,X=shell.kinematics(model,holdout,model['h'],model['variable'])
    defect=float(np.linalg.norm(N@Z)/max(1.,np.linalg.norm(N)))
    if defect>1e-9:raise ValueError('complete_trace_holdout')
    return Z,ranks,defect


def interface_recovery(model,h,u,order=7):
    plane=np.array([[.17,.23],[.61,.47],[.87,.79]])
    p=np.column_stack((plane,np.zeros(len(plane))))
    N,G,B,det,X=shell.kinematics(model,p,h,True)
    F=np.eye(3)+np.einsum('na,pnb->pab',u,G)
    E=np.einsum('kij,pij->pk',MB,(F.transpose(0,2,1)@F-np.eye(3))/2)
    # Evaluate the material on both sides at the identical interface point.
    result=[]
    for theta,factor in [(.23,1.),(-.41,.7)]:
        _,J,_,H,Ha=shell.geometry.reference(plane,p[:,2],h,model['parameters'],True)
        normal=np.cross(J[:,:,0],J[:,:,1]);normal/=np.linalg.norm(normal,axis=1)[:,None]
        fiber=np.tile([np.cos(theta),np.sin(theta),.1],(len(p),1))
        axis=fiber-normal*np.sum(fiber*normal,axis=1)[:,None];axis/=np.linalg.norm(axis,axis=1)[:,None]
        A=np.stack((axis,np.cross(normal,axis),normal),axis=2)
        T=np.einsum('kij,pia,lab,pjb->pkl',MB,A,MB,A)
        Cg=T@(kernel.material()*factor)@T.transpose(0,2,1)
        estar=np.einsum('kij,ij->k',MB,(FSTAR.T@FSTAR-np.eye(3))/2)
        s0=np.einsum('kij,ij->k',MB,S0)
        e0=estar-np.linalg.solve(Cg,np.tile(s0,(len(p),1))[:,:,None])[:,:,0]
        stress=np.einsum('pij,pj->pi',Cg,E-e0)
        local=np.einsum('pji,pj->pi',T,E-e0)
        independent=np.einsum('pij,pj->pi',T,np.einsum('ij,pj->pi',kernel.material()*factor,local))
        defect=float(np.linalg.norm(stress-independent)/max(1.,np.linalg.norm(stress)))
        result.append(dict(stress=stress,defect=defect))
    difference=float(np.linalg.norm(result[0]['stress']-result[1]['stress']))
    if difference<=1e-8:raise ValueError('layer_recovery_negative_control_not_excited')
    return dict(max_transform_error=max(a['defect'] for a in result),two_sided_difference=difference,
                lower=result[0]['stress'].tolist(),upper=result[1]['stress'].tolist())
