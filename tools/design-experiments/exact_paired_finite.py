"""Same full 3D energy, lossless 24-jet QR compression along thickness.

No local singular value is discarded. The 42 strain rows (7 thickness nodes)
factor through 24 kinematic jets, so a full QR with 24 columns reduces the
stream while retaining its column space exactly. Observation uses the full
physical gradient on a different quadrature.
"""
import numpy as np
from scipy.linalg import qr,svd,solve_triangular
import exact_paired_ancf as shell
import pg_ritz_targets as target
import pg_ritz_coarse_v2 as coarse
from coefficient_compatible_finite import C,CH
MB=np.zeros((6,3,3));MB[0,0,0]=MB[1,1,1]=MB[2,2,2]=1
for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2)


def beta(model,h,variable):
    plane,wp=shell.geometry.quadrature(model['n'],13);z,wz=np.polynomial.legendre.leggauss(7);root=None
    for start in range(0,len(plane),256):
        xy=plane[start:start+256];p=np.column_stack((np.repeat(xy,len(z),axis=0),np.tile(z,len(xy))))
        E,U=target.fields(p,h,model['parameters'],variable);Eg,Ug=coarse.fields(model,shell,p,h,variable)
        det=shell.geometry.reference(p[:,:2],p[:,2],h,model['parameters'],variable)[2]
        weight=np.sqrt((wp[start:start+256,None]*wz).ravel()*det)[:,None,None]/h**1.5
        block=np.column_stack(((CH@Eg*weight).reshape(-1,8),(CH@E*weight).reshape(-1,3)))
        root=qr(block if root is None else np.vstack((root,block)),mode='r',check_finite=False)[0][:11]
    return solve_triangular(root[:8,:8],root[:8,8:],check_finite=False)


def points(plane,z):return np.column_stack((np.repeat(plane,len(z),axis=0),np.tile(z,len(plane))))


def kinematics_from_fields(model,plane,z,h,variable,F):
    p=points(plane,z);X,J,det,H,Ha=shell.geometry.reference(p[:,:2],p[:,2],h,model['parameters'],variable)
    back=np.repeat(np.arange(len(plane)),len(z));parts=[a[:,back] for a in F]
    zz=H*p[:,2]/2;za=Ha*p[:,2,None]/2
    N=parts[0][0]+zz[:,None]*parts[0][1]+zz[:,None]**2/2*parts[0][2]
    G=np.empty((len(p),model['ns'],3))
    for a in range(2):G[:,:,a]=parts[a+1][0]+zz[:,None]*parts[a+1][1]+zz[:,None]**2/2*parts[a+1][2]+za[:,a,None]*(parts[0][1]+zz[:,None]*parts[0][2])
    G[:,:,2]=H[:,None]/2*(parts[0][1]+zz[:,None]*parts[0][2])
    return p,N,G@np.linalg.inv(J),det,X


def universal(model,plane,z,wz,h,variable):
    _,A,normal,dn,g=shell.geometry.graph(plane,model['parameters'])
    p=points(plane,z);_,_,det,H,Ha=shell.geometry.reference(p[:,:2],p[:,2],h,model['parameters'],variable)
    zz=(H*p[:,2]/2).reshape(len(plane),len(z));J=np.empty((len(plane),len(z),3,3))
    J[:,:,:,:2]=A[:,None]+zz[:,:,None,None]*dn[:,None];J[:,:,:,2]=normal[:,None]
    inv=np.linalg.inv(J);factors=[np.ones_like(zz)]*3+[zz]*3+[zz**2/2]*2;directions=[0,1,2,0,1,2,0,1]
    T=np.empty((len(plane),len(z),6,24))
    for j,(factor,direction) in enumerate(zip(factors,directions)):
        T[:,:,:,3*j:3*j+3]=(CH@np.einsum('kij,pzj->pzki',MB,inv[:,:,direction,:]))*factor[:,:,None,None]
    weight=np.sqrt(wz[None,:]*det.reshape(len(plane),len(z)))/h**1.5
    scales=np.repeat([1.,1.,1.,1/h,1/h,1/h,1/h**2,1/h**2],3)
    T=(T*weight[:,:,None,None]).reshape(len(plane),len(z)*6,24)*scales
    Q,R=np.linalg.qr(T,mode='reduced')
    return Q,R,scales,weight,p


def solve(model,h,variable=False,order=15,beta_value=None):
    b0=beta(model,h,variable) if beta_value is None else beta_value
    plane,wp=shell.geometry.quadrature(model['n'],order);z,wz=np.polynomial.legendre.leggauss(7)
    root=None;ns=model['ns'];dofs=3*ns;massroot=None;compression=0.
    for start in range(0,len(plane),192):
        xy=plane[start:start+192];weights=wp[start:start+192];F=shell.fields(model,xy)
        Q,R,scales,weight,p=universal(model,xy,z,wz,h,variable)
        scalar_jets=[F[1][0],F[2][0],F[0][1],F[1][1],F[2][1],F[0][2],F[1][2],F[2][2]]
        jets=np.zeros((len(xy),24,dofs))
        for a in range(8):
            for j in range(3):jets[:,3*a+j,j::3]=scalar_jets[a]/scales[3*a+j]
        A=(R@jets)*np.sqrt(weights)[:,None,None]
        E,U=target.fields(p,h,model['parameters'],variable);Eg,Ug=coarse.fields(model,shell,p,h,variable)
        b=(CH@(E-Eg@b0)).reshape(len(xy),len(z),6,3)*weight[:,:,None,None]
        b=b.reshape(len(xy),len(z)*6,3);projected=Q.transpose(0,2,1)@b
        compression=max(compression,float(np.linalg.norm(b-Q@projected)/max(1.,np.linalg.norm(b))))
        block=np.column_stack((A.reshape(-1,dofs),(projected*np.sqrt(weights)[:,None,None]).reshape(-1,3)))
        root=qr(block if root is None else np.vstack((root,block)),mode='r',check_finite=False)[0][:dofs+3]
        # Four points integrate N^T N det(J) exactly in thickness (degree <=6).
        zm,wm=np.polynomial.legendre.leggauss(4)
        pm,N,G,det,X=kinematics_from_fields(model,xy,zm,h,variable,F)
        massblock=N*np.sqrt((weights[:,None]*wm).ravel()*det)[:,None]
        massroot=qr(massblock if massroot is None else np.vstack((massroot,massblock)),mode='r',check_finite=False)[0][:ns]
    scale=1/np.linalg.norm(root[:,:dofs],axis=0)
    Q,T=qr(root[:,:dofs]*scale,mode='economic',check_finite=False)
    ss=svd(T,compute_uv=False,check_finite=False)
    ranks=[int(sum(ss>tol*ss[0])) for tol in [1e-9,1e-10,1e-11]]
    if ranks!=[dofs]*3:raise ValueError('compressed_full_stiffness_rank:'+str(ranks))
    y=solve_triangular(T,Q.T@root[:,dofs:],check_finite=False);coef=scale[:,None]*y
    ms=1/np.linalg.norm(massroot,axis=0);sm=svd(massroot*ms,compute_uv=False,check_finite=False)
    mass_ranks=[int(sum(sm>tol*sm[0])) for tol in [1e-9,1e-10,1e-11]]
    if mass_ranks!=[ns]*3:raise ValueError('mass_square_root_rank:'+str(mass_ranks))
    mmin=float(sm[-1]**2)
    return dict(coef=coef,beta=b0,scale=scale,root=T,rcond=float(ss[-1]/ss[0]),mass_min=mmin,
                local_load_complement=compression,ranks=ranks,mass_root_ranks=mass_ranks,massroot=massroot)


def observe(model,h,sol,variable=False,order=17):
    plane,wp=shell.geometry.quadrature(model['n'],order);z,wz=np.polynomial.legendre.leggauss(9)
    error=np.zeros(3);reference=np.zeros(3);pe=np.zeros(3);pn=np.zeros(3);gradient=np.zeros((model['ns'],3,3))
    q=sol['coef'].reshape(model['ns'],3,3)
    for start in range(0,len(plane),256):
        xy=plane[start:start+256];F=shell.fields(model,xy)
        p,N,G,det,X=kinematics_from_fields(model,xy,z,h,variable,F)
        E,U=target.fields(p,h,model['parameters'],variable);Eg,Ug=coarse.fields(model,shell,p,h,variable)
        du=np.einsum('pnb,nat->pabt',G,q,optimize=True)
        strain=np.einsum('kab,pabt->pkt',MB,du,optimize=True)
        e=strain+Eg@sol['beta']-E;u=(N@sol['coef'].reshape(model['ns'],9)).reshape(len(p),3,3)+Ug@sol['beta']
        physical=(wp[start:start+256,None]*wz).ravel()*det;weight=physical/h**3
        stress=C@e
        error+=np.einsum('p,pkt,pkt->t',weight,e,stress)
        reference+=np.einsum('p,pkt,kl,plt->t',weight,E,C,E)
        pe+=np.einsum('p,pat,pat->t',physical,u-U,u-U);pn+=np.einsum('p,pat,pat->t',physical,U,U)
        tensor=np.einsum('kab,pkt->pabt',MB,stress)
        gradient+=np.einsum('p,pnb,pabt->nat',weight,G,tensor,optimize=True)
    residual=np.linalg.norm(sol['scale'][:,None]*gradient.reshape(-1,3),axis=0)/np.sqrt(3*model['ns']*reference)
    return dict(error=np.sqrt(error/reference).tolist(),position=np.sqrt(pe/pn).tolist(),stationarity=residual.tolist())
