"""Adopted PG shell in a fixed physical frame; shared verification utilities."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from functools import lru_cache
import numpy as np
from scipy.linalg import qr,svd
import mass_paired_ancf as shell
import exact_paired_finite as engine
import pg_ritz_mechanics as material
engine.shell=shell

@lru_cache(maxsize=8)
def flat_model(n=1,h=.1,span=1.,width=1.):
    # Isolated design-process instrumentation. Production preparation returns
    # this immutable factor directly; no module mutation is a product API.
    captured=[];original=shell.cholesky
    def capture(M):
        L=original(M);captured.append(L);return L
    shell.cholesky=capture
    try:model=shell.prepare(n,(0.,0.,0.,0.,0.,0.,0.,span/width),h/width,False,'extended')
    finally:shell.cholesky=original
    model=dict(model,physical_scale=width,physical_h=h,span=span,width=width)
    model['source_to_normalized']=captured[-1].inv().transpose()
    return model
def fields(model,p):
    return shell.kinematics(model,p,model['h'],False)
def data(model,E=1000.,rho=2.,origin=None,Q=None,order=None):
    origin=np.zeros(3) if origin is None else np.array(origin);Q=np.eye(3) if Q is None else np.array(Q);scale=model['physical_scale']
    if order is None:order=2*model['degree']+3
    xy,wp=shell.geometry.quadrature(model['n'],order);z,wz=np.polynomial.legendre.leggauss(4);blocks={k:[] for k in ['p','N','G','weights','X']}
    for begin in range(0,len(xy),128):
        pp,N,G,det,X=engine.kinematics_from_fields(model,xy[begin:begin+128],z,model['h'],False,shell.fields(model,xy[begin:begin+128]))
        values=dict(p=pp,N=N,G=G@Q.T/scale,weights=(wp[begin:begin+128,None]*wz).ravel()*det*scale**3,X=origin+scale*X@Q.T)
        for k in blocks:blocks[k].append(values[k])
    out={k:np.concatenate(v) for k,v in blocks.items()};N=out['N'];w=out['weights'];M=N.T@((w*rho)[:,None]*N)
    translation=np.linalg.solve(M,N.T@(w*rho));q0=scale*model['reference_coordinates']@Q.T+translation[:,None]*origin
    out.update(M=M,C=np.broadcast_to(np.eye(6)*E,(len(N),6,6)),e0=np.zeros((len(N),6)),q0=q0,translation=translation,model=model,E=E,rho=rho,origin=origin,Q=Q)
    return out
def linear_stiffness(data):
    G=data['G'];w=data['weights'];ns=G.shape[1];Gram=np.empty((3,3,ns,ns))
    for a in range(3):
        for c in range(a,3):
            Gram[a,c]=G[:,:,a].T@(w[:,None]*G[:,:,c]);Gram[c,a]=Gram[a,c].T
    trace=sum(Gram[a,a] for a in range(3));K=np.zeros((3*ns,3*ns));mu=data['E']/2
    for a in range(3):
        for c in range(3):K[a::3,c::3]=mu*((trace if a==c else 0)+Gram[c,a])
    return (K+K.T)/2
def point(model,p,origin=None,Q=None):
    N,G,B,det,X=fields(model,p);scale=model['physical_scale'];origin=np.zeros(3) if origin is None else origin;Q=np.eye(3) if Q is None else Q
    return N,G@np.array(Q).T/scale,np.array(origin)+scale*X@np.array(Q).T
def independent_rows(A,exact_rank=None):
    norms=np.linalg.norm(A,axis=1);keep=norms>1e-30;normalized=A[keep]/norms[keep,None];s=svd(normalized,compute_uv=False)
    ranks=[int(sum(s>tol*s[0])) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1 or exact_rank is not None and ranks[0]!=exact_rank:raise ValueError('unstable_trace_rank:'+str(ranks)+' exact='+str(exact_rank))
    _,_,indices=qr(normalized.T,pivoting=True,mode='economic');rows=np.flatnonzero(keep)[indices[:ranks[0]]];T=A[rows]/norms[rows,None]
    return T,rows,1/norms[rows],ranks
def face(model,axis=0,side=0.,order=None):
    if order is None:order=model['degree']+3
    y,w=np.polynomial.legendre.leggauss(order);edge=np.concatenate([(i+(y+1)/2)/(2*model['n']) for i in range(2*model['n'])]);ew=np.tile(w/(4*model['n']),2*model['n']);z,wz=np.polynomial.legendre.leggauss(4)
    if axis==2:
        xy,weights=shell.geometry.quadrature(model['n'],order);p=np.column_stack((xy,np.full(len(xy),2*side-1)))
    else:
        p=np.zeros((len(edge)*len(z),3));p[:,axis]=side;p[:,1-axis]=np.repeat(edge,len(z));p[:,2]=np.tile(z,len(edge));weights=(ew[:,None]*wz).ravel()
    X,J,det,H,Ha=shell.geometry.reference(p[:,:2],p[:,2],model['h'],model['parameters'],False);scale=model['physical_scale']
    area=(2*side-1)*np.cross(J[:,:,(axis+1)%3],J[:,:,(axis+2)%3])*scale**2
    return p,weights,area
