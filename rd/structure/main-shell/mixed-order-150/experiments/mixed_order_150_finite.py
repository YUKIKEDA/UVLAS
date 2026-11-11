"""Square-root condensed finite-thickness screen; no product API."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import solve_triangular
import mixed_order_150 as m
import q3_quartic_approximation as qrtools
import developable_general_reference as reference
INS=np.array([j*5+i for j in range(1,4) for i in range(1,4)]+[25+j*4+i for j in range(1,3) for i in range(1,3)]+[45])
BND=np.array([i for i in range(50) if i not in INS]);PERM=(3*np.r_[INS,BND][:,None]+np.arange(3)).reshape(150)
def h1(mesh,u,order,independent=False):
    p,w=m.quadrature(order,2);Nr,gr=m.shape(p,0,independent);total=area=0.
    for ids in mesh['cells']:
        Js=np.einsum('na,pnb->pab',mesh['ref'][ids[:25]],gr);G=Js.transpose(0,2,1)@Js;weight=w*np.sqrt(np.linalg.det(G));area+=weight.sum()
        for field in [1,2]:
            N,g=m.shape(p,field,independent);value=N@u[ids[m.SLICES[field]]];du=np.einsum('na,pnb->pab',u[ids[m.SLICES[field]]],g)
            total+=weight@(np.sum(value*value,axis=1)+np.einsum('pai,pij,paj->p',du,np.linalg.inv(G),du))
    return float(np.sqrt(total/area))
def assemble(n,h):
    mesh=m.mesh(.4,n,.2);free=np.concatenate(mesh['free']);interior=set(int(node) for ids in mesh['cells'] for node in ids[INS])
    boundary=np.array([node for node in free if node not in interior]);lookup={int(v):i for i,v in enumerate(boundary)}
    normalization=h**1.5;points,w=m.quadrature(9);Nr,_=m.shape(points,0)
    A=np.zeros((109*len(mesh['cells']),3*len(boundary)));b=np.zeros(len(A));locals=[];norm=0.;minlocal=1.;colnorm=np.zeros_like(mesh['ref'])
    for cell,ids in enumerate(mesh['cells']):
        q=mesh['ref'][ids].reshape(150);_,J,B=m.raw(q,points,h*points[:,2]/2);weight=w*np.linalg.det(J)*h/2
        x=(Nr@mesh['ref'][ids[:25]])[:,0];star=np.einsum('kij,pij->pk',m.MB,reference.strain('quartic',.4,x,h*points[:,2]/2))
        W=(np.einsum('ij,pjn->pin',m.CH,B)*np.sqrt(weight)[:,None,None]).reshape(-1,150)/normalization
        target=(star@m.CH.T*np.sqrt(weight)[:,None]).reshape(-1)/normalization;norm+=target@target
        colnorm[ids]+=np.sum(W*W,axis=0).reshape(50,3);scale=1/np.linalg.norm(W[:,PERM],axis=0)
        R=qrtools.rfactor(np.column_stack((W[:,PERM]*scale,target)));minlocal=min(minlocal,qrtools.rcond(R[:42,:42]))
        tail=R[42:,42:150]/scale[42:];sl=slice(109*cell,109*(cell+1));b[sl]=R[42:,150]
        for j,node in enumerate(ids[BND]):
            if int(node) in lookup:A[sl,3*lookup[int(node)]:3*lookup[int(node)]+3]=tail[:,3*j:3*j+3]
        locals.append(dict(ids=ids,R=R,scale=scale))
    return dict(mesh=mesh,free=free,boundary=boundary,A=A,b=b,local=locals,norm=norm,normalization=normalization,minlocal=minlocal,colnorm=colnorm)
def solve(a):
    A,b=a['A'],a['b'];scale=1/np.linalg.norm(A,axis=0);R=qrtools.rfactor(np.column_stack((A*scale,b)));nb=A.shape[1]
    T=R[:nb,:nb];y=solve_triangular(T,R[:nb,-1],check_finite=False);u=np.zeros_like(a['mesh']['ref']);u[a['boundary']]=(scale*y).reshape(-1,3)
    for cell in a['local']:
        ids,Rl,s=cell['ids'],cell['R'],cell['scale'];yb=u[ids[BND]].reshape(-1)/s[42:]
        yi=solve_triangular(Rl[:42,:42],Rl[:42,150]-Rl[:42,42:150]@yb,check_finite=False);u[ids[INS]]=(s[:42]*yi).reshape(-1,3)
    gradient=np.zeros_like(u);error=0.
    for cell in a['local']:
        ids,Rl,s=cell['ids'],cell['R'],cell['scale'];v=u[ids].reshape(-1)[PERM]/s;r=Rl[:,:150]@v-Rl[:,150];error+=r@r
        g=Rl[:,:150].T@r/s;back=np.zeros(150);back[PERM]=g;gradient[ids]+=back.reshape(50,3)
    fs=1/np.sqrt(a['colnorm'][a['free']]);kkt=float(np.linalg.norm(fs*gradient[a['free']])/(np.sqrt(3*len(a['free']))*np.sqrt(a['norm'])))
    return u,dict(error=float(np.sqrt(error/a['norm'])),global_rcond=qrtools.rcond(T),local_rcond=a['minlocal'],kkt=kkt,
        linear_residual=float(np.linalg.norm(T@y-R[:nb,-1])/np.linalg.norm(R[:nb,-1])),h1=h1(a['mesh'],u,9))
def observe(mesh,h,u,independent=True):
    p,w=m.quadrature(11);Nr,_=m.shape(p,0,independent);error=norm=0.
    for ids in mesh['cells']:
        q=mesh['ref'][ids].reshape(150);_,J,B=m.raw(q,p,h*p[:,2]/2,independent);weight=w*np.linalg.det(J)*h/2
        x=(Nr@mesh['ref'][ids[:25]])[:,0];star=reference.strain('quartic',.4,x,h*p[:,2]/2)
        observed=np.einsum('kij,pk->pij',m.MB,B@u[ids].reshape(150));delta=observed-star
        def energy(E):return .3/(1.3*.4)*np.trace(E,axis1=1,axis2=2)**2+2/2.6*np.einsum('pij,pij->p',E,E)
        error+=weight@energy(delta);norm+=weight@energy(star)
    return dict(error=float(np.sqrt(error/norm)),h1=h1(mesh,u,11,independent))
