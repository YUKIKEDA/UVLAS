"""Square-root interior condensation and shared Q3 quartic approximation."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from functools import lru_cache
import itertools
import numpy as np
from scipy.linalg import qr,solve_triangular
from scipy.linalg.lapack import get_lapack_funcs
import q3_curved_element as element
import developable_general_reference as reference
import check_q9_curved_best_approximation_bounds as rounding
INS=np.array([5,6,9,10]);BND=np.array([i for i in range(16) if i not in INS]);PERM=(9*np.r_[INS,BND][:,None]+np.arange(9)).reshape(144)
def rfactor(A):return qr(A,mode='r',overwrite_a=True,check_finite=False)[0][:A.shape[1]].copy()
def rcond(A):
    value,info=get_lapack_funcs('trcon',(A,))(np.asfortranarray(A))
    if info:raise ValueError('trcon_info_'+str(info))
    return float(value)
def ref_u(xy):
    x=xy[:,0];u=np.zeros((len(x),9));u[:,:3]=np.column_stack((-.4*x**5/15,np.zeros(len(x)),x**4/12));u[:,3:6]=np.column_stack((-x**3/3,np.zeros(len(x)),-.4*x**4/3));return u
@lru_cache(None)
def surface(order):
    p,w=np.polynomial.legendre.leggauss(order);pts=np.array(list(itertools.product(p,p)));N,g=element.surface(pts);weights=np.array([a*b for a,b in itertools.product(w,w)]);return N,g,weights
def reference_h1(order):
    p,w=np.polynomial.legendre.leggauss(order);x=(p+1)/2;w=w/2;D=1+.4**2*x*x;area=w@np.sqrt(D);value=w@(((x**6+.4**2*x**8)/9+(x**4+16*.4**2*x**6/9)/D)*np.sqrt(D));return float(np.sqrt(value/area)),float(area)
def h1_matrix(xy,cells,coords,free,order):
    N,g,w=surface(order);M=np.zeros((len(xy),len(xy)));area=0.
    for i,j,ids in cells:
        J=np.einsum('na,pnb->pab',coords[ids,0],g);G=J.transpose(0,2,1)@J;weight=w*np.sqrt(np.linalg.det(G));area+=weight.sum();local=np.einsum('p,pi,pj->ij',weight,N,N)+np.einsum('p,pia,pab,pjb->ij',weight,g,np.linalg.inv(G),g);M[np.ix_(ids,ids)]+=local
    return M[np.ix_(free,free)]/area,float(area)
def h1(u,free,M):return float(np.sqrt(np.sum(u[free,3:9]*(M@u[free,3:9]))))
def assemble(n,h,normalization):
    xy,cells,coords,_=element.mesh(.4,n,.2);interp=ref_u(xy);free=np.flatnonzero(np.arange(len(xy))%(3*n+1));interior=set(int(v) for _,_,ids in cells for v in ids[INS]);boundary=np.array([v for v in free if v not in interior]);lookup={int(v):i for i,v in enumerate(boundary)};A=np.zeros((109*len(cells),9*len(boundary)));b=np.zeros(109*len(cells));local=[];colnorm=np.zeros((len(xy),9));points,w=element.quadrature(9);refnorm=0.;interpolation=0.;minlocal=1.
    for cell,(i,j,ids) in enumerate(cells):
        _,J,B=element.evaluate(coords[ids].reshape(144),h,points);weight=w*np.linalg.det(J);base,_=element.bilinear(points,xy[ids][[0,3,15,12]]);star=np.einsum('kij,pij->pk',element.MB,reference.strain('quartic',.4,base[:,0],h*points[:,2]/2));W=(np.einsum('ij,pjk->pik',element.CH,B)*np.sqrt(weight)[:,None,None]).reshape(-1,144)/normalization;e=(star@element.CH.T*np.sqrt(weight)[:,None]).reshape(-1)/normalization;refnorm+=e@e;colnorm[ids]+=np.sum(W*W,axis=0).reshape(16,9);scale=1/np.linalg.norm(W[:,PERM],axis=0);R=rfactor(np.column_stack((W[:,PERM]*scale,e)));minlocal=min(minlocal,rcond(R[:36,:36]));tail=R[36:,36:144]/scale[36:];sl=slice(109*cell,109*(cell+1));b[sl]=R[36:,144]
        for a,node in enumerate(ids[BND]):
            if int(node) in lookup:
                index=lookup[int(node)];A[sl,9*index:9*index+9]=tail[:,9*a:9*a+9]
        interpolation+=np.linalg.norm(W@interp[ids].reshape(144)-e)**2;local.append(dict(ids=ids,R=R,scale=scale))
    M,area=h1_matrix(xy,cells,coords,free,9)
    return dict(xy=xy,cells=cells,coords=coords,free=free,boundary=boundary,local=local,A=A,b=b,colnorm=colnorm,M=M,area=area,normalization=normalization,reference_ratio=float(refnorm),interpolation=float(np.sqrt(interpolation/refnorm)),local_rcond=minlocal)
def reconstruct(a,ub):
    u=np.zeros((len(a['xy']),9));u[a['boundary']]=ub.reshape(-1,9)
    for cell in a['local']:
        ids,R,scale=cell['ids'],cell['R'],cell['scale'];yB=u[ids[BND]].reshape(-1)/scale[36:];yI=solve_triangular(R[:36,:36],R[:36,144]-R[:36,36:144]@yB,check_finite=False);u[ids[INS]]=(scale[:36]*yI).reshape(4,9)
    return u
def solve(a):
    A,b=a['A'],a['b'];scale=1/np.linalg.norm(A,axis=0);R=rfactor(np.column_stack((A*scale,b)));nb=A.shape[1];T=R[:nb,:nb];condition=rcond(T);y=solve_triangular(T,R[:nb,-1],check_finite=False);ub=scale*y;u=reconstruct(a,ub);linear=float(np.linalg.norm(T@y-R[:nb,-1])/np.linalg.norm(R[:nb,-1]));gradient=np.zeros_like(u);error=0.
    for cell in a['local']:
        ids,Rl,s=cell['ids'],cell['R'],cell['scale'];v=u[ids].reshape(-1)[PERM]/s;r=Rl[:,:144]@v-Rl[:,144];error+=r@r;g=Rl[:,:144].T@r/s;back=np.zeros(144);back[PERM]=g;gradient[ids]+=back.reshape(16,9)
    fs=1/np.sqrt(a['colnorm'][a['free']]);kkt=float(np.linalg.norm(fs*gradient[a['free']])/(np.sqrt(9*len(a['free']))*np.sqrt(a['reference_ratio'])));fit=dict(error=float(np.sqrt(error/a['reference_ratio'])),boundary_minimum_model=float(abs(R[nb,nb])/np.sqrt(a['reference_ratio'])),global_rcond=condition,linear_residual=linear,kkt=kkt,h1=h1(u,a['free'],a['M']))
    return u,fit,dict(triangle=R,scale=scale,boundary_solution=ub)
def observe(a,h,u,order=11):
    points,w=element.quadrature(order);total=error=upper2=0.
    for i,j,ids in a['cells']:
        _,J,B=element.evaluate(a['coords'][ids].reshape(144),h,points);weight=w*np.linalg.det(J);xy,_=element.bilinear(points,a['xy'][ids][[0,3,15,12]]);star=np.einsum('kij,pij->pk',element.MB,reference.strain('quartic',.4,xy[:,0],h*points[:,2]/2));A=(np.einsum('ij,pjk->pik',element.CH,B)*np.sqrt(weight)[:,None,None]).reshape(-1,144);b=(star@element.CH.T*np.sqrt(weight)[:,None]).reshape(-1);value,ve=rounding.dot_rows(A,u[ids].reshape(144));r=value-b;_,upper=rounding.norms(r,ve+rounding.gamma(3)*(abs(value)+abs(b)));total+=b@b;error+=r@r;upper2+=upper*upper
    margin=rounding.gamma(len(points)*6*len(a['cells'])+10);M,area=h1_matrix(a['xy'],a['cells'],a['coords'],a['free'],11);N,g,ws=surface(11);position=position_ref=0.
    for i,j,ids in a['cells']:
        xy=N@a['xy'][ids];x=xy[:,0];star=np.column_stack((-.4*x**5/15,np.zeros(len(x)),x**4/12));J=np.einsum('na,pnb->pab',a['coords'][ids,0],g);weight=ws*np.sqrt(np.linalg.det(J.transpose(0,2,1)@J));r=N@u[ids,:3]-star;position+=weight@np.sum(r*r,axis=1);position_ref+=weight@np.sum(star*star,axis=1)
    return dict(error=float(np.sqrt(error/total)),upper=float(np.sqrt(upper2*(1+margin)/(total*(1-margin)))),reference_norm=float(total),h1=h1(u,a['free'],M),position_error=float(np.sqrt(position/position_ref)),area=area)
