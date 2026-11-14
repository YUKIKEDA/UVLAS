"""MMS-1 leading bending energy, all mixed q01 and shared q2."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import svd,lstsq,cholesky,solve_triangular
import mixed_mean_stress_thin_v6 as mixed
old=mixed.old;f=mixed.f;c=mixed.c;m=mixed.m;fit=mixed.fit
operators=old.operators;nullspace=old.nullspace;rank_info=old.rank_info;even_internal_basis=old.even_internal_basis
KINDS=['quartic','quintic','exponential']
def leading_rows(op,kappa):
    roots,_=old.leading_rows(op,kappa);s=op['state'];N,_=m.shape(op['p'],0);x=(N@s['q'].reshape(75,3)[:25])[:,0];a=np.column_stack((np.ones(len(x)),np.zeros(len(x)),kappa*x));targets=np.zeros((roots[0].shape[0],3))
    for j,curvature in enumerate([x*x,20*x**3,np.exp(x)]):
        E=-op['p'][:,2,None,None]/2*curvature[:,None,None]*a[:,:,None]*a[:,None,:]/(1+kappa*kappa*x*x)[:,None,None];mandel=np.einsum('kij,pij->pk',m.MB,E);targets[:len(op['p'])*6,j]=(mandel@m.CH.T*np.sqrt(op['w']*op['d0'])[:,None]).reshape(-1)
    return roots,targets
def local(q,kappa):
    cell=mixed.local(q,kappa);op=operators(q);s=op['state'];nc=cell['nc'];ni=cell['ni'];roots,b=leading_rows(op,kappa);Z=cell['Z'];W01=np.column_stack(roots[:2])@Z;Pi=even_internal_basis(s);Wi=roots[2][:,:ni]@Pi;Wb=roots[2][:,ni:]
    U,sv,V=svd(Wi,full_matrices=False,check_finite=False);r=rank_info(sv)
    if len(set(r))!=1:raise ValueError('MMS1_q2_internal_rank_ambiguous:'+str(r))
    Q=U[:,:r[0]];combined=np.column_stack((W01,Wb,b));combined-=Q@(Q.T@combined);R=f.qrtools.rfactor(np.asfortranarray(combined));cell.update(R=R,reference=np.sum(b*b,axis=0),q2_internal_ranks=r);return cell
def recover(cell,y,q2b,kappa):
    op=operators(cell['q']);roots,b=leading_rows(op,kappa);q01=cell['Z']@y;nc=cell['nc'];ni=cell['ni'];rhs=b-roots[0]@q01[:nc]-roots[1]@q01[nc:]-roots[2][:,ni:]@q2b;Pi=even_internal_basis(op['state']);q,_,rank,sv=lstsq(roots[2][:,:ni]@Pi,rhs,cond=1e-10,lapack_driver='gelsd',check_finite=False)
    if rank_info(sv)!=cell['q2_internal_ranks']:raise ValueError('MMS1_q2_recovery_rank')
    return np.stack((q01[:nc],q01[nc:],np.vstack((Pi@q,q2b))))
def geometry_difference(q,k):
    q=q.reshape(75,3);corners=q[[0,4,20,24]];xp=mixed.stress.p.bilinear(corners[:,0]);yp=mixed.stress.p.bilinear(corners[:,1]);points=np.array([[x,y] for y in m.LOBATTO for x in m.LOBATTO]);x=sum(float(v)*points[:,0]**i*points[:,1]**j for (i,j),v in xp.items());y=sum(float(v)*points[:,0]**i*points[:,1]**j for (i,j),v in yp.items());expected=np.vstack((np.column_stack((x,y,k*x*x/2)),np.column_stack((-k*x,np.zeros(25),np.ones(25))),np.zeros((25,3))));return float(np.max(abs(expected-q)))
def constraints_and_position(op,values,k):
    s=op['state'];ns=s['nscalar'];v=values.reshape(3,ns,3,3);norm=max(1.,np.linalg.norm(v));N0=op['A'][0]@values[0];zero=float(np.linalg.norm(N0)/norm);jump=0.
    for face in s['faces']:
        jump=max(jump,float(np.linalg.norm(np.einsum('pn,nak->pak',face['jump'],v[0]))/norm))
        if face['axis']==2:jump=max(jump,float(np.linalg.norm(np.einsum('pn,nak->pak',face['jump'],v[1]))/norm))
    count=len(op['p']);nz=len(np.unique(op['p'][:,2]));np2=count//nz;_,wz=np.polynomial.legendre.leggauss(nz);E=(op['A'][1]@values[0]+op['A'][0]@values[1])/op['d0'][:,None,None];E=E.reshape(np2,nz,6,3);mean=np.einsum('z,pzkt->pkt',wz/2,E);fluct=float(np.linalg.norm(E-mean[:,None])/norm);plane,wp=m.quadrature(int(np.sqrt(np2)),2);L,_=mixed.stress.evaluate(s['q'],k,plane,nz==9);weight=2*wp*op['d0'].reshape(np2,nz)[:,0];gram=np.einsum('p,pir,ij,pjs->rs',weight,L,np.linalg.inv(m.C),L,optimize=True);rhs=np.einsum('p,pkr,pkt->rt',weight,L,mean,optimize=True);work=float(np.linalg.norm(solve_triangular(cholesky(gram,lower=True),rhs,lower=True))/norm)
    Q,R,coef=mixed.scalar_position(s);a=solve_triangular(R,Q.T@v[0].reshape(ns,-1));p=np.column_stack((plane,np.zeros(len(plane))));Lc,_=c.scalar_basis(p,nz==9);position=(Lc@coef@a).reshape(np2,3,3);target=mixed.target(s['q'],p,k);N,g=m.shape(p,0,nz==9);J=np.einsum('na,pnb->pab',s['q'].reshape(75,3)[:25],g);measure=wp*np.linalg.norm(np.cross(J[:,:,0],J[:,:,1]),axis=1);err=np.einsum('p,pak,pak->k',measure,position-target,position-target);ref=np.einsum('p,pak,pak->k',measure,target,target)
    return dict(zero=zero,jump=jump,work=work,fluctuation=fluct),err,ref
def observe(cells,values,k,independent=False):
    elastic=np.zeros(3);vertical=np.zeros(3);reference=np.zeros(3);poserr=np.zeros(3);posref=np.zeros(3);constraints=[]
    for cell,v in zip(cells,values):
        op=operators(cell['q'],15 if independent else 13,9 if independent else 7,independent);roots,b=leading_rows(op,k);res=roots[0]@v[0]+roots[1]@v[1]+roots[2]@v[2]-b;nr=6*len(op['p']);elastic+=np.sum(res[:nr]**2,axis=0);vertical+=np.sum(res[nr:]**2,axis=0);reference+=np.sum(b*b,axis=0);check,e,r=constraints_and_position(op,v,k);constraints.append(check);poserr+=e;posref+=r
    return dict(error=np.sqrt((elastic+vertical)/reference).tolist(),elastic=np.sqrt(elastic/reference).tolist(),stabilization=np.sqrt(vertical/reference).tolist(),position_error=np.sqrt(poserr/posref).tolist(),constraints={key:max(x[key] for x in constraints) for key in constraints[0]})
