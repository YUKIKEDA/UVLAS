"""G-EM1 linear field constructor; no production or nonlinear material API."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import itertools
import numpy as np
from numpy.polynomial.legendre import leggauss
import check_q9_curved_strain_screen as raw
MB=raw.EB.copy();MB[3:]*=np.sqrt(2.)
def mandel(E):return np.einsum('kij,pijn->pkn',MB,E,optimize=True)
def from_mandel(e):return np.einsum('kij,kn->ijn',MB,e,optimize=True)
def source(q,width,offset,points):
    J,B=raw.kinematics(q,width,offset,points)
    E=np.concatenate((raw.tensor(B),np.broadcast_to(MB.transpose(1,2,0),(len(J),3,3,6))),axis=3)
    return J,E
def defect(state,points):
    J,_=raw.kinematics(state['q'],state['width'],state['offset'],points)
    K=raw.covariant(J,np.broadcast_to(MB.transpose(1,2,0),(len(J),3,3,6)))[:,[0,1],2,:]
    return K-np.einsum('pai,ik->pak',raw.trace.interpolation(points),state['constant_samples'],optimize=True)
def edge_fit(state,axis,sign):
    p,w=leggauss(5);points=np.array([[x,sign,z] if axis==0 else [sign,x,z] for x,z in itertools.product(p,p)]);weight=np.array([a*b for a,b in itertools.product(w,w)])
    J,E=source(state['q'],state['width'],state['offset'],points);Ae=J[:,:,[axis,2]];G=np.einsum('pia,pib->pab',Ae,Ae);L=np.linalg.cholesky(G);U=np.linalg.solve(L,Ae.transpose(0,2,1)).transpose(0,2,1)
    T=np.einsum('pia,pijn,pjb->pabn',U,E,U,optimize=True);Y=np.stack((T[:,0,0],T[:,1,1],np.sqrt(2.)*T[:,0,1]),axis=1)
    weight=weight*np.sqrt(np.linalg.det(G));weight=np.sqrt(weight/weight.sum());Y=(weight[:,None,None]*Y).reshape(75,87);D=Y[:,81:];u,s,v=np.linalg.svd(D,full_matrices=False);ranks=[int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]];rank=ranks[1]
    if rank<3:raise ValueError('insufficient_edge_observability')
    coefficients=(v[:rank].T/s[:rank])@(u[:,:rank].T@Y);null=np.eye(6)-v[:rank].T@v[:rank]
    test=np.array([[x,sign,z] if axis==0 else [sign,x,z] for x,z in itertools.product(leggauss(9)[0],[-.8,.2,.7])]);K=defect(state,test)[:,axis,:];Jt,_=raw.kinematics(state['q'],state['width'],state['offset'],test)
    base=raw.covariant(Jt,np.broadcast_to(MB.transpose(1,2,0),(len(Jt),3,3,6)))[:,axis,2,:]
    missing=float(np.linalg.norm(np.einsum('pk,kj->pj',K,null))/max(np.linalg.norm(base),np.finfo(float).tiny))
    return dict(coefficients=coefficients,ranks=ranks,condition=float(s[0]/s[rank-1]),missing_mode_defect=missing,singular_values=s)
def r1(state,points):
    points=np.asarray(points);J,E=source(state['q'],state['width'],state['offset'],points);cov=raw.covariant(J,E);S=raw.trace.interpolation(points);replacement=np.einsum('pai,in->pan',S,state['samples'],optimize=True)
    replacement+=np.einsum('pak,kn->pan',defect(state,points),state['mean'],optimize=True)
    for axis in [0,1]:
        for sign in [-1,1]:
            edge=points.copy();edge[:,1-axis]=sign;K=defect(state,edge)[:,axis,:];c=state['edges'][axis,sign]['coefficients']-state['mean'];blend=(1+sign*points[:,1-axis])/2
            replacement[:,axis,:]+=blend[:,None]*np.einsum('pk,kn->pn',K,c,optimize=True)
    cov[:,0,2]=cov[:,2,0]=replacement[:,0];cov[:,1,2]=cov[:,2,1]=replacement[:,1];inverse=np.linalg.inv(J)
    result=np.einsum('pai,pabn,pbj->pijn',inverse,cov,inverse,optimize=True)
    return J,E,result
def prepare(q,width,offset,order):
    p,w=leggauss(order);points=np.array(list(itertools.product(p,p,p)));weights=np.array([a*b*c for a,b,c in itertools.product(w,w,w)]);ties,components=raw.trace.sampling_points();Js,Es=source(q,width,offset,ties);cov=raw.covariant(Js,Es)
    state=dict(q=q,width=width,offset=offset,samples=np.array([cov[i,c,2] for i,c in enumerate(components)]),constant_samples=np.array([cov[i,c,2,81:] for i,c in enumerate(components)]),edges={})
    J,E=source(q,width,offset,points);measure=weights*np.linalg.det(J);V=measure.sum();state['mean']=np.einsum('p,pkn->kn',measure,mandel(E))/V
    for axis,sign in itertools.product([0,1],[-1,1]):state['edges'][axis,sign]=edge_fit(state,axis,sign)
    _,_,first=r1(state,points);state['moment']=np.einsum('p,pkn->kn',measure,mandel(first-E))/V;b0=(1-points[:,0]**2)*(1-points[:,1]**2);state['bubble_mean']=float(measure@b0/V)
    if state['bubble_mean']<=0:raise ValueError('nonpositive_bubble_integral')
    return state
def evaluate(state,points):
    points=np.asarray(points);J,E,first=r1(state,points);bubble=(1-points[:,0]**2)*(1-points[:,1]**2)/state['bubble_mean'];result=first-bubble[:,None,None,None]*from_mandel(state['moment'])[None]
    return J,mandel(E),mandel(result)
