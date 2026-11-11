"""Independent M-Q3 trial; not a published Regge/MITC element implementation."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import itertools
import numpy as np
from scipy.linalg import svd,svdvals
import q3_curved_element as element
DM=39;DG=78;DIM=DM+DG
def lift(J,cov):
    inv=np.linalg.inv(J);tensor=np.einsum('pai,pabn,pbj->pijn',inv,cov,inv,optimize=True);return np.einsum('kij,pijn->pkn',element.MB,tensor,optimize=True)
def data(q,h,points):
    points=np.asarray(points);_,J,B=element.evaluate(q,h,points);J=J.copy();J[:,:,2]*=2/h;E=np.concatenate((B,np.broadcast_to(np.eye(6),(len(points),6,6))),axis=2);mid=points.copy();mid[:,2]=0.;_,J0,B0=element.evaluate(q,h,mid);J0=J0.copy();J0[:,:,2]*=2/h;E0=np.concatenate((B0,np.broadcast_to(np.eye(6),(len(points),6,6))),axis=2)
    HC=np.einsum('pia,kij,pjb->pabk',J,element.MB,J,optimize=True);HC0=np.einsum('pia,kij,pjb->pabk',J0,element.MB,J0,optimize=True);gm=np.stack((HC0[:,0,0],HC0[:,1,1],HC0[:,0,1]),axis=1);gs=np.stack((HC[:,0,2],HC[:,1,2]),axis=1);m=np.einsum('pak,pkn->pan',gm,E0,optimize=True);g=np.einsum('pak,pkn->pan',gs,E,optimize=True)
    ls=np.polynomial.legendre.legvander(points[:,0],3);lt=np.polynomial.legendre.legvander(points[:,1],3);lz=np.polynomial.legendre.legvander(points[:,2],2);M=np.zeros((len(points),3,DIM));S=np.zeros((len(points),2,DIM));column=0
    for component,imax,jmax in [(0,2,3),(1,3,2),(2,2,2)]:
        for i,j in itertools.product(range(imax+1),range(jmax+1)):M[:,component,column]=ls[:,i]*lt[:,j];column+=1
    M[:,:,33:39]=gm;column=DM
    for component,imax,jmax in [(0,2,3),(1,3,2)]:
        for i,j,k in itertools.product(range(imax+1),range(jmax+1),range(3)):S[:,component,column]=ls[:,i]*lt[:,j]*lz[:,k];column+=1
    S[:,:,DM+72:]=gs
    basis=np.zeros((len(points),3,3,DIM));basis[:,0,0]=M[:,0];basis[:,1,1]=M[:,1];basis[:,0,1]=basis[:,1,0]=M[:,2];basis[:,0,2]=basis[:,2,0]=S[:,0];basis[:,1,2]=basis[:,2,1]=S[:,1]
    raw=np.zeros((len(points),3,3,150));raw[:,0,0]=m[:,0];raw[:,1,1]=m[:,1];raw[:,0,1]=raw[:,1,0]=m[:,2];raw[:,0,2]=raw[:,2,0]=g[:,0];raw[:,1,2]=raw[:,2,1]=g[:,1]
    return dict(J=J,E=E,W=lift(J,basis),F=lift(J,raw),M=M,S=S,m=m,g=g)
def constraints(q,h,volume,weight):
    p,w=np.polynomial.legendre.leggauss(9);L=np.polynomial.legendre.legvander(p,2);rows=[];rhs=[]
    for tangent,sign in itertools.product([0,1],[-1.,1.]):
        points=np.zeros((9,3));points[:,tangent]=p;points[:,1-tangent]=sign;d=data(q,h,points)
        for i in range(3):rows.append(np.einsum('p,pn->n',w*L[:,i],d['M'][:,tangent]));rhs.append(np.einsum('p,pn->n',w*L[:,i],d['m'][:,tangent]))
        face=np.zeros((81,3));face[:,tangent]=np.repeat(p,9);face[:,1-tangent]=sign;face[:,2]=np.tile(p,9);d=data(q,h,face)
        for i,j in itertools.product(range(3),repeat=2):
            wf=np.outer(w*L[:,i],w*L[:,j]).reshape(-1);rows.append(np.einsum('p,pn->n',wf,d['S'][:,tangent]));rhs.append(np.einsum('p,pn->n',wf,d['g'][:,tangent]))
    rows.extend(np.einsum('p,pkn->kn',weight,volume['W']));rhs.extend(np.einsum('p,pkn->kn',weight,volume['F']));return np.array(rows),np.array(rhs)
def spectrum_ranks(s):return [int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]]
def prepare(q,h):
    points,w=element.quadrature(9);d=data(q,h,points);measure=w*np.linalg.det(d['J'])*h/2;weight=measure/measure.sum();W=(np.einsum('kl,pln->pkn',element.CH,d['W'])*np.sqrt(weight)[:,None,None]).reshape(-1,DIM);norm=np.linalg.norm(W,axis=0);active=norm>1e-14*norm.max();U,s,V=svd(W[:,active]/norm[active],full_matrices=False,check_finite=False);ranks=spectrum_ranks(s);r=ranks[1];T=np.zeros((DIM,r));T[active]=(V[:r].T/s[:r])/norm[active,None];del U;AW=W@T;Y=(np.einsum('kl,pln->pkn',element.CH,d['F'])*np.sqrt(weight)[:,None,None]).reshape(-1,150);coef0=AW.T@Y;F,g=constraints(q,h,d,weight);F=F@T;rownorm=np.linalg.norm(F,axis=1);rowactive=rownorm>1e-13*rownorm.max();Fn=F[rowactive]/rownorm[rowactive,None];gn=g[rowactive]/rownorm[rowactive,None];Uc,sc,Vc=svd(Fn,full_matrices=False,check_finite=False);cranks=spectrum_ranks(sc);rc=cranks[1];res=gn-Fn@coef0;incompatible=res-Uc[:,:rc]@(Uc[:,:rc].T@res);compat=float(np.linalg.norm(incompatible)/max(1.,np.linalg.norm(gn)));zero_rhs=float(np.linalg.norm(g[~rowactive])/max(1.,np.linalg.norm(g)));coef=coef0+(Vc[:rc].T/sc[:rc])@(Uc[:,:rc].T@res);coeff=T@coef;R=d['E']-d['F']+np.einsum('pkn,nj->pkj',d['W'],coeff,optimize=True);constraint_error=float(np.linalg.norm(F@coef-g)/max(1.,np.linalg.norm(g)));constant=float(np.max(abs(R[:,:,144:]-np.eye(6))));work_raw=np.einsum('p,pkn->kn',weight,d['E'][:,:,:144]);work_mod=np.einsum('p,pkn->kn',weight,R[:,:,:144]);work=float(np.linalg.norm(work_mod-work_raw)/np.linalg.norm(work_raw));raw=(np.einsum('kl,pln->pkn',element.CH,d['E'][:,:,:144])*np.sqrt(weight)[:,None,None]).reshape(-1,144);mod=(np.einsum('kl,pln->pkn',element.CH,R[:,:,:144])*np.sqrt(weight)[:,None,None]).reshape(-1,144);scale=1/np.linalg.norm(raw,axis=0);sr=svdvals(mod*scale,check_finite=False);branks=spectrum_ranks(sr);rigid=0.
    for axis in np.eye(3):
        v=np.zeros((16,3,3));v[:,0]=axis
        for mode in [v.reshape(144),np.cross(np.broadcast_to(axis,(16,3,3)),q.reshape(16,3,3)).reshape(144)]:rigid=max(rigid,float(np.linalg.norm((mod*scale)@(mode/scale))/(sr[0]*np.linalg.norm(mode/scale))))
    checks=dict(space_rank=ranks==[r]*3,constraint_rank=cranks==[rc]*3,compatible=compat<=1e-10,zero_row=zero_rhs<=1e-10,constant=constant<=1e-10,work=work<=1e-10,moments=constraint_error<=1e-10,rigid=rigid<=1e-10,local_rank=branks==[138]*3)
    return dict(q=q,h=h,coeff=coeff,space_transform=T,space_singular=s,constraint_singular=sc,B_singular=sr,constraint_matrix=F,constraint_rhs=g,log=dict(space_ranks=ranks,constraint_ranks=cranks,B_ranks=branks,compatibility=compat,zero_row_rhs=zero_rhs,constant_error=constant,work_error=work,moment_error=constraint_error,rigid_error=rigid,checks={key:bool(v) for key,v in checks.items()}))
def evaluate(state,points):
    d=data(state['q'],state['h'],points);return d['E']-d['F']+np.einsum('pkn,nj->pkj',d['W'],state['coeff'],optimize=True)
