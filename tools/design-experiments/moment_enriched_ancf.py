"""ME-1: one volume position moment and two zero-normal reference fluxes."""
import numpy as np
from scipy.linalg import svd,qr,solve_triangular
import reference_work_ancf_v3 as old
f=old;m=old.m;c=old.c;qrtools=old.qrtools;ranks=old.ranks;volume_points=old.volume_points;geometry=old.geometry;basis=old.basis;PID=old.PID
VID=old.VID+[(5,5,0)];VI=np.array([c.LOOK[v] for v in VID]);WID=old.WID+[PID.index((5,5,0))];NV=76
WP0,WG0=c.potential_coefficients();WP=np.zeros((len(c.IDS),105));WG=np.zeros((len(c.IDS),3,105));WP[:,:104]=WP0;WG[:,:,:104]=WG0
WP[c.LOOK[5,5,0],104]=1/11
for axis in [0,1]:
    for k in [4,2,0]:
        key=[5,5,0];key[axis]=k;WG[c.LOOK[tuple(key)],axis,104]=np.sqrt((2*5+1)*(2*k+1))/11

def extended_mixed(q,independent=False):
    S,_,meta,modal=c.mixed_coefficients(q,independent);n=S.shape[-1];B=np.zeros((len(c.IDS),3,2))
    for axis in [0,1]:
        for degree,sign in [(6,1),(4,-1)]:
            key=[5,5,0];key[axis]=degree;B[c.LOOK[tuple(key)],axis,axis]=sign/(11*c.NORM[c.LOOK[tuple(key)]])
    flat=S.reshape(-1,n);res=B.reshape(-1,2);scale=np.max(np.linalg.norm(res,axis=0));res/=scale;res[c.RTROWS]=0
    if n>280:res-=flat[:,280:]@(flat[:,280:].T@res)
    U,sv,V=svd(res,full_matrices=False,lapack_driver='gesvd' if independent else 'gesdd',check_finite=False);rr=c.parent.classify(sv)
    if rr!=[2]*3:raise ValueError('bubble_flux_pair_rank:'+str(rr))
    extra=(res@V.T)/sv;coeff=np.concatenate((S,extra.reshape(len(c.IDS),3,2)),axis=2)
    return coeff,c.derivative(coeff),dict(geometry=meta,bubble_ranks=rr),modal

def prepare(q,h,kind='ME1',order=13,zorder=7,independent=False):
    if h<=0 or kind!='ME1':raise ValueError('invalid_ME1_input')
    p,w=volume_points(order,zorder);mu=w/8;_,J,det,inv,_=geometry(q,h,p,independent)
    if det.min()<=0:raise ValueError('reference_orientation')
    Scoef,divcoef,meta,modal=extended_mixed(q,independent);L,_=c.scalar_basis(p,independent);S=np.einsum('pk,kjr->pjr',L,Scoef,optimize=True);V=L[:,VI];div=L@divcoef
    state=dict(q=q,h=h,kind=kind,p=p,w=w,mu=mu,weight=w*det,volume=float(w@det),J=J,invJ=inv,det=det,S=S,Scoeff=Scoef,div=div,divcoeff=divcoef,V=V,minimum_J=float(det.min()),mixed_completion=meta,ranks=[S.shape[-1]]*3)
    faces=[];offset=NV
    for axis,sign in [(2,-1),(2,1)]+old.prev.old.FACES:
        fp,fw=old.face_points(axis,sign,order,zorder);_,_,_,_,Cf=geometry(q,h,fp,independent);Lf,_=c.scalar_basis(fp,independent);Mc,fc=c.face_coefficients(q,axis,sign,independent);M=Lf@Mc
        Sf=np.einsum('pk,kjr->pjr',Lf,Scoef,optimize=True);area=sign*Cf[:,axis,:];length=np.linalg.norm(area,axis=1);Sn=c.restrict(Scoef[:,axis,:],axis,sign)
        faces.append(dict(axis=axis,sign=sign,points=fp,w=fw,phi=Lf[:,VI],hat=M,completion=fc,offset=offset,n=M.shape[1],S=Sf,Mcoeff=Mc,Sncoeff=Sn,area=area,measure=fw*length,normal=area/length[:,None],area_norm=length));offset+=M.shape[1]
    ns=offset;nt=S.shape[-1];D=np.zeros((nt,ns));D[:,:NV]=-divcoef[VI].T
    for face in faces:D[:,face['offset']:face['offset']+face['n']]=face['sign']*.5*face['Sncoeff'].T@face['Mcoeff']
    flat=Scoef.reshape(-1,nt);A=flat.T@WG.reshape(-1,105);U,sv,Vh=svd(A,full_matrices=False,lapack_driver='gesvd' if independent else 'gesdd',check_finite=False);rr=ranks(sv)
    if rr!=[105]*3:raise ValueError('potential_rank:'+str(rr))
    if independent:
        Q,T=qr(A,mode='economic',check_finite=False);H=solve_triangular(T,Q.T,check_finite=False)
    else:H=(Vh.T/sv)@U.T
    residual=WG.reshape(-1,105).copy();residual[c.RTROWS]=0;extra=flat[:,280:];residual-=extra@(extra.T@residual)
    bodycoef=(flat+residual@H).reshape(len(c.IDS),3,nt);body=np.einsum('pk,kjr->pjr',L,bodycoef,optimize=True);Gwref=body@D;Gw=np.einsum('pjn,pjk->pkn',Gwref,inv,optimize=True)
    Pvol=WP[VI];Vcoeff=np.zeros((len(c.IDS),NV));Vcoeff[VI,np.arange(NV)]=1;vp=np.zeros_like(WP);vp[VI]=Pvol
    X0=np.zeros((ns,3));X0[:75]=c.reference_position(modal,h);Xcoef=np.zeros((len(c.IDS),3));Xcoef[VI]=X0[:NV]
    for face in faces:
        axis,sign=face['axis'],face['sign'];Mc=face['Mcoeff'];M=face['hat'];E=M@(Mc.T@c.restrict(WP-vp,axis,sign));projectedV=M@(Mc.T@c.restrict(Vcoeff,axis,sign))
        raw=np.zeros((len(M),ns));raw[:,:NV]=projectedV;raw[:,face['offset']:face['offset']+face['n']]=-M
        face['E']=E;face['jump_raw']=raw.copy();face['jump']=raw+E@H@D;X0[face['offset']:face['offset']+face['n']]=Mc.T@c.restrict(Xcoef,axis,sign)
    B=np.einsum('kij,pnj->pkni',m.MB,Gw.transpose(0,2,1)).reshape(len(p),6,3*ns);elastic=(np.einsum('ij,pjn->pin',m.CH,B)*np.sqrt(w*det)[:,None,None]).reshape(-1,3*ns)
    stab=np.vstack([np.kron(face['jump']*np.sqrt(2.5*face['measure'])[:,None],np.eye(3)) for face in faces]);F0=np.einsum('na,pjn->paj',X0,Gw)
    state.update(faces=faces,D=D,H=H,A=A,Pvol=Pvol,Sbody=body,bodycoeff=bodycoef,Gw=Gw,Gwref=Gwref,B=B,root=np.vstack((elastic,stab)),elastic_rows=len(elastic),tau=2.5,X0=X0,nscalar=ns,ni=NV+sum(x['n'] for x in faces[:2]),
        potential_ranks=rr,gram_error=float(np.linalg.norm(flat.T@flat-np.eye(nt))),reference_gradient_error=float(np.linalg.norm(F0-np.eye(3))/np.sqrt(3*len(p))))
    return state

def reduce(state,fun):
    u=np.zeros((state['nscalar'],fun(state['p']).shape[1]));u[:NV]=state['V'].T@(state['mu'][:,None]*fun(state['p']))
    for face in state['faces']:u[face['offset']:face['offset']+face['n']]=face['hat'].T@((face['w']/4)[:,None]*fun(face['points']))
    return u

def nonlinear(state,current,v1,v2):
    # Same exact U/f/K and face dual as RW-3, with the potential dimension dynamic.
    ns=state['nscalar'];u=current.reshape(ns,3)-state['X0'];Gw=state['Gw'];w=state['weight'];F=np.eye(3)+np.einsum('na,pjn->paj',u,Gw)
    E=np.einsum('kij,pij->pk',m.MB,(F.transpose(0,2,1)@F-np.eye(3))/2);stress=E@m.C.T;St=np.einsum('kij,pk->pij',m.MB,stress);P=F@St;B=np.einsum('kij,pai,pjn->pkna',m.MB,F,Gw).reshape(-1,6,3*ns)
    d1=np.einsum('na,pjn->paj',v1.reshape(ns,3),Gw);d2=np.einsum('na,pjn->paj',v2.reshape(ns,3),Gw);sec=np.einsum('kij,pij->pk',m.MB,(d1.transpose(0,2,1)@d2+d2.transpose(0,2,1)@d1)/2)
    U=np.einsum('p,pk,pk->',w,E,stress)/2;force=np.einsum('p,pjn,paj->na',w,Gw,P);K=np.einsum('p,pk,kl,pl->',w,B@v1,m.C,B@v2)+np.einsum('p,pk,pk->',w,stress,sec)
    physical=np.einsum('pjr,pjk->pkr',state['Sbody'],state['invJ'],optimize=True);bulk=np.einsum('p,pjr,paj->ar',w,physical,P,optimize=True);eta=np.zeros((3,state['H'].shape[0]),dtype=F.dtype);jumps=[]
    for face in state['faces']:
        jump=face['jump']@u;jumps.append(jump);ww=state['tau']*face['measure'];U+=np.einsum('p,pa,pa->',ww,jump,jump)/2;force+=np.einsum('p,pn,pa->na',ww,face['jump'],jump)
        K+=np.einsum('p,pa,pa->',ww,face['jump']@v1.reshape(ns,3),face['jump']@v2.reshape(ns,3));eta+=np.einsum('p,pr,pa->ar',ww,face['E'],jump)
    dual=(bulk+eta@state['H'])/8;fluxes=[face['sign']*np.einsum('pr,ar->pa',face['S'][:,face['axis'],:],dual)/face['area_norm'][:,None]-state['tau']*jump for face,jump in zip(state['faces'],jumps)]
    return U,force.reshape(-1),K,dict(min_F=float(np.min(np.linalg.det(F.real))),fluxes=fluxes)

boundary_map=old.boundary_map
def evaluate_gradient(state,p,independent=True):
    L,_=c.scalar_basis(p,independent);body=np.einsum('pk,kjr->pjr',L,state['bodycoeff'],optimize=True);_,_,_,inv,_=geometry(state['q'],state['h'],p,independent)
    return np.einsum('pjn,pjk->pkn',body@state['D'],inv,optimize=True)
