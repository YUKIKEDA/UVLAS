"""RW-1 research kernel: reference work, six faces, geometric Piola completion."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import itertools
import numpy as np
from scipy.linalg import svd,qr,solve_triangular
import hybrid_gradient_ancf as old
import polynomial_trace_ancf as poly
m=old.m;qrtools=old.qrtools;ranks=old.ranks;volume_points=old.volume_points
FACES=old.FACES+[(2,-1),(2,1)]
VID=list(itertools.product(range(5),range(5),range(3)))
PID=poly.POT_INDICES
WID=[i for i,v in enumerate(PID) if i and v[:2]!=(5,5)]
RTID=[(a,v) for a,degs in enumerate([(5,4,2),(4,5,2),(4,4,3)]) for v in itertools.product(*(range(d+1) for d in degs))]

def basis(p,ids,independent=False):
    vs=[];ds=[]
    for a in range(3):
        deg=max(v[a] for v in ids);vs.append(np.polynomial.legendre.legvander(p[:,a],deg));dd=[]
        for j in range(deg+1):
            c=np.eye(deg+1)[j]
            if independent:
                cc=np.polynomial.legendre.leg2poly(c);dd.append(np.polynomial.polynomial.polyval(p[:,a],np.polynomial.polynomial.polyder(cc)))
            else:dd.append(np.polynomial.legendre.legval(p[:,a],np.polynomial.legendre.legder(c)))
        ds.append(np.array(dd).T)
    val=np.stack([vs[0][:,i]*vs[1][:,j]*vs[2][:,k] for i,j,k in ids],axis=1)
    grad=np.stack([np.stack([ds[0][:,i]*vs[1][:,j]*vs[2][:,k],vs[0][:,i]*ds[1][:,j]*vs[2][:,k],vs[0][:,i]*vs[1][:,j]*ds[2][:,k]],axis=1) for i,j,k in ids],axis=2)
    return val,grad

def geometry(q,h,p,independent=False):
    phi,g,_=m.tables(p,h*p[:,2]/2,independent);J=np.einsum('na,pnb->pab',q.reshape(75,3),g);J[:,:,2]*=h/2
    det=np.linalg.det(J);inv=np.linalg.inv(J);return phi@q.reshape(75,3),J,det,inv,det[:,None,None]*inv

def rt(p,independent=False):
    ids=[v for _,v in RTID];v,g=basis(p,ids,independent);factor=np.array([np.sqrt(np.prod([2*k+1 for k in v])) for v in ids])
    S=np.zeros((len(p),3,280));div=np.zeros((len(p),280))
    for j,(a,_) in enumerate(RTID):S[:,a,j]=v[:,j]*factor[j];div[:,j]=g[:,a,j]*factor[j]
    return S,div

def classify(s):return [int(np.count_nonzero(s>t)) for t in [1e-9,1e-10,1e-11]]

def complete(base,extra,mu,independent=False):
    # Inputs base are orthonormal for the reference probability measure.
    weight=mu.reshape((-1,)+(1,)*(extra.ndim-1));scale=np.sqrt(np.sum(weight*extra**2,axis=tuple(range(extra.ndim-1))));scale=np.maximum(scale,1e-30)
    if extra.ndim==2:scale[:]=max(scale)
    C=extra/scale;bw=(base*np.sqrt(weight)).reshape(-1,base.shape[-1]);cw=(C*np.sqrt(weight)).reshape(-1,3);proj=bw.T@cw;res=C-base@proj
    W=(res*np.sqrt(mu).reshape((-1,)+(1,)*(res.ndim-1))).reshape(-1,3)
    U,s,V=svd(W,full_matrices=False,lapack_driver='gesvd' if independent else 'gesdd',check_finite=False);rr=classify(s)
    if len(set(rr))!=1:raise ValueError('geometric_rank_ambiguous:'+str(rr))
    r=rr[0];lift=V[:r].T/s[:r] if r else np.zeros((3,0))
    S=np.concatenate((base,res@lift),axis=-1)
    return S,dict(scale=scale,projection=proj,lift=lift,ranks=rr,singular=s)

def apply_complete(base,extra,c):return np.concatenate((base,(extra/c['scale']-base@c['projection'])@c['lift']),axis=-1)

def mixed(state,p,independent=False):
    S,d=rt(p,independent);*_,C=geometry(state['q'],state['h'],p,independent);c=state['mixed_completion']
    return apply_complete(S,C,c),np.column_stack((d,-d@c['projection']@c['lift']))

def face_base(p,axis,independent=False):
    ids=[v for v in VID if v[axis]==0];v,_=basis(p,ids,independent);return v*np.sqrt([np.prod([2*k+1 for k in v]) for v in ids])

def face_values(state,face,p,independent=False):
    L=face_base(p,face['axis'],independent);*_,C=geometry(state['q'],state['h'],p,independent)
    return apply_complete(L,C[:,face['axis'],:],face['completion'])

def face_points(axis,sign,order,zorder):
    axes=[a for a in range(3) if a!=axis];grids=[np.polynomial.legendre.leggauss(zorder if a==2 else order) for a in axes]
    p=np.zeros((len(grids[0][0])*len(grids[1][0]),3));p[:,axis]=sign
    p[:,axes]=np.array(list(itertools.product(*[g[0] for g in grids])))
    w=np.array([a*b for a,b in itertools.product(*[g[1] for g in grids])]);return p,w

def prepare(q,h,kind='RW1',order=13,zorder=7,independent=False):
    if kind!='RW1' or h<=0:raise ValueError('invalid_RW1_input')
    p,w=volume_points(order,zorder);mu=w/8;X,J,det,inv,C=geometry(q,h,p,independent)
    if det.min()<=0:raise ValueError('reference_orientation')
    V,_=basis(p,VID,independent);V*=np.sqrt([np.prod([2*k+1 for k in v]) for v in VID])
    S0,d0=rt(p,independent);S,c=complete(S0,C,mu,independent);div=np.column_stack((d0,-d0@c['projection']@c['lift']))
    state=dict(q=q,h=h,kind=kind,p=p,w=w,mu=mu,weight=w*det,volume=float(w@det),J=J,invJ=inv,det=det,mixed_completion=c,S=S,div=div,V=V,
        minimum_J=float(det.min()),ranks=[S.shape[2]]*3)
    faces=[];offset=75
    # Local ordering: volume, bottom/top, then the four globally shared faces.
    for axis,sign in [(2,-1),(2,1)]+old.FACES:
        fp,fw=face_points(axis,sign,order,zorder);_,Jf,df,iv,Cf=geometry(q,h,fp,independent)
        L=face_base(fp,axis,independent);M,fc=complete(L,Cf[:,axis,:],fw/4,independent)
        Vf,_=basis(fp,VID,independent);Vf*=np.sqrt([np.prod([2*k+1 for k in v]) for v in VID])
        Sf,_=mixed(state,fp,independent);area=sign*Cf[:,axis,:];length=np.linalg.norm(area,axis=1)
        faces.append(dict(axis=axis,sign=sign,points=fp,w=fw,phi=Vf,hat=M,completion=fc,offset=offset,n=M.shape[1],S=Sf,
            area=area,measure=fw*length,normal=area/length[:,None],area_norm=length))
        offset+=M.shape[1]
    ns=offset;nt=S.shape[2];D=np.zeros((nt,ns));D[:,:75]=-np.einsum('p,pr,pn->rn',mu,div,V,optimize=True)
    for face in faces:
        sl=slice(face['offset'],face['offset']+face['n']);D[:,sl]=np.einsum('p,pr,pn->rn',face['w']/8,face['sign']*face['S'][:,face['axis'],:],face['hat'],optimize=True)
    Wphi,Wg=basis(p,[PID[j] for j in WID],independent);A=np.einsum('p,pjr,pjn->rn',mu,S,Wg,optimize=True)
    U,sv,Vh=svd(A,full_matrices=False,lapack_driver='gesvd' if independent else 'gesdd',check_finite=False)
    rr=ranks(sv)
    if rr!=[104]*3:raise ValueError('potential_rank:'+str(rr))
    if independent:
        Q,T=qr(A,mode='economic',check_finite=False);H=solve_triangular(T,Q.T,check_finite=False)
    else:H=(Vh.T/sv)@U.T
    Pvol=V.T@(mu[:,None]*Wphi);body=S+(Wg-S@A)@H;Gwref=body@D;Gw=np.einsum('pjn,pjk->pkn',Gwref,inv,optimize=True)
    X0=np.zeros((ns,3));X0[:75]=V.T@(mu[:,None]*X)
    for face in faces:
        fp=face['points'];M=face['hat'];Vf=face['phi'];Wf,_=basis(fp,[PID[j] for j in WID],independent)
        E=M@(M.T@((face['w']/4)[:,None]*(Wf-Vf@Pvol)))
        jump=np.zeros((len(fp),ns));jump[:,:75]=Vf;jump[:,face['offset']:face['offset']+face['n']]=-M
        face['E']=E;face['jump_raw']=jump.copy();face['jump']=jump+E@H@D
        Xf,*_=geometry(q,h,fp,independent);X0[face['offset']:face['offset']+face['n']]=M.T@((face['w']/4)[:,None]*Xf)
    B=np.einsum('kij,pnj->pkni',m.MB,Gw.transpose(0,2,1)).reshape(len(p),6,3*ns)
    elastic=(np.einsum('ij,pjn->pin',m.CH,B)*np.sqrt(w*det)[:,None,None]).reshape(-1,3*ns)
    stabilization=np.vstack([np.kron(face['jump']*np.sqrt(2.5*face['measure'])[:,None],np.eye(3)) for face in faces])
    root=np.vstack((elastic,stabilization));F0=np.einsum('na,pjn->paj',X0,Gw)
    state.update(faces=faces,D=D,H=H,A=A,Pvol=Pvol,Sbody=body,Gw=Gw,Gwref=Gwref,B=B,root=root,elastic_rows=len(elastic),tau=2.5,X0=X0,nscalar=ns,
        ni=75+sum(f['n'] for f in faces if f['axis']==2),potential_ranks=rr,gram_error=float(np.linalg.norm(np.einsum('p,pjr,pjs->rs',mu,S,S,optimize=True)-np.eye(nt))),
        reference_gradient_error=float(np.linalg.norm(F0-np.eye(3))/np.sqrt(3*len(p))))
    return state

def nonlinear(state,current,v1,v2):
    ns=state['nscalar'];u=current.reshape(ns,3)-state['X0'];Gw=state['Gw'];w=state['weight'];F=np.eye(3)+np.einsum('na,pjn->paj',u,Gw)
    E=np.einsum('kij,pij->pk',m.MB,(F.transpose(0,2,1)@F-np.eye(3))/2);stress=E@m.C.T;St=np.einsum('kij,pk->pij',m.MB,stress);P=F@St
    B=np.einsum('kij,pai,pjn->pkna',m.MB,F,Gw).reshape(-1,6,3*ns)
    d1=np.einsum('na,pjn->paj',v1.reshape(ns,3),Gw);d2=np.einsum('na,pjn->paj',v2.reshape(ns,3),Gw)
    sec=np.einsum('kij,pij->pk',m.MB,(d1.transpose(0,2,1)@d2+d2.transpose(0,2,1)@d1)/2)
    U=np.einsum('p,pk,pk->',w,E,stress)/2;force=np.einsum('p,pjn,paj->na',w,Gw,P)
    tangent=np.einsum('p,pk,kl,pl->',w,B@v1,m.C,B@v2)+np.einsum('p,pk,pk->',w,stress,sec)
    physical_body=np.einsum('pjr,pjk->pkr',state['Sbody'],state['invJ'],optimize=True)
    bulk=np.einsum('p,pjr,paj->ar',w,physical_body,P,optimize=True);eta=np.zeros((3,104),dtype=F.dtype);jumps=[]
    for face in state['faces']:
        jump=face['jump']@u;jumps.append(jump);ww=state['tau']*face['measure'];U+=np.einsum('p,pa,pa->',ww,jump,jump)/2
        force+=np.einsum('p,pn,pa->na',ww,face['jump'],jump);tangent+=np.einsum('p,pa,pa->',ww,face['jump']@v1.reshape(ns,3),face['jump']@v2.reshape(ns,3))
        eta+=np.einsum('p,pr,pa->ar',ww,face['E'],jump)
    dual=(bulk+eta@state['H'])/8;fluxes=[]
    for face,jump in zip(state['faces'],jumps):
        tr=face['sign']*np.einsum('pr,ar->pa',face['S'][:,face['axis'],:],dual)/face['area_norm'][:,None]-state['tau']*jump;fluxes.append(tr)
    return U,force.reshape(-1),tangent,dict(min_F=float(np.min(np.linalg.det(F.real))),fluxes=fluxes)

def reduce(state,fun):
    u=np.zeros((state['nscalar'],fun(state['p']).shape[1]));u[:75]=state['V'].T@(state['mu'][:,None]*fun(state['p']))
    for face in state['faces']:u[face['offset']:face['offset']+face['n']]=face['hat'].T@((face['w']/4)[:,None]*fun(face['points']))
    return u

def evaluate_gradient(state,p,independent=True):
    S,_=mixed(state,p,independent);_,Wg=basis(p,[PID[j] for j in WID],independent);body=S+(Wg-S@state['A'])@state['H'];_,_,_,inv,_=geometry(state['q'],state['h'],p,independent)
    return np.einsum('pjn,pjk->pkn',body@state['D'],inv,optimize=True)

def boundary_map(state):
    result=np.zeros((state['nscalar'],3))
    for face in state['faces']:result[face['offset']:face['offset']+face['n']]+=np.einsum('p,pn,pj->nj',face['w'],face['hat'],face['area'])
    return result

def volume_nodal(state,coeff):
    factor=np.sqrt([np.prod([2*k+1 for k in v]) for v in VID]);return poly.nodal_volume(coeff*factor[:,None],state['h'])
