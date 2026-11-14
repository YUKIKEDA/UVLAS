"""RW-3: same RW-2 spaces, analytic modal moments and stable basis evaluation."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import svd,qr,solve_triangular
import reference_work_ancf_v2 as prev
import reference_work_coefficients as c
m=prev.m;qrtools=prev.qrtools;ranks=prev.ranks;volume_points=prev.volume_points
VID=prev.VID;PID=prev.PID;WID=prev.WID;FACES=prev.FACES;basis=prev.basis
geometry=prev.geometry;face_points=prev.face_points;boundary_map=prev.boundary_map
reduce=prev.reduce;volume_nodal=prev.volume_nodal;nonlinear=prev.nonlinear
VI=np.array([c.LOOK[v] for v in VID]);WP,WG=c.potential_coefficients()

def mixed(state,p,independent=False):
    L,_=c.scalar_basis(p,independent)
    return np.einsum('pk,kjr->pjr',L,state['Scoeff'],optimize=True),L@state['divcoeff']

def face_values(state,face,p,independent=False):
    L,_=c.scalar_basis(p,independent);return L@face['Mcoeff']

def prepare(q,h,kind='RW3',order=13,zorder=7,independent=False):
    if kind!='RW3' or h<=0:raise ValueError('invalid_RW3_input')
    p,w=volume_points(order,zorder);mu=w/8;X,J,det,inv,C=geometry(q,h,p,independent)
    if det.min()<=0:raise ValueError('reference_orientation')
    Scoef,divcoef,meta,modal=c.mixed_coefficients(q,independent);L,_=c.scalar_basis(p,independent);S=np.einsum('pk,kjr->pjr',L,Scoef,optimize=True);div=L@divcoef;V=L[:,VI]
    state=dict(q=q,h=h,kind=kind,p=p,w=w,mu=mu,weight=w*det,volume=float(w@det),J=J,invJ=inv,det=det,mixed_completion=meta,S=S,div=div,V=V,
        minimum_J=float(det.min()),ranks=[S.shape[2]]*3,Scoeff=Scoef,divcoeff=divcoef)
    faces=[];offset=75
    for axis,sign in [(2,-1),(2,1)]+prev.old.FACES:
        fp,fw=face_points(axis,sign,order,zorder);_,Jf,df,iv,Cf=geometry(q,h,fp,independent);Lf,_=c.scalar_basis(fp,independent)
        Mc,fc=c.face_coefficients(q,axis,sign,independent);M=Lf@Mc;Sf=np.einsum('pk,kjr->pjr',Lf,Scoef,optimize=True);area=sign*Cf[:,axis,:];length=np.linalg.norm(area,axis=1)
        Sn=c.restrict(Scoef[:,axis,:],axis,sign)
        faces.append(dict(axis=axis,sign=sign,points=fp,w=fw,phi=Lf[:,VI],hat=M,completion=fc,offset=offset,n=M.shape[1],S=Sf,Mcoeff=Mc,Sncoeff=Sn,
            area=area,measure=fw*length,normal=area/length[:,None],area_norm=length))
        offset+=M.shape[1]
    ns=offset;nt=S.shape[2];D=np.zeros((nt,ns));D[:,:75]=-divcoef[VI].T
    for face in faces:D[:,face['offset']:face['offset']+face['n']]=face['sign']*.5*face['Sncoeff'].T@face['Mcoeff']
    Sflat=Scoef.reshape(-1,nt);Gflat=WG.reshape(-1,104);A=Sflat.T@Gflat
    U,sv,Vh=svd(A,full_matrices=False,lapack_driver='gesvd' if independent else 'gesdd',check_finite=False);rr=ranks(sv)
    if rr!=[104]*3:raise ValueError('potential_rank:'+str(rr))
    if independent:
        Q,T=qr(A,mode='economic',check_finite=False);H=solve_triangular(T,Q.T,check_finite=False)
    else:H=(Vh.T/sv)@U.T
    residual=Gflat.copy();residual[c.RTROWS]=0;extra=Sflat[:,280:]
    if extra.shape[1]:residual-=extra@(extra.T@residual)
    bodycoef=(Sflat+residual@H).reshape(len(c.IDS),3,nt);body=np.einsum('pk,kjr->pjr',L,bodycoef,optimize=True)
    Gwref=body@D;Gw=np.einsum('pjn,pjk->pkn',Gwref,inv,optimize=True);Pvol=WP[VI]
    X0=np.zeros((ns,3));X0[:75]=c.reference_position(modal,h);Xcoef=np.zeros((len(c.IDS),3));Xcoef[VI]=X0[:75]
    volume_p=np.zeros_like(WP);volume_p[VI]=Pvol
    for face in faces:
        Mc=face['Mcoeff'];axis=face['axis'];sign=face['sign'];M=face['hat'];Vf=face['phi']
        correction=c.restrict(WP-volume_p,axis,sign);E=M@(Mc.T@correction)
        jump=np.zeros((len(M),ns));jump[:,:75]=Vf;jump[:,face['offset']:face['offset']+face['n']]=-M
        face['E']=E;face['jump_raw']=jump.copy();face['jump']=jump+E@H@D
        X0[face['offset']:face['offset']+face['n']]=Mc.T@c.restrict(Xcoef,axis,sign)
    B=np.einsum('kij,pnj->pkni',m.MB,Gw.transpose(0,2,1)).reshape(len(p),6,3*ns)
    elastic=(np.einsum('ij,pjn->pin',m.CH,B)*np.sqrt(w*det)[:,None,None]).reshape(-1,3*ns)
    stabilization=np.vstack([np.kron(face['jump']*np.sqrt(2.5*face['measure'])[:,None],np.eye(3)) for face in faces]);root=np.vstack((elastic,stabilization))
    F0=np.einsum('na,pjn->paj',X0,Gw)
    state.update(faces=faces,D=D,H=H,A=A,Pvol=Pvol,Sbody=body,bodycoeff=bodycoef,Gw=Gw,Gwref=Gwref,B=B,root=root,elastic_rows=len(elastic),tau=2.5,X0=X0,nscalar=ns,
        ni=75+sum(face['n'] for face in faces if face['axis']==2),potential_ranks=rr,gram_error=float(np.linalg.norm(Sflat.T@Sflat-np.eye(nt))),
        reference_gradient_error=float(np.linalg.norm(F0-np.eye(3))/np.sqrt(3*len(p))))
    return state

def evaluate_gradient(state,p,independent=True):
    L,_=c.scalar_basis(p,independent);body=np.einsum('pk,kjr->pjr',L,state['bodycoeff'],optimize=True);_,_,_,inv,_=geometry(state['q'],state['h'],p,independent)
    return np.einsum('pjn,pjk->pkn',body@state['D'],inv,optimize=True)
