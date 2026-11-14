"""Reference-fixed mixed weak gradient and shared position trace for ANCF EL2."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import itertools
import numpy as np
from scipy.linalg import svd,qr,solve_triangular
import qem_el2_q4 as m
import q3_quartic_approximation as qrtools
FACES=[(0,-1),(0,1),(1,-1),(1,1)]

def ranks(s):
    if not len(s) or not s[0]:return [0,0,0]
    return [int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]]

def volume_points(order,zorder):
    plane,wp=m.quadrature(order,2);z,wz=np.polynomial.legendre.leggauss(zorder)
    return np.array([[s,t,c] for (s,t),c in itertools.product(plane,z)]),np.array([a*b for a,b in itertools.product(wp,wz)])

def edge_nodes(axis,sign):
    indices=[j*5+(0 if sign<0 else 4) for j in range(5)] if axis==0 else [(0 if sign<0 else 20)+i for i in range(5)]
    return np.array([25*f+i for f in range(3) for i in indices])

def reference_state(q):return np.vstack((q.reshape(75,3),*[q.reshape(75,3)[edge_nodes(*face)] for face in FACES]))

def generators(q,h,p,kind,independent=False):
    degree=4 if kind=='S4' else 3
    tables=[np.polynomial.legendre.legvander(p[:,i],degree if i<2 else 2) for i in range(3)]
    phi=np.einsum('pi,pj,pk->pijk',*tables).reshape(len(p),-1);ns=phi.shape[1]
    A=np.zeros((len(p),3,3*ns+(3 if kind=='S3+' else 0)))
    for j in range(3):A[:,j,j*ns:(j+1)*ns]=phi
    _,g,_=m.tables(p,h*p[:,2]/2,independent);J=np.einsum('na,pnb->pab',q.reshape(75,3),g);G=g@np.linalg.inv(J)
    if kind=='S3+':
        L=np.polynomial.legendre.legval(m.LOBATTO,[0,0,0,0,1]);highest=np.outer(L,L).reshape(-1)
        for f,scale in enumerate([1.,2/h,8/h**2]):A[:,:,-3+f]=np.einsum('n,pnj->pj',highest,G[:,m.SLICES[f],:])*scale
    return A,J,G

def face_data(q,h,axis,sign,order,zorder,independent=False):
    t,wt=np.polynomial.legendre.leggauss(order);z,wz=np.polynomial.legendre.leggauss(zorder)
    p=np.zeros((order*zorder,3));p[:,axis]=sign;p[:,1-axis]=np.repeat(t,zorder);p[:,2]=np.tile(z,order)
    w=np.outer(wt,wz).reshape(-1);phi,J,_=m.raw(q,p,h*p[:,2]/2,independent)
    area=sign*np.cross(J[:,:,(axis+1)%3],J[:,:,(axis+2)%3])*h/2
    return p,w,phi,area

def prepare(q,h,kind,order=13,zorder=7,independent=False):
    if h<=0:raise ValueError('positive_physical_thickness_required')
    p,w=volume_points(order,zorder);A,J,G=generators(q,h,p,kind,independent);det=np.linalg.det(J)
    if np.min(det)<=0:raise ValueError('invalid_reference_jacobian')
    weight=w*det*h/2;V0=weight.sum();mu=weight/V0
    W=(A*np.sqrt(mu)[:,None,None]).reshape(-1,A.shape[2]);norm=np.linalg.norm(W,axis=0);Wn=W/norm
    if independent:
        s=svd(Wn,compute_uv=False,lapack_driver='gesvd',check_finite=False);U,R,pivot=qr(Wn,mode='economic',pivoting=True,check_finite=False)
        coeff=np.zeros((A.shape[2],A.shape[2]));coeff[pivot]=solve_triangular(R,np.eye(len(R)),check_finite=False);coeff/=norm[:,None]
    else:
        U,s,V=svd(Wn,full_matrices=False,check_finite=False);coeff=(V.T/s)/norm[:,None]
    rr=ranks(s)
    if rr!=[A.shape[2]]*3:raise ValueError('mixed_space_rank_failed:'+str(rr)+' expected '+str(A.shape[2]))
    S=A@coeff;gram=np.einsum('p,pjr,pjs->rs',mu,S,S,optimize=True)
    residual=float(np.linalg.norm(W-U@(U.T@W))/np.linalg.norm(W));gram_error=float(np.linalg.norm(gram-np.eye(len(gram))))
    if max(residual,gram_error)>1e-10:raise ValueError('mixed_space_gram_failed:'+str((residual,gram_error)))
    D=np.zeros((S.shape[2],135));D[:,:75]=np.einsum('p,pjr,pnj->rn',mu,S,G,optimize=True);faces=[]
    for face,(axis,sign) in enumerate(FACES):
        fp,fw,phi,area=face_data(q,h,axis,sign,order,zorder,independent);Af,_,_=generators(q,h,fp,kind,independent);Sf=Af@coeff
        normal_length=np.linalg.norm(area,axis=1);traction=np.einsum('pjr,pj->pr',Sf,area)
        D[:,:75]-=np.einsum('p,pr,pn->rn',fw,traction,phi,optimize=True)/V0
        hat=phi[:,edge_nodes(axis,sign)];offset=75+15*face
        D[:,offset:offset+15]+=np.einsum('p,pr,pn->rn',fw,traction,hat,optimize=True)/V0
        jump=np.zeros((len(fp),135));jump[:,:75]=phi;jump[:,offset:offset+15]=-hat
        faces.append(dict(points=fp,w=fw,phi=phi,hat=hat,area=area,measure=fw*normal_length,normal=area/normal_length[:,None],S=Sf,jump=jump,offset=offset))
    Gw=S@D;B=np.einsum('kij,pnj->pkni',m.MB,Gw.transpose(0,2,1)).reshape(len(p),6,405)
    elastic=(np.einsum('ij,pjn->pin',m.CH,B)*np.sqrt(weight)[:,None,None]).reshape(-1,405)
    tau=2.5;stabilization=np.vstack([np.kron(face['jump']*np.sqrt(tau*face['measure'])[:,None],np.eye(3)) for face in faces])
    root=np.vstack((elastic,stabilization));X0=reference_state(q);F0=np.einsum('na,pjn->paj',X0,Gw)
    return dict(q=q,h=h,kind=kind,p=p,w=w,J=J,weight=weight,mu=mu,volume=V0,S=S,coeff=coeff,D=D,Gw=Gw,faces=faces,
        B=B,root=root,elastic_rows=len(elastic),tau=tau,X0=X0,ranks=rr,gram_error=gram_error,span_error=residual,
        reference_gradient_error=float(np.linalg.norm(F0-np.eye(3))/np.sqrt(3*len(p))),minimum_J=float(det.min()))

def nonlinear(state,current,v1,v2):
    displacement=current.reshape(135,3)-state['X0'];Gw=state['Gw'];weight=state['weight']
    F=np.eye(3)+np.einsum('na,pjn->paj',displacement,Gw);E=np.einsum('kij,pij->pk',m.MB,(F.transpose(0,2,1)@F-np.eye(3))/2)
    stress=E@m.C.T;Stensor=np.einsum('kij,pk->pij',m.MB,stress);P=F@Stensor
    B=np.einsum('kij,pai,pjn->pkna',m.MB,F,Gw).reshape(-1,6,405)
    dF1=np.einsum('na,pjn->paj',v1.reshape(135,3),Gw);dF2=np.einsum('na,pjn->paj',v2.reshape(135,3),Gw)
    second=np.einsum('kij,pij->pk',m.MB,(dF1.transpose(0,2,1)@dF2+dF2.transpose(0,2,1)@dF1)/2)
    U=np.einsum('p,pk,pk->',weight,E,stress)/2;force=np.einsum('p,pjn,paj->na',weight,Gw,P)
    tangent=np.einsum('p,pk,kl,pl->',weight,B@v1,m.C,B@v2)+np.einsum('p,pk,pk->',weight,stress,second)
    pcoef=np.einsum('p,pjr,paj->ar',state['mu'],state['S'],P,optimize=True);fluxes=[]
    for face in state['faces']:
        jump=face['jump']@displacement;weighted=state['tau']*face['measure'];U+=np.einsum('p,pa,pa->',weighted,jump,jump)/2
        force+=np.einsum('p,pn,pa->na',weighted,face['jump'],jump)
        tangent+=np.einsum('p,pa,pa->',weighted,face['jump']@v1.reshape(135,3),face['jump']@v2.reshape(135,3))
        Pface=np.einsum('pjr,ar->paj',face['S'],pcoef);fluxes.append(np.einsum('paj,pj->pa',Pface,face['normal'])-state['tau']*jump)
    return U,force.reshape(-1),tangent,dict(min_F=float(np.min(np.linalg.det(F.real))),fluxes=fluxes)

def boundary_map(state,order=None,zorder=None,independent=False):
    q,h=state['q'],state['h'];order=order or int(np.sqrt(len(state['p'])/len(np.unique(state['p'][:,2]))));zorder=zorder or len(np.unique(state['p'][:,2]))
    exterior=np.zeros((135,3))
    for face in state['faces']:
        exterior[face['offset']:face['offset']+15]+=np.einsum('p,pn,pj->nj',face['w'],face['hat'],face['area'])
    plane,w=m.quadrature(order,2)
    for sign in [-1,1]:
        p=np.column_stack((plane,np.full(len(plane),sign)));phi,J,_=m.raw(q,p,h*p[:,2]/2,independent);area=sign*np.cross(J[:,:,0],J[:,:,1])
        exterior[:75]+=np.einsum('p,pn,pj->nj',w,phi,area)
    return exterior

def mesh(k,n,amp):
    original=m.mesh(k,n,amp);edge_lookup={};connect=[];trace_ref=[];fixed=[]
    for cell,ids in enumerate(original['cells']):
        j,i=divmod(cell,n);keys=[(0,i,j),(0,i+1,j),(1,j,i),(1,j+1,i)];local=[]
        q=original['ref'][ids]
        for face,key in enumerate(keys):
            value=q[edge_nodes(*FACES[face])]
            if key not in edge_lookup:
                edge_lookup[key]=len(trace_ref);trace_ref.append(value)
                if key[0]==0 and key[1]==0:fixed.extend(range(15*edge_lookup[key],15*(edge_lookup[key]+1)))
            else:
                if np.linalg.norm(trace_ref[edge_lookup[key]]-value)>1e-12:raise ValueError('reference_trace_mismatch')
            local.extend(range(15*edge_lookup[key],15*(edge_lookup[key]+1)))
        connect.append(local)
    trace_ref=np.vstack(trace_ref);free=np.setdiff1d(np.arange(len(trace_ref)),fixed)
    return dict(original=original,trace_ref=trace_ref,connect=np.array(connect),free=free,fixed=np.array(fixed),n=n,kappa=k,amplitude=amp)
