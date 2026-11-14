"""C1-WM: conforming absolute-position kinematics and work-completed mean strains."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import svd,qr,solve_triangular
import mixed_mean_stress_finite as old
m=old.m;f=old.f;KINDS=old.KINDS
def rule(order,zorder):
    p,w=f.base.volume_points(order,zorder);plane,wp=m.quadrature(order,2);z,wz=np.polynomial.legendre.leggauss(zorder)
    return p,w,plane,wp,z,wz
def spline(t,n):
    knots=np.r_[np.zeros(3),np.arange(1,n)/n,np.ones(3)];B=BSpline(knots,np.eye(n+2),2,extrapolate=False)
    return B(t),B.derivative()(t)
def lagrange(x,p):
    nodes=np.linspace(-1,1,p+1);v=[];d=[]
    for i in range(p+1):
        poly=np.poly1d([1.])
        for j in range(p+1):
            if j!=i:poly*=np.poly1d([1.,-nodes[j]])/(nodes[i]-nodes[j])
        v.append(poly(x));d.append(np.polyder(poly)(x))
    return np.array(v).T,np.array(d).T
def surface(n,cell,plane,geo=None):
    j,i=divmod(cell,n);s=(i+(plane[:,0]+1)/2)/n;t=(j+(plane[:,1]+1)/2)/n
    bs,ds=spline(s,n);bt,dt=spline(t,n)
    R=np.einsum('pj,pi->pji',bt,bs).reshape(len(plane),-1)
    Rx=np.einsum('pj,pi->pji',bt,ds).reshape(len(plane),-1)/(2*n);Ry=np.einsum('pj,pi->pji',dt,bs).reshape(len(plane),-1)/(2*n)
    if geo is not None and geo['rank']:
        p3=np.column_stack((plane,np.zeros(len(plane))));N,g=m.shape(p3,0);q=geo['refs'][cell].reshape(75,3)[:25]
        X=N@q;Xx=g[:,:,0]@q;Xy=g[:,:,1]@q;V=geo['V'];sc=geo['s']
        R,Rx,Ry=[np.column_stack((a,(b-a@geo['coeff'])@V/sc)) for a,b in [(R,X),(Rx,Xx),(Ry,Xy)]]
    result=[(R,Rx,Ry)]
    for degree in [2,1]:
        a,da=lagrange(plane[:,0],degree);b,db=lagrange(plane[:,1],degree);count=(degree*n+1)**2
        val=np.zeros((len(plane),count));gx=np.zeros_like(val);gy=np.zeros_like(val)
        for y in range(degree+1):
            for x in range(degree+1):
                node=(degree*j+y)*(degree*n+1)+degree*i+x
                val[:,node]=a[:,x]*b[:,y];gx[:,node]=da[:,x]*b[:,y];gy[:,node]=a[:,x]*db[:,y]
        result.append((val,gx,gy))
    return result
def geometry_space(k,n,amp,order):
    grid=m.mesh(k,n,amp);refs=[grid['ref'][ids].reshape(-1) for ids in grid['cells']];plane,wp=m.quadrature(order,2);A=[];X=[]
    for cell,q in enumerate(refs):
        R=surface(n,cell,plane)[0][0];N,_=m.shape(np.column_stack((plane,np.zeros(len(plane)))),0)
        A.append(R*np.sqrt(wp)[:,None]);X.append((N@q.reshape(75,3)[:25])*np.sqrt(wp)[:,None])
    A=np.vstack(A);X=np.vstack(X);coef=np.linalg.lstsq(A,X,rcond=1e-13)[0];res=X-A@coef;U,s,V=svd(res,full_matrices=False)
    rank=int(np.count_nonzero(s>1e-11*max(1.,np.linalg.norm(X))))
    return dict(grid=grid,refs=refs,coeff=coef,V=V[:rank].T,s=s[:rank],rank=rank)
def kinematics(n,cell,q,h,p,geo):
    fields=surface(n,cell,p[:,:2],geo);z=h*p[:,2]/2;ns=sum(v.shape[1] for v,gx,gy in fields)
    N=np.zeros((len(p),ns));G=np.zeros((len(p),ns,3));off=0
    for power,(v,gx,gy) in enumerate(fields):
        count=v.shape[1];fac=np.ones(len(p)) if power==0 else z if power==1 else z*z/2
        derivative=np.zeros(len(p)) if power==0 else np.full(len(p),h/2) if power==1 else z*h/2
        N[:,off:off+count]=fac[:,None]*v;G[:,off:off+count,0]=fac[:,None]*gx;G[:,off:off+count,1]=fac[:,None]*gy;G[:,off:off+count,2]=derivative[:,None]*v;off+=count
    X,J,det,inv,_=f.base.geometry(q,h,p)
    physical=np.einsum('pni,pij->pnj',G,inv,optimize=True)
    B=np.einsum('kij,pnj->pkni',m.MB,physical).reshape(len(p),6,3*ns)
    return N,physical,B,det,X
def build(k,n,amp,h,order=7,zorder=5):
    geo=geometry_space(k,n,amp,order);p,w,plane,wp,z,wz=rule(order,zorder);nr=(n+2)**2+geo['rank'];nd=(2*n+1)**2;nc=(n+1)**2;ns=nr+nd+nc;cols=3*ns
    nodeB=np.zeros(((n+1)**2,6,cols));nodeStar=np.zeros(((n+1)**2,6,3));hits=np.zeros((n+1)**2);cell_data=[];corners=np.array([[-1.,-1.],[1.,-1.],[-1.,1.],[1.,1.]])
    for cell,q in enumerate(geo['refs']):
        j,i=divmod(cell,n);pp=np.array([[a,b,c] for a,b in corners for c in z]);_,_,B,_,_=kinematics(n,cell,q,h,pp,geo)
        E,_=old.reference(q,pp,h,k);barB=np.einsum('z,pzkn->pkn',wz/2,B.reshape(4,zorder,6,cols));barE=np.einsum('z,pzkt->pkt',wz/2,E.reshape(4,zorder,6,3))
        ids=np.array([j*(n+1)+i,j*(n+1)+i+1,(j+1)*(n+1)+i,(j+1)*(n+1)+i+1]);nodeB[ids]+=barB;nodeStar[ids]+=barE;hits[ids]+=1;cell_data.append(ids)
    nodeB/=hits[:,None,None];nodeStar/=hits[:,None,None]
    H=[];hb=[];raw=[];ab=[];target=[];at=[];mean_weights=[];mass=np.zeros((ns,ns));geometry_error=0.;constraints=[]
    for cell,q in enumerate(geo['refs']):
        N,G,B,det,X=kinematics(n,cell,q,h,p,geo);E,u=old.reference(q,p,h,k)
        if det.min()<=0:raise ValueError('C1WM_reference_orientation')
        gamma=(1/det.reshape(len(plane),zorder))@(wz/2);factor=(1/(gamma[:,None]*det.reshape(len(plane),zorder))).reshape(-1)
        bbar=np.einsum('z,pzkn->pkn',wz/2,B.reshape(len(plane),zorder,6,cols));estar=np.einsum('z,pzkt->pkt',wz/2,E.reshape(len(plane),zorder,6,3))
        high=B-factor[:,None,None]*np.repeat(bbar,zorder,axis=0);he=E-factor[:,None,None]*np.repeat(estar,zorder,axis=0)
        H.append((np.einsum('ij,pjn->pin',m.CH,high)*np.sqrt(w*det)[:,None,None]).reshape(-1,cols)/h**1.5);hb.append((np.einsum('ij,pjt->pit',m.CH,he)*np.sqrt(w*det)[:,None,None]).reshape(-1,3)/h**1.5)
        weights=2*wp/gamma;mean_weights.append(weights);a,_=lagrange(plane[:,0],1);b,_=lagrange(plane[:,1],1);interp=np.einsum('pj,pi->pji',b,a).reshape(len(plane),4)
        assumed=np.einsum('pc,ckn->pkn',interp,nodeB[cell_data[cell]]);ae=np.einsum('pc,ckt->pkt',interp,nodeStar[cell_data[cell]])
        weight=np.sqrt(weights)[:,None,None]
        raw.append((np.einsum('ij,pjn->pin',m.CH,bbar)*weight).reshape(-1,cols));ab.append((np.einsum('ij,pjn->pin',m.CH,assumed)*weight).reshape(-1,cols))
        target.append((np.einsum('ij,pjt->pit',m.CH,estar)*weight).reshape(-1,3));at.append((np.einsum('ij,pjt->pit',m.CH,ae)*weight).reshape(-1,3))
        mass+=N.T@((w*det)[:,None]*N)
        if cell%n==0:
            edge=np.array([[-1.,t,c] for t in np.linspace(-1,1,7) for c in [-1.,0.,1.]])
            constraints.append(kinematics(n,cell,q,h,edge,geo)[0])
    raw=np.vstack(raw);ab=np.vstack(ab);target=np.vstack(target);at=np.vstack(at);weights=np.concatenate(mean_weights)
    constants=(np.sqrt(weights)[:,None,None]*np.broadcast_to(m.CH,(len(weights),6,6))).reshape(-1,6)
    W=np.column_stack((raw[:,3*nr:3*(nr+nd)],constants));scale=1/np.maximum(np.linalg.norm(W,axis=0),1e-300);Q,s,V=svd(W*scale,full_matrices=False,check_finite=False)
    ranks=[int(np.count_nonzero(s>tol)) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('C1WM_work_rank_ambiguous:'+str(ranks))
    Q=Q[:,:ranks[0]];mean=ab+Q@(Q.T@(raw-ab));bmean=at+Q@(Q.T@(target-at))
    work_defect=float(np.linalg.norm(constants.T@(mean-raw))/max(1.,np.linalg.norm(constants)*np.linalg.norm(raw)))
    A=np.vstack(H+[mean/h**1.5]);bb=np.vstack(hb+[bmean/h**1.5]);C=np.vstack(constraints);_,sc,Vc=svd(C,full_matrices=True,check_finite=False);cr=[int(np.count_nonzero(sc>tol)) for tol in [1e-9,1e-10,1e-11]]
    if len(set(cr))!=1:raise ValueError('C1WM_boundary_rank_ambiguous:'+str(cr))
    Zs=Vc[cr[0]:].T;Z=np.kron(Zs,np.eye(3));Mscale=1/np.sqrt(np.diag(mass));mass_min=float(np.linalg.eigvalsh((mass*Mscale)*Mscale[:,None])[0])
    return dict(A=A,b=bb,Z=Z,geo=geo,n=n,k=k,amp=amp,h=h,order=order,zorder=zorder,nr=nr,nd=nd,nc=nc,ns=ns,work_ranks=ranks,boundary_ranks=cr,work_defect=work_defect,mass_scaled_min=mass_min,reference_norm=np.sum(bb*bb,axis=0))
def solve(a,free_rank=True):
    A=a['A'];b=a['b'];Z=a['Z'];AZ=A@Z;scale=1/np.linalg.norm(AZ,axis=0);R=old.qrtools.rfactor(np.asfortranarray(np.column_stack((AZ*scale,b))));n=AZ.shape[1];sv=svd(R[:,:n],compute_uv=False,check_finite=False);ranks=[int(np.count_nonzero(sv>tol)) for tol in [1e-9,1e-10,1e-11]]
    if ranks!=[n]*3:return dict(passed_rank=False,ranks=ranks,columns=n,smallest=sv[-8:].tolist())
    y=solve_triangular(R[:n,:n],R[:n,n:],check_finite=False);coef=Z@(scale[:,None]*y);res=A@coef-b
    kkt=np.linalg.norm(scale[:,None]*(AZ.T@res),axis=0)/np.maximum(1e-30,np.sqrt(n)*np.linalg.norm(b,axis=0));linear=np.linalg.norm(R[:n,:n]@y-R[:n,n:],axis=0)/np.linalg.norm(R[:n,n:],axis=0)
    out=dict(passed_rank=True,coef=coef,error=np.sqrt(np.sum(res*res,axis=0)/a['reference_norm']).tolist(),kkt=kkt.tolist(),linear=linear.tolist(),rcond=old.qrtools.rcond(R[:n,:n]),ranks=ranks,columns=n)
    if free_rank:
        scale0=1/np.linalg.norm(A,axis=0);T=old.qrtools.rfactor(np.asfortranarray(A*scale0));s=svd(T,compute_uv=False,check_finite=False);out['free_nullities']=[int(np.count_nonzero(s<=tol)) for tol in [1e-9,1e-10,1e-11]]
    return out
def position_error(a,coef,order=15,zorder=9):
    p,w,_,_,_,_=rule(order,zorder);err=np.zeros(3);norm=np.zeros(3)
    for cell,q in enumerate(a['geo']['refs']):
        N,_,_,det,_=kinematics(a['n'],cell,q,a['h'],p,a['geo']);_,target=old.reference(q,p,a['h'],a['k']);u=(N@coef.reshape(a['ns'],9)).reshape(len(p),3,3)
        err+=np.einsum('p,pat,pat->t',w*det,u-target,u-target);norm+=np.einsum('p,pat,pat->t',w*det,target,target)
    return np.sqrt(err/norm).tolist()
