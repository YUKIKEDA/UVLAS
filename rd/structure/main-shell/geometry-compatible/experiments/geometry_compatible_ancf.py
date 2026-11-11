"""GC-ANCF: geometry-coordinate primitive/derivative closure, full 3D energy."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import svd
import c1_work_mean as base
old=base.old;m=base.m;f=base.f;KINDS=base.KINDS
def profiles(x,n,axis,kind):
    lo=0. if axis==0 else -.5;s=x-lo;knots=np.r_[np.zeros(6),np.arange(1,n)/n,np.ones(6)]
    B=BSpline(knots,np.eye(n+5),5,extrapolate=False);B1=B.derivative();B2=B.derivative(2);integ=B.antiderivative()
    if kind==0:
        a=B(s);return np.column_stack((a,s[:,None]*a-integ(s)+integ(0))),np.column_stack((B1(s),s[:,None]*B1(s)))
    return np.column_stack((B1(s),s[:,None]*B1(s))),np.column_stack((B2(s),B1(s)+s[:,None]*B2(s)))
def plane_rule(q,n,order):
    # Split bilinear physical-coordinate knot crossings before quadrature.
    g,w=np.polynomial.legendre.leggauss(order);corners=q.reshape(75,3)[[0,4,20,24],:2]
    yc=[-1.,1.]
    for axis in [0,1]:
        levels=np.arange(1,n)/n-(.5 if axis==1 else 0.)
        for a,b in [(0,2),(1,3)]:
            delta=corners[b,axis]-corners[a,axis]
            if abs(delta)>1e-14:
                yc.extend(2*(level-corners[a,axis])/delta-1 for level in levels if min(corners[a,axis],corners[b,axis])+1e-13<level<max(corners[a,axis],corners[b,axis])-1e-13)
    yc=np.unique(yc);points=[];weights=[]
    for y0,y1 in zip(yc[:-1],yc[1:]):
        for t,wt in zip((y0+y1)/2+(y1-y0)*g/2,w*(y1-y0)/2):
            edge=np.array([corners[0]*(1-t)/2+corners[2]*(1+t)/2,corners[1]*(1-t)/2+corners[3]*(1+t)/2]);xc=[-1.,1.]
            for axis in [0,1]:
                delta=edge[1,axis]-edge[0,axis]
                if abs(delta)>1e-14:
                    levels=np.arange(1,n)/n-(.5 if axis==1 else 0.)
                    xc.extend(2*(level-edge[0,axis])/delta-1 for level in levels if min(edge[:,axis])+1e-13<level<max(edge[:,axis])-1e-13)
            xc=np.unique(xc)
            for x0,x1 in zip(xc[:-1],xc[1:]):
                points.extend((x,t) for x in (x0+x1)/2+(x1-x0)*g/2);weights.extend(w*(x1-x0)*wt/2)
    return np.array(points),np.array(weights)
def pool(q,plane,n,kind):
    N,g=m.shape(np.column_stack((plane,np.zeros(len(plane)))),0);X=N@q.reshape(75,3)[:25];D=np.einsum('pna,nc->pca',g,q.reshape(75,3)[:25]);values=[];gx=[];gy=[]
    for axis in [0,1]:
        v,d=profiles(X[:,axis],n,axis,kind);values.append(v);gx.append(d*D[:,axis,0,None]);gy.append(d*D[:,axis,1,None])
    return np.column_stack(values),np.column_stack(gx),np.column_stack(gy)
def prepare_geometry(k,n,amp,order):
    geo=base.geometry_space(k,n,amp,order);geo['n']=n;completions=[]
    for kind in [0,1]:
        A=[];P=[]
        for cell,q in enumerate(geo['refs']):
            plane,w=plane_rule(q,n,order);R=base.surface(n,cell,plane,geo)[kind][0];v,_,_=pool(q,plane,n,kind)
            A.append(R*np.sqrt(w)[:,None]);P.append(v*np.sqrt(w)[:,None])
        A=np.vstack(A);P=np.vstack(P);coef=np.linalg.lstsq(A,P,rcond=1e-13)[0];res=P-A@coef;U,s,V=svd(res,full_matrices=False,check_finite=False)
        ranks=[int(np.count_nonzero(s>tol*max(1.,np.linalg.norm(P)))) for tol in [1e-9,1e-10,1e-11]]
        if len(set(ranks))!=1:raise ValueError('GC_enrichment_rank_ambiguous:'+str(ranks))
        rank=ranks[0];completions.append(dict(coeff=coef,V=V[:rank].T,s=s[:rank],rank=rank,ranks=ranks))
    geo['enrich']=completions;return geo
def fields(n,cell,q,plane,geo):
    result=base.surface(n,cell,plane,geo)
    for kind,comp in enumerate(geo['enrich']):
        v,gx,gy=pool(q,plane,n,kind);R,Rx,Ry=result[kind]
        result[kind]=tuple(np.column_stack((a,(b-a@comp['coeff'])@comp['V']/comp['s'])) for a,b in [(R,v),(Rx,gx),(Ry,gy)])
    return result
def kinematics(n,cell,q,h,p,geo):
    fs=fields(n,cell,q,p[:,:2],geo);z=h*p[:,2]/2;ns=sum(v.shape[1] for v,gx,gy in fs);N=np.zeros((len(p),ns));G=np.zeros((len(p),ns,3));off=0
    for power,(v,gx,gy) in enumerate(fs):
        size=v.shape[1];fac=np.ones(len(p)) if power==0 else z if power==1 else z*z/2;df=np.zeros(len(p)) if power==0 else np.full(len(p),h/2) if power==1 else z*h/2
        N[:,off:off+size]=fac[:,None]*v;G[:,off:off+size,0]=fac[:,None]*gx;G[:,off:off+size,1]=fac[:,None]*gy;G[:,off:off+size,2]=df[:,None]*v;off+=size
    X,J,det,inv,_=f.base.geometry(q,h,p);physical=np.einsum('pni,pij->pnj',G,inv,optimize=True);B=np.einsum('kij,pnj->pkni',m.MB,physical).reshape(len(p),6,3*ns)
    return N,physical,B,det,X
def volume_rule(q,n,order,zorder):
    plane,wp=plane_rule(q,n,order);z,wz=np.polynomial.legendre.leggauss(zorder)
    return np.column_stack((np.repeat(plane,zorder,axis=0),np.tile(z,len(plane)))),(wp[:,None]*wz).reshape(-1)
def build(k,n,amp,h,order=7,zorder=5,geo=None):
    geo=prepare_geometry(k,n,amp,order) if geo is None else geo
    fs=fields(n,0,geo['refs'][0],np.zeros((1,2)),geo);nr,nd,nc=[v.shape[1] for v,gx,gy in fs];ns=nr+nd+nc;cols=3*ns;rows=[];rhs=[];constraints=[];mass=np.zeros((ns,ns))
    for cell,q in enumerate(geo['refs']):
        p,w=volume_rule(q,n,order,zorder);N,G,B,det,X=kinematics(n,cell,q,h,p,geo);E,_=old.reference(q,p,h,k)
        if det.min()<=0:raise ValueError('GC_reference_orientation')
        weight=np.sqrt(w*det)[:,None,None];rows.append((np.einsum('ij,pjn->pin',m.CH,B)*weight).reshape(-1,cols)/h**1.5);rhs.append((np.einsum('ij,pjt->pit',m.CH,E)*weight).reshape(-1,3)/h**1.5);mass+=N.T@((w*det)[:,None]*N)
        if cell%n==0:
            edge=np.array([[-1.,t,c] for t in np.linspace(-1,1,13) for c in [-1.,0.,1.]])
            constraints.append(kinematics(n,cell,q,h,edge,geo)[0])
    A=np.vstack(rows);b=np.vstack(rhs);C=np.vstack(constraints);_,s,V=svd(C,full_matrices=True,check_finite=False);ranks=[int(np.count_nonzero(s>tol)) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('GC_boundary_rank_ambiguous:'+str(ranks))
    Zs=V[ranks[0]:].T;Z=np.kron(Zs,np.eye(3));ms=1/np.sqrt(np.diag(mass));mass_min=float(np.linalg.eigvalsh((mass*ms)*ms[:,None])[0])
    return dict(A=A,b=b,Z=Z,geo=geo,n=n,k=k,amp=amp,h=h,order=order,zorder=zorder,nr=nr,nd=nd,nc=nc,ns=ns,mass_scaled_min=mass_min,boundary_ranks=ranks,reference_norm=np.sum(b*b,axis=0))
solve=base.solve
def position_error(a,coef,order=15,zorder=9):
    error=np.zeros(3);norm=np.zeros(3)
    for cell,q in enumerate(a['geo']['refs']):
        p,w=volume_rule(q,a['n'],order,zorder);N,G,B,det,X=kinematics(a['n'],cell,q,a['h'],p,a['geo']);_,target=old.reference(q,p,a['h'],a['k']);u=(N@coef.reshape(a['ns'],9)).reshape(len(p),3,3)
        error+=np.einsum('p,pat,pat->t',w*det,u-target,u-target);norm+=np.einsum('p,pat,pat->t',w*det,target,target)
    return np.sqrt(error/norm).tolist()
