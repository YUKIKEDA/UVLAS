"""Load-driven PG linear solver; manufactured displacements are observation only."""
import numpy as np
from scipy.linalg import svd,solve_triangular
import paired_geometry_ancf as paired
import pg_continuum_loads as loads
import geometry_compatible_ancf_validation_stream as validation
f=paired.f;m=f.m

def xy(q,p):
    N,G=m.shape(np.column_stack((p,np.zeros(len(p)))),0);Q=q.reshape(75,3)[:25];X=N@Q;J=np.einsum('pni,na->pai',G,Q)
    return X,J,J[:,0,0]*J[:,1,1]-J[:,0,1]*J[:,1,0]

def components(a,cell,q,p):
    fs=paired.base.surface(a['n'],cell,p,a['geo']);R,D,C=[v for v,_,_ in fs];p3=np.column_stack((p,np.zeros(len(p))));er,eg=paired.extra(a['n'],cell,q,2.,p3,a['geo']);ed=eg[:,:,2];out=np.zeros((3,len(p),a['ns']));nr,nd,nc=a['nr'],a['nd'],a['nc'];out[0,:,:nr]=R;out[1,:,nr:nr+nd]=D;out[2,:,nr+nd:nr+nd+nc]=C;out[0,:,nr+nd+nc:]=er;out[1,:,nr+nd+nc:]=ed
    return out

def prepare(k,n,amp,h,order=13):
    geo=f.c.prepare_geometry(k,n,amp,order);fs=paired.base.surface(n,0,np.zeros((1,2)),geo);nr,nd,nc=[v.shape[1] for v,_,_ in fs];ns=nr+nd+nc+geo['enrich'][0]['rank'];a=dict(k=k,n=n,amp=amp,h=h,geo=geo,nr=nr,nd=nd,nc=nc,ns=ns);constraints=[];rows=[]
    for cell,q in enumerate(geo['refs']):
        p,w=f.c.volume_rule(q,n,order,7);R=None
        for pp,ww in f.chunks(p,w):
            _,_,B,det,_=f.c.kinematics(n,cell,q,h,pp,geo);A=(np.einsum('ij,pjn->pin',m.CH,B)*np.sqrt(ww*det)[:,None,None]/h**1.5).reshape(-1,3*ns);R=f.old.qrtools.rfactor(np.asfortranarray(A if R is None else np.vstack((R,A))))
        rows.append(R)
        if cell%n==0:
            edge=np.array([[-1.,t,z] for t in np.linspace(-1,1,33) for z in [-1.,0.,1.]]);constraints.append(f.shape_values(a,cell,q,edge))
    C=np.vstack(constraints);_,s,V=svd(C,full_matrices=True,check_finite=False);ranks=[int(sum(s>tol)) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('support_rank:'+str(ranks))
    Zs=V[ranks[0]:].T;Z=np.kron(Zs,np.eye(3));A=np.vstack(rows);AZ=A@Z;scale=1/np.linalg.norm(AZ,axis=0);R=f.old.qrtools.rfactor(np.asfortranarray(AZ*scale));sv=svd(R,compute_uv=False,check_finite=False);fr=[int(sum(sv>tol)) for tol in [1e-9,1e-10,1e-11]]
    if fr!=[Z.shape[1]]*3:raise ValueError('stiffness_rank:'+str(fr))
    a.update(A=A,Zs=Zs,Z=Z,scale=scale,R=R[:Z.shape[1]],rcond=f.old.qrtools.rcond(R[:Z.shape[1]]),support_ranks=ranks);return a

def external(a,mode,order=15,count=14,beta=None):
    """Integrate only body/traction moments against scalar position functions."""
    F=np.zeros((a['ns'],3,3));fixed=np.zeros_like(F);coarse=np.zeros((8,3));k,h,n=a['k'],a['h'],a['n']
    def point_data(q,p):
        X,J,jac=xy(q,p);x=X[:,0];v,vp=loads.curvature(x,mode);K,Kx,U=loads.global_basis(x,k)
        if beta is not None:v=v-K@beta;vp=vp-Kx@beta
        W,Fx,Fy=loads.moments(x,k,h,v,vp,count);return J,jac,U,W,Fx,Fy
    for cell,q in enumerate(a['geo']['refs']):
        p,w=f.c.plane_rule(q,n,order);J,jac,U,W,_,_=point_data(q,p);S=components(a,cell,q,p);weight=w*jac
        F+=np.einsum('p,lpn,lpat->nat',weight,S,W,optimize=True);coarse+=np.einsum('p,lpan,lpat->nt',weight,U,W,optimize=True)
        row,col=divmod(cell,n)
        for axis,sign,boundary in [(0,-1,col==0),(0,1,col==n-1),(1,-1,row==0),(1,1,row==n-1)]:
            if not boundary:continue
            # Split physical curvature knots on boundary edges, using the established rule.
            p3,ww=validation.face_rule(q,n,axis,sign,order,1);p=p3[:,:2];J,jac,U,_,Fx,Fy=point_data(q,p);S=components(a,cell,q,p);traction=sign*(Fx if axis==0 else Fy);weight=ww/2*J[:,1-axis,1-axis]
            term=np.einsum('p,lpn,lpat->nat',weight,S,traction,optimize=True)
            if axis==0 and sign<0:fixed+=term
            else:F+=term;coarse+=np.einsum('p,lpan,lpat->nt',weight,U,traction,optimize=True)
    return F.reshape(3*a['ns'],3)/h**3,fixed.reshape(3*a['ns'],3)/h**3,coarse/h**3

def coarse_factor(k,h):
    x,wx=np.polynomial.legendre.leggauss(40);x=(x+1)/2;wx=wx/2;t,wt=np.polynomial.legendre.leggauss(9);t=t/2;wt=wt/2;g=1+k*k*x*x;weight=m.C[0,0]*g*g*np.sum(wt[None,:]*t[None,:]**2/(g[:,None]-k*h*t[None,:]),axis=1);K,_,_=loads.global_basis(x,k);return f.old.qrtools.rfactor(np.asfortranarray(np.sqrt(wx*weight)[:,None]*K))

def solve(a,mode):
    F,fixed,L=external(a,mode);Rb=coarse_factor(a['k'],a['h']);beta=solve_triangular(Rb,solve_triangular(Rb.T,L,lower=True,check_finite=False),check_finite=False);g,_,_=external(a,mode,beta=beta)
    rhs=a['scale'][:,None]*(a['Z'].T@g);z=solve_triangular(a['R'].T,rhs,lower=True,check_finite=False);y=solve_triangular(a['R'],z,check_finite=False);coef=a['Z']@(a['scale'][:,None]*y)
    return dict(coef=coef,beta=beta,F=F,fixed=fixed,residual_force=g,linear=float(max(np.linalg.norm(a['R'].T@z-rhs,axis=0)/np.maximum(np.linalg.norm(rhs,axis=0),1e-30))))

def observe(a,s,mode,order=15):
    import geometry_compatible_ancf_validation_stream as reference
    reference.select(mode);error=np.zeros(3);norm=np.zeros(3);pe=np.zeros(3);pn=np.zeros(3);internal=np.zeros_like(s['F']);energy=np.zeros(3);position_work=np.zeros(3)
    for cell,q in enumerate(a['geo']['refs']):
        p,w=f.c.volume_rule(q,a['n'],order,9)
        for pp,ww in f.chunks(p,w):
            N,G,B,det,_=f.c.kinematics(a['n'],cell,q,a['h'],pp,a['geo']);Eg,Ug=loads.global_fields(q,pp,a['h'],a['k'],m);E,U=reference.f.old.reference(q,pp,a['h'],a['k']);strain=B@s['coef']+Eg@s['beta'];u=(N@s['coef'].reshape(a['ns'],9)).reshape(len(pp),3,3)+Ug@s['beta'];delta=strain-E;stress=m.C@strain;weight=ww*det/a['h']**3
            error+=np.einsum('p,pkt,kl,plt->t',weight,delta,m.C,delta);norm+=np.einsum('p,pkt,kl,plt->t',weight,E,m.C,E);pe+=np.einsum('p,pat,pat->t',ww*det,u-U,u-U);pn+=np.einsum('p,pat,pat->t',ww*det,U,U);internal+=np.einsum('p,pkn,pkt->nt',weight,B,stress,optimize=True);energy+=np.einsum('p,pkt,pkt->t',weight,strain,stress)
    residual=internal-s['F'];g=a['scale'][:,None]*(a['Z'].T@residual);stationarity=np.linalg.norm(g,axis=0)/np.sqrt(a['Z'].shape[1]*norm)
    return dict(error=np.sqrt(error/norm).tolist(),position=np.sqrt(pe/pn).tolist(),stationarity=stationarity.tolist(),energy=energy.tolist()),residual
