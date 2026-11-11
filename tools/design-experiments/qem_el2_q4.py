"""Research-only EL2 Q4 Lobatto realization and same-space Gauss control."""
import itertools
import numpy as np
from scipy.linalg import svd
from scipy.optimize import brentq
import q3_curved_element as old
MB=old.MB;C=old.C;CH=old.CH
DEGREES=[4,4,4];SLICES=[slice(0,25),slice(25,50),slice(50,75)]
LOBATTO=np.array([-1.,-np.sqrt(3/7),0.,np.sqrt(3/7),1.])
LOBATTO_W=2/(20*np.polynomial.legendre.legval(LOBATTO,[0,0,0,0,1])**2)
NODES=[LOBATTO.copy() for _ in range(3)]

def shape(points,field,independent=False):
    nodes=NODES[field];x,y=points[:,:2].T
    def cardinal(t):
        if independent:
            coeff=np.linalg.inv(np.vander(nodes,len(nodes),increasing=True))
            powers=np.arange(len(nodes));values=t[:,None]**powers@coeff
            deriv=np.column_stack((np.zeros(len(t)),*[i*t**(i-1) for i in range(1,len(nodes))]))@coeff
            return values,deriv
        values=[];derivatives=[]
        for i,a in enumerate(nodes):
            others=[b for j,b in enumerate(nodes) if j!=i];den=np.prod([a-b for b in others])
            values.append(np.prod([t-b for b in others],axis=0)/den)
            derivatives.append(sum(np.prod([t-c for k,c in enumerate(others) if k!=j],axis=0) for j in range(len(others)))/den)
        return np.array(values).T,np.array(derivatives).T
    ls,ds=cardinal(x);lt,dt=cardinal(y)
    return np.einsum('pi,pj->pij',lt,ls).reshape(len(points),-1),np.stack((
        np.einsum('pi,pj->pij',lt,ds).reshape(len(points),-1),np.einsum('pi,pj->pij',dt,ls).reshape(len(points),-1)),axis=2)

def tables(points,z,independent=False):
    z=np.broadcast_to(z,(len(points),));phi=np.zeros((len(points),75),dtype=np.result_type(z,float))
    g=np.zeros((len(points),75,3),dtype=phi.dtype);g1=np.zeros_like(g)
    for f,sl in enumerate(SLICES):
        N,d=shape(points,f,independent);factor=[1.,1.,.5][f]
        phi[:,sl]=N*factor*z[:,None]**f;g[:,sl,:2]=d*factor*z[:,None,None]**f
        if f:g[:,sl,2]=N*factor*f*z[:,None]**(f-1);g1[:,sl,:2]=d*factor*f*z[:,None,None]**(f-1)
        if f==2:g1[:,sl,2]=N
    return phi,g,g1

def raw(q,points,z,independent=False):
    phi,g,_=tables(points,z,independent);J=np.einsum('na,pnb->pab',q.reshape(75,3),g)
    G=g@np.linalg.inv(J);B=np.einsum('kij,pnj->pkni',MB,G).reshape(len(points),6,225)
    return phi,J,B

def expansion(q,points):
    _,g,g1=tables(points,0.);coords=q.reshape(75,3);J0=np.einsum('na,pnb->pab',coords,g);J1=np.einsum('na,pnb->pab',coords,g1)
    Ji=np.linalg.inv(J0);G0=g@Ji;G1=g1@Ji-G0@J1@Ji
    return J0,np.einsum('kij,pnj->pkni',MB,G0).reshape(-1,6,225),np.einsum('kij,pnj->pkni',MB,G1).reshape(-1,6,225)

def quadrature(order,dimension=3):
    p,w=np.polynomial.legendre.leggauss(order)
    return np.array(list(itertools.product(p,repeat=dimension))),np.array([np.prod(t) for t in itertools.product(w,repeat=dimension)])

def mesh(k,n,amp):
    vertices,_,_=old.old_grid.mesh(n,amp);vertices=vertices.copy();vertices[:,:,1]-=.5
    xy=[];connect=[];offset=0;field_ids=[]
    for field,degree in enumerate(DEGREES):
        side=degree*n+1;coords=np.full((side*side,2),np.nan);ids_all=[]
        for j,i in itertools.product(range(n),range(n)):
            ids=np.array([(degree*j+b)*side+degree*i+a for b,a in itertools.product(range(degree+1),repeat=2)])
            points=np.array([[s,t] for t,s in itertools.product(NODES[field],repeat=2)])
            corners=vertices[[j,j,j+1,j+1],[i,i+1,i+1,i]];values,_=old.bilinear(points,corners)
            existing=np.all(np.isfinite(coords[ids]),axis=1)
            if np.any(existing):assert np.max(abs(coords[ids][existing]-values[existing]))<1e-13
            coords[ids]=values;ids_all.append(ids+offset)
        xy.append(coords);field_ids.append(np.arange(offset,offset+len(coords)));offset+=len(coords);connect.append(ids_all)
    ref=np.zeros((offset,3));ref[field_ids[0]]=np.column_stack((xy[0],k*xy[0][:,0]**2/2))
    ref[field_ids[1]]=np.column_stack((-k*xy[1][:,0],np.zeros(len(xy[1])),np.ones(len(xy[1]))))
    cells=np.array([np.concatenate([connect[f][cell] for f in range(3)]) for cell in range(n*n)])
    free=[ids[np.arange(len(ids))%(degree*n+1)!=0] for ids,degree in zip(field_ids,DEGREES)]
    return dict(xy=xy,ref=ref,cells=cells,field_ids=field_ids,free=free)

def displacement(mesh,kind,k):
    u=np.zeros_like(mesh['ref']);xr=mesh['xy'][0][:,0];xd=mesh['xy'][1][:,0]
    if kind=='quadratic':U=-k*xr**3/3;W=xr*xr/2;dp=xd
    elif kind=='cubic':U=-k*xr**4/8;W=xr**3/6;dp=xd*xd/2
    elif kind=='exponential':U=-k*((xr-1)*np.exp(xr)+1-xr*xr/2);W=np.exp(xr)-1-xr;dp=np.exp(xd)-1
    else:U=-k*xr**5/15;W=xr**4/12;dp=xd**3/3
    u[mesh['field_ids'][0]]=np.column_stack((U,np.zeros(len(U)),W))
    u[mesh['field_ids'][1]]=np.column_stack((-dp,np.zeros(len(dp)),-k*xd*dp))
    return u

def target(kind,k,x):
    v={'quadratic':np.ones_like(x),'cubic':x,'quartic':x*x,'exponential':np.exp(x)}[kind];a=np.column_stack((np.ones_like(x),np.zeros_like(x),k*x))
    return np.einsum('kij,pij->pk',MB,-v[:,None,None]*a[:,:,None]*a[:,None,:]/(1+k*k*x*x)[:,None,None])

def reference_h1(kind,k,order):
    p,w=np.polynomial.legendre.leggauss(order);x=(p+1)/2;w=w/2;D=1+k*k*x*x
    dp,dd={'quadratic':(x,np.ones_like(x)),'cubic':(x*x/2,x),'quartic':(x**3/3,x*x),'exponential':(np.exp(x)-1,np.exp(x))}[kind]
    value=dp*dp*D+(dd*dd+k*k*(dp+x*dd)**2)/D
    return float(np.sqrt((w@(value*np.sqrt(D)))/(w@np.sqrt(D))))

def boundary(q,h,faces=None,independent=True):
    if faces is None:faces=list(itertools.product(range(3),[-1.,1.]))
    plane,w=quadrature(9,2);force=np.zeros((75,3,6))
    for axis,sign in faces:
        points=np.zeros((len(plane),3));points[:,axis]=sign;points[:,[i for i in range(3) if i!=axis]]=plane
        phi,J,_=raw(q,points,h*points[:,2]/2,independent)
        J=J.copy();J[:,:,2]*=h/2
        normal=sign*np.cross(J[:,:,(axis+1)%3],J[:,:,(axis+2)%3]);traction=np.einsum('kab,pb->pak',MB,normal)
        force+=np.einsum('p,pn,pak->nak',w,phi,traction)
    return force.reshape(225,6)

def native_plane():
    return np.array([[s,t] for t,s in itertools.product(LOBATTO,repeat=2)]),np.array([a*b for a,b in itertools.product(LOBATTO_W,repeat=2)])

def sampled_expansion(q,points,conservative=False):
    native,wn=native_plane();J,B0,B1=expansion(q,native);N,_=shape(points,0)
    R0=np.einsum('pn,nki->pki',N,B0);R1=np.einsum('pn,nki->pki',N,B1)
    if conservative:
        pg,wg=quadrature(9,2);Jg,G0,G1=expansion(q,pg);Ng,_=shape(pg,0)
        mu=wg*np.linalg.det(Jg);mu/=mu.sum()
        delta=np.einsum('p,pki->ki',mu,G0-np.einsum('pn,nki->pki',Ng,B0))
        R0+=delta
    Jp,_,_=expansion(q,points)
    return Jp,R0,R1

def native_volume(order_z=5):
    plane,wp=native_plane();z,wz=np.polynomial.legendre.leggauss(order_z)
    points=np.array([[s,t,c] for (s,t),c in itertools.product(plane,z)])
    weights=np.array([a*b for a,b in itertools.product(wp,wz)])
    return points,weights

def project_fields(q,h,Eraw,Enative,order=9,order_z=5,conservative=True):
    plane,wp=quadrature(order,2);z,wz=np.polynomial.legendre.leggauss(order_z)
    points=np.array([[s,t,c] for (s,t),c in itertools.product(plane,z)]);w=np.array([a*b for a,b in itertools.product(wp,wz)])
    _,J,_=raw(q,points,h*points[:,2]/2);mu=w*np.linalg.det(J);mu/=mu.sum()
    N,_=shape(plane,0);values=Enative.reshape(25,order_z,6,-1)
    interpolated=np.einsum('pn,nzki->pzki',N,values).reshape(len(points),6,-1)
    correction=np.einsum('p,pki->ki',mu,Eraw-interpolated) if conservative else np.zeros_like(interpolated[0])
    return points,mu,interpolated+correction

def local_system(q,h,conservative=True):
    plane,wp=quadrature(9,2);z,wz=np.polynomial.legendre.leggauss(5)
    p=np.array([[s,t,c] for (s,t),c in itertools.product(plane,z)]);w=np.array([a*b for a,b in itertools.product(wp,wz)])
    pn,wn=native_volume();_,J,B=raw(q,p,h*p[:,2]/2);_,Jn,Bn=raw(q,pn,h*pn[:,2]/2)
    _,mu,R=project_fields(q,h,B,Bn,conservative=conservative)
    return p,w,J,B,R,pn,wn,Jn,Bn

def nonlinear(q,qcur,h,v1,v2,conservative=True):
    p,w,J,B,R,pn,wn,Jn,Bn=local_system(q,h,conservative)
    def data(points):
        _,g,_=tables(points,h*points[:,2]/2);J0=np.einsum('na,pnb->pab',q.reshape(75,3),g);G=g@np.linalg.inv(J0)
        F=np.eye(3)+np.einsum('na,pnb->pab',(qcur-q).reshape(75,3),G)
        E=np.einsum('kij,pij->pk',MB,(F.transpose(0,2,1)@F-np.eye(3))/2)
        Bcur=np.einsum('kij,pai,pnj->pkna',MB,F,G).reshape(len(points),6,225)
        d1=np.einsum('na,pnb->pab',v1.reshape(75,3),G);d2=np.einsum('na,pnb->pab',v2.reshape(75,3),G)
        second=np.einsum('kij,pij->pk',MB,(d1.transpose(0,2,1)@d2+d2.transpose(0,2,1)@d1)/2)
        return E,Bcur,second
    e,b,d2=data(p);en,bn,d2n=data(pn)
    _,mu,er=project_fields(q,h,e[:,:,None],en[:,:,None],conservative=conservative);er=er[:,:,0]
    _,_,br=project_fields(q,h,b,bn,conservative=conservative)
    _,_,d2r=project_fields(q,h,d2[:,:,None],d2n[:,:,None],conservative=conservative)
    stress=er@C.T;energy=np.einsum('p,pk,pk->',mu,er,stress)/2;force=np.einsum('p,pki,pk->i',mu,br,stress)
    material=np.einsum('p,pk,kl,pl->',mu,br@v1,C,br@v2)
    geometric=np.einsum('p,pk,pk->',mu,stress,d2r[:,:,0])
    return energy,force,material+geometric

def assemble(mesh,k,order,method='raw'):
    p,w=quadrature(order,2);free_rd=np.r_[mesh['free'][0],mesh['free'][1]];free_c=mesh['free'][2]
    rd_lookup={int(v):i for i,v in enumerate(free_rd)};c_lookup={int(v):i for i,v in enumerate(free_c)}
    rows=len(mesh['cells'])*len(p)*6;A0=np.zeros((rows,3*len(free_rd)));A1=np.zeros_like(A0);D1=np.zeros((rows,3*len(free_c)))
    targets={kind:np.zeros(rows) for kind in ['quadratic','cubic','quartic','exponential']}
    dc=np.r_[mesh['free'][1],mesh['free'][2]];dc_lookup={int(v):i for i,v in enumerate(dc)};M=np.zeros((len(dc),len(dc)));area=0.;zero_c=0.
    for cell,ids in enumerate(mesh['cells']):
        q=mesh['ref'][ids].reshape(225)
        J,B0,B1=expansion(q,p) if method=='raw' else sampled_expansion(q,p,method=='conservative')
        weight=np.sqrt(w*np.linalg.det(J));sl=slice(cell*len(p)*6,(cell+1)*len(p)*6)
        W0=np.einsum('ij,pjn->pin',CH,B0)*weight[:,None,None];W1=np.einsum('ij,pjn->pin',CH,B1)*weight[:,None,None]
        zero_c=max(zero_c,float(np.linalg.norm(B0[:,:,150:])))
        for a,node in enumerate(ids):
            if int(node) in rd_lookup:
                j=rd_lookup[int(node)];A0[sl,3*j:3*j+3]=W0[:,:,3*a:3*a+3].reshape(-1,3);A1[sl,3*j:3*j+3]=W1[:,:,3*a:3*a+3].reshape(-1,3)
            elif int(node) in c_lookup:
                j=c_lookup[int(node)];D1[sl,3*j:3*j+3]=W1[:,:,3*a:3*a+3].reshape(-1,3)
        Nr,gr=shape(p,0);Jsurf=np.einsum('na,pnb->pab',mesh['ref'][ids[:25]],gr);G=Jsurf.transpose(0,2,1)@Jsurf;aw=w*np.sqrt(np.linalg.det(G));area+=aw.sum()
        for field in [1,2]:
            N,g=shape(p,field);local=np.einsum('p,pi,pj->ij',aw,N,N)+np.einsum('p,pia,pab,pjb->ij',aw,g,np.linalg.inv(G),g)
            nodes=ids[SLICES[field]];keep=[a for a,node in enumerate(nodes) if int(node) in dc_lookup];where=[dc_lookup[int(nodes[a])] for a in keep]
            M[np.ix_(where,where)]+=local[np.ix_(keep,keep)]
        x=(Nr@mesh['ref'][ids[:25]])[:,0]
        for kind in targets:targets[kind][sl]=(target(kind,k,x)@CH.T*weight[:,None]).reshape(-1)
    return dict(A0=A0,A1=A1,D1=D1,targets=targets,M=M/area,area=area,free_rd=free_rd,free_c=free_c,zero_c=zero_c)

def metric_factor(a,N,mesh):
    nd=len(mesh['free'][1])*3;nc=len(mesh['free'][2])*3;nk=N.shape[1];nr=len(mesh['free'][0])*3
    P=np.zeros((nd+nc,nk+nc));P[:nd,:nk]=N[nr:];P[nd:,nk:]=np.eye(nc)
    L=np.linalg.cholesky(np.kron(a['M'],np.eye(3)))
    return L.T@P

def solve_ball(T,b,F,radius):
    U,s,V=svd(F,full_matrices=True,check_finite=False);del U
    ranks=lambda s:[int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]]
    rh=ranks(s)
    if rh!=[rh[1]]*3:raise ValueError('H1_metric_rank_uncertain:'+str(rh))
    r=rh[1];W=V[:r].T/s[:r];Z=V[r:].T;Ap=T@W;Af=T@Z
    Uf,ss,Vf=svd(Af,full_matrices=False,check_finite=False) if Z.shape[1] else (np.zeros((len(b),0)),np.array([]),np.zeros((0,0)))
    rr=ranks(ss) if len(ss) else [0]*3
    if rr!=[rr[1]]*3:raise ValueError('unpenalized_rank_uncertain:'+str(rr))
    fr=rr[1];Uf=Uf[:,:fr];Vf=Vf[:fr];ss=ss[:fr]
    A=Ap-Uf@(Uf.T@Ap);e=b-Uf@(Uf.T@b);U,sv,Vh=svd(A,full_matrices=False,check_finite=False);rr=ranks(sv)
    if rr!=[rr[1]]*3:raise ValueError('curvature_rank_uncertain:'+str(rr))
    ra=rr[1];rhs=U[:,:ra].T@e;sv=sv[:ra];Vh=Vh[:ra]
    def z(lam):return Vh.T@(sv*rhs/(sv*sv+lam))
    uncon=z(0.);lam=0.
    if np.linalg.norm(uncon)>radius:
        hi=1.
        while np.linalg.norm(z(hi))>radius:hi*=10
        lam=brentq(lambda x:np.linalg.norm(z(x))/radius-1,0.,hi,xtol=1e-16,rtol=1e-13)
    zp=z(lam);unpen=Vf.T@((Uf.T@(b-Ap@zp))/ss) if fr else np.zeros(Z.shape[1]);q=W@zp+Z@unpen;res=T@q-b
    kkt=max(float(np.linalg.norm(Af.T@res)/max(1.,np.linalg.norm(Af)*np.linalg.norm(b))),
        float(np.linalg.norm(Ap.T@res+lam*zp)/max(1.,np.linalg.norm(Ap)*np.linalg.norm(b)+lam*np.linalg.norm(zp))))
    return q,dict(error=float(np.linalg.norm(res)/np.linalg.norm(b)),h1=float(np.linalg.norm(F@q)),multiplier=float(lam),kkt=kkt,
        H1_ranks=rh,H1_nullity=Z.shape[1],unpenalized_rank=fr,unpenalized_nonuniqueness=Z.shape[1]-fr)
