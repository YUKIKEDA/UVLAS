"""Independent fields and boundary work for the GC-ANCF trial."""
import numpy as np
import geometry_compatible_ancf_trial as f
c=f.c;m=f.m
STANDARD_REFERENCE=f.old.reference
def rational_reference(q,p,h,k):
    N,_=m.shape(p,0);x=(N@q.reshape(75,3)[:25])[:,0];z=h*p[:,2]/2;alpha=np.array([2.,4.,6.]);xx=x[:,None];a=alpha[None,:]
    i0=np.arctan(a*xx)/a;i1=np.log1p(a*a*xx*xx)/(2*a*a);i2=xx/(a*a)-np.arctan(a*xx)/(a*a*a)
    b=xx*i0-i1;primitive=(xx*xx*i0-i2)/2;u=np.zeros((len(p),3,3));u[:,0]=-k*primitive-z[:,None]*i0;u[:,2]=b-k*xx*z[:,None]*i0
    tangent=np.column_stack((np.ones(len(x)),np.zeros(len(x)),k*x));curvature=1/(1+a*a*xx*xx)
    E=np.einsum('kij,pijt->pkt',m.MB,-z[:,None,None,None]*tangent[:,:,None,None]*tangent[:,None,:,None]*curvature[:,None,None,:]/(1+k*k*x*x-k*z)[:,None,None,None])
    return E,u
def select(mode):f.old.reference=STANDARD_REFERENCE if mode=='standard' else rational_reference
def interface(a,sol):
    n=a['n'];h=a['h'];coef=sol['coef'];beta=sol['beta'];sharing=0.;fixed=0.;constraint=0.;zvals=[-1.,-.31,0.,.47,1.]
    for cell,q in enumerate(a['geo']['refs']):
        j,i=divmod(cell,n)
        for axis,neighbor in [(0,cell+1 if i<n-1 else None),(1,cell+n if j<n-1 else None)]:
            if neighbor is None:continue
            p=np.array([[s,t,z] for t in np.linspace(-1,1,17) for z in zvals for s in [1.]]) if axis==0 else np.array([[s,1.,z] for s in np.linspace(-1,1,17) for z in zvals])
            other=p.copy();other[:,axis]=-1.;N=c.kinematics(n,cell,q,h,p,a['geo'])[0];M=c.kinematics(n,neighbor,a['geo']['refs'][neighbor],h,other,a['geo'])[0]
            sharing=max(sharing,float(np.linalg.norm(N-M)/max(1.,np.linalg.norm(N))))
        if i==0:
            p=np.array([[-1.,t,z] for t in np.linspace(-1,1,57) for z in zvals]);N=c.kinematics(n,cell,q,h,p,a['geo'])[0];_,U=f.bending(q,p,h,a['k'],n)
            value=(N@coef.reshape(a['ns'],9)).reshape(len(p),3,3)+U@beta;fixed=max(fixed,float(np.max(abs(value))))
            constraint=max(constraint,float(np.linalg.norm(N@a['Zs'])/max(1.,np.linalg.norm(N))))
    return dict(sharing=sharing,fixed=fixed,boundary_constraint=constraint)
def reference_coefficients(a):
    n=a['n'];geo=a['geo'];q0=np.zeros((a['ns'],3));nb=(n+2)**2;q0[:nb]=geo['coeff']
    if geo['rank']:q0[nb:nb+geo['rank']]=geo['s'][:,None]*geo['V'].T
    for cell,q in enumerate(geo['refs']):
        j,i=divmod(cell,n)
        for field,degree in [(1,2),(2,1)]:
            points=np.array([[s,t,0.] for t in np.linspace(-1,1,degree+1) for s in np.linspace(-1,1,degree+1)]);N,_=m.shape(points,0);value=N@q.reshape(75,3)[25*field:25*(field+1)]
            offset=a['nr'] if field==1 else a['nr']+a['nd'];ids=np.array([(degree*j+b)*(degree*n+1)+degree*i+x for b in range(degree+1) for x in range(degree+1)])
            q0[offset+ids]=value
    return q0
def anchor_consistency(a,sol,order=15,zorder=9):
    import curvature_work_profiles as cp
    error=np.zeros(3);reference=np.zeros(3)
    for cell,q in enumerate(a['geo']['refs']):
        p,w=c.volume_rule(q,a['n'],order,zorder);N,dN=m.shape(p,0);coords=q.reshape(75,3)[:25];x=(N@coords)[:,0];xp=dN[:,:,0]@coords[:,0];yp=dN[:,:,1]@coords[:,0];z=a['h']*p[:,2]/2;k=a['k']
        curv,I=cp.evaluate(x,a['n'],0,[k,0.]);i0=I[0];ux=np.zeros((len(p),3,curv.shape[1]));uz=np.zeros_like(ux)
        ux[:,0]=-k*x[:,None]*i0-z[:,None]*curv;ux[:,2]=i0-k*z[:,None]*(i0+x[:,None]*curv);uz[:,0]=-i0;uz[:,2]=-k*x[:,None]*i0
        parent=np.stack((ux*xp[:,None,None],ux*yp[:,None,None],uz*a['h']/2),axis=2);_,_,det,inv,_=f.f.base.geometry(q,a['h'],p)
        grad=np.einsum('pain,pij->pajn',parent,inv,optimize=True);actual=np.einsum('kab,pabn->pkn',m.MB,grad)
        E,_=f.bending(q,p,a['h'],k,a['n']);delta=(actual-E)@sol['beta'];target,_=f.old.reference(q,p,a['h'],k)
        error+=np.einsum('p,pkt,kl,plt->t',w*det,delta,m.C,delta);reference+=np.einsum('p,pkt,kl,plt->t',w*det,target,m.C,target)
    return np.sqrt(error/reference).tolist()
def face_rule(q,n,axis,sign,order,zorder):
    if axis==2:
        plane,w=c.plane_rule(q,n,order);return np.column_stack((plane,np.full(len(plane),sign))),w
    g,w=np.polynomial.legendre.leggauss(order);z,wz=np.polynomial.legendre.leggauss(zorder);p0=np.zeros((2,3));p0[:,axis]=sign;p0[:,1-axis]=[-1,1];N,_=m.shape(p0,0);xy=N@q.reshape(75,3)[:25];cuts=[-1.,1.]
    for coordinate in [0,1]:
        delta=xy[1,coordinate]-xy[0,coordinate]
        if abs(delta)>1e-13:
            levels=np.arange(1,2*n)/(2*n)-(.5 if coordinate else 0.)
            cuts.extend(2*(v-xy[0,coordinate])/delta-1 for v in levels if min(xy[:,coordinate])+1e-13<v<max(xy[:,coordinate])-1e-13)
    cuts=np.unique(cuts);points=[];weights=[]
    for left,right in zip(cuts[:-1],cuts[1:]):
        for t,wt in zip((left+right)/2+(right-left)*g/2,w*(right-left)/2):
            for zz,ww in zip(z,wz):
                p=np.zeros(3);p[axis]=sign;p[1-axis]=t;p[2]=zz;points.append(p);weights.append(wt*ww)
    return np.array(points),np.array(weights)
def mechanics(a,order=15,zorder=9):
    q0=reference_coefficients(a);F=np.array([[1.03,.02,-.01],[.01,.98,.03],[0.,-.02,1.01]]);inc=q0@(F-np.eye(3)).T;E0=np.einsum('kij,ij->k',m.MB,(F.T@F-np.eye(3))/2);S0=np.einsum('kij,k->ij',m.MB,m.C@E0);P0=F@S0
    theta=.37;R=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1.]])
    rigid=q0@(R-np.eye(3)).T;force=np.zeros_like(q0);boundary=np.zeros_like(q0);geometry=0.;affine=0.;rigid_error=0.;derivative=0.;csderivative=0.
    d1=np.sin(np.arange(q0.size).reshape(q0.shape)*.19)*1e-4;d2=np.cos(np.arange(q0.size).reshape(q0.shape)*.23)*1e-4;epsilon=1e-30
    for cell,q in enumerate(a['geo']['refs']):
        p,w=c.volume_rule(q,a['n'],order,zorder);N,G,B,det,X=c.kinematics(a['n'],cell,q,a['h'],p,a['geo']);weight=w*det
        geometry=max(geometry,float(np.max(abs(N@q0-X))));Fp=np.eye(3)+np.einsum('na,pnb->pab',inc,G);Fr=np.eye(3)+np.einsum('na,pnb->pab',rigid,G)
        affine=max(affine,float(np.max(abs(Fp-F))));rigid_error=max(rigid_error,float(np.max(abs(Fr.transpose(0,2,1)@Fr-np.eye(3)))))
        E=np.einsum('kij,pij->pk',m.MB,(Fp.transpose(0,2,1)@Fp-np.eye(3))/2);S=np.einsum('kij,pk->pij',m.MB,E@m.C.T);P=Fp@S
        force+=np.einsum('p,pnb,pab->na',weight,G,P,optimize=True)
        df1=np.einsum('na,pnb->pab',d1,G);df2=np.einsum('na,pnb->pab',d2,G)
        e1=np.einsum('kij,pij->pk',m.MB,(Fp.transpose(0,2,1)@df1+df1.transpose(0,2,1)@Fp)/2)
        e2=np.einsum('kij,pij->pk',m.MB,(Fp.transpose(0,2,1)@df2+df2.transpose(0,2,1)@Fp)/2)
        sec=(df1.transpose(0,2,1)@df2+df2.transpose(0,2,1)@df1)/2
        derivative+=np.einsum('p,pk,kl,pl->',weight,e1,m.C,e2)+np.einsum('p,pij,pij->',weight,S,sec)
        Fc=Fp+1j*epsilon*df2;Ec=np.einsum('kij,pij->pk',m.MB,(Fc.transpose(0,2,1)@Fc-np.eye(3))/2);Pc=Fc@np.einsum('kij,pk->pij',m.MB,Ec@m.C.T)
        fc=np.einsum('p,pnb,pab->na',weight,G,Pc,optimize=True);csderivative+=float(np.imag(np.sum(fc*d1))/epsilon)
        for axis in range(3):
            for sign in [-1,1]:
                pp,ww=face_rule(q,a['n'],axis,sign,order,zorder);NN=c.kinematics(a['n'],cell,q,a['h'],pp,a['geo'])[0];_,_,dd,inv,_=f.f.base.geometry(q,a['h'],pp);area=sign*dd[:,None]*inv[:,axis,:]
                boundary+=np.einsum('p,pn,ab,pb->na',ww,NN,P0,area,optimize=True)
    return dict(reference_geometry=geometry,affine_gradient=affine,rigid_strain=rigid_error,constant_stress_work=float(np.linalg.norm(force-boundary)/max(1e-30,np.linalg.norm(force),np.linalg.norm(boundary))),nonlinear_tangent=float(abs(derivative-csderivative)/max(1e-30,abs(derivative),abs(csderivative))))
