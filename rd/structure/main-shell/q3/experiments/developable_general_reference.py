"""Continuous linear bending reference and manufactured continuum loads."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from functools import lru_cache
from decimal import Decimal,localcontext
import itertools
import numpy as np
LD=np.longdouble
def constants(dtype=float):
    nu=dtype('.3');mu=1/(2*(1+nu));lam=nu/((1+nu)*(1-2*nu));return lam,mu
def positions(kind,k,x,y,z):
    zero=np.zeros_like(x);X=np.stack((x-k*x*z,y,k*x*x/2+z),axis=-1)
    if kind=='quadratic':W=x*x/2;Wp=x;U=-k*x**3/3
    else:W=x**4/12;Wp=x**3/3;U=-k*x**5/15
    u=np.stack((U-z*Wp,zero,W-k*x*z*Wp),axis=-1);return X,u
def direct(kind,k,points):
    J=[];Q=[]
    for axis in range(3):
        p=points.astype(complex);p[:,axis]+=1e-25j;X,u=positions(kind,k,*p.T);J.append(X.imag/1e-25);Q.append(u.imag/1e-25)
    J=np.stack(J,axis=2);Q=np.stack(Q,axis=2);G=Q@np.linalg.inv(J);return J,(G+G.transpose(0,2,1))/2
def inverse(k,x,z):
    j=1+k*k*x*x-k*z;out=np.zeros((len(x),3,3),dtype=np.result_type(x,z));out[:,0,0]=1/j;out[:,0,2]=k*x/j;out[:,1,1]=1;out[:,2,0]=-k*x/j;out[:,2,2]=(1-k*z)/j;return out,j
def strain(kind,k,x,z):
    a=np.stack((np.ones_like(x),np.zeros_like(x),k*x),axis=1);v=np.ones_like(x) if kind=='quadratic' else x*x;j=1+k*k*x*x-k*z;return -z[:,None,None]*v[:,None,None]/j[:,None,None]*a[:,:,None]*a[:,None,:]
def stress(kind,k,x,z):
    E=strain(kind,k,x,z);lam,mu=constants(E.dtype.type);return lam*np.trace(E,axis1=1,axis2=2)[:,None,None]*np.eye(3,dtype=E.dtype)[None]+2*mu*E
def body(kind,k,x,z):
    lam,mu=constants(x.dtype.type);D=1+k*k*x*x;j=D-k*z;v=np.ones_like(x) if kind=='quadratic' else x*x;vp=np.zeros_like(x) if kind=='quadratic' else 2*x;a=np.stack((np.ones_like(x),np.zeros_like(x),k*x),axis=1);ap=np.broadcast_to(np.array([0,0,k],dtype=x.dtype),a.shape);T=lam*D[:,None,None]*np.eye(3,dtype=x.dtype)[None]+2*mu*a[:,:,None]*a[:,None,:];Tx=lam*(2*k*k*x)[:,None,None]*np.eye(3,dtype=x.dtype)[None]+2*mu*(ap[:,:,None]*a[:,None,:]+a[:,:,None]*ap[:,None,:]);f=-z*v/j;fx=-z*(vp*j-v*2*k*k*x)/j**2;fz=-v*D/j**2;Sx=fx[:,None,None]*T+f[:,None,None]*Tx;Sz=fz[:,None,None]*T;Ji,_=inverse(k,x,z);return -np.einsum('pij,pj->pi',Sx,Ji[:,0,:])-np.einsum('pij,pj->pi',Sz,Ji[:,2,:])
def energy(E):
    lam,mu=constants(E.dtype.type);return lam*np.trace(E,axis1=1,axis2=2)**2+2*mu*np.einsum('pij,pij->p',E,E)
@lru_cache(None)
def decimal_gauss(n):
    with localcontext() as ctx:
        ctx.prec=80;nodes=[];weights=[]
        def poly(x):
            pm=Decimal(1);p=x
            if n==1:return p,pm
            for degree in range(2,n+1):pm,p=p,((2*degree-1)*x*p-(degree-1)*pm)/degree
            return p,pm
        for initial in np.polynomial.legendre.leggauss(n)[0]:
            x=Decimal.from_float(float(initial))
            for _ in range(16):
                p,pm=poly(x);dp=n*(x*p-pm)/(x*x-1);step=p/dp;x-=step
                if abs(step)<Decimal('1e-75'):break
            p,pm=poly(x);dp=n*(x*p-pm)/(x*x-1);nodes.append(+x);weights.append(2/((1-x*x)*dp*dp))
        return tuple(nodes),tuple(weights)
def gauss(n):
    p,w=decimal_gauss(n);return np.array([LD(str(x)) for x in p]),np.array([LD(str(x)) for x in w])
def basis(points):
    x,y,z=points.T;columns=[];grad=[]
    for a,b,c in itertools.product(range(4),range(4),range(3)):
        factor=LD('.5') if c==2 else LD(1);columns.append(factor*x**a*y**b*z**c);row=[]
        for index,power in enumerate([a,b,c]):
            powers=[a,b,c];powers[index]-=1;row.append(factor*power*x**powers[0]*y**powers[1]*z**powers[2] if power else np.zeros(len(x),dtype=LD))
        grad.append(np.stack(row,axis=1))
    return np.array(columns).T,np.array(grad).transpose(1,0,2)
def weak(kind,k0,h0,orders):
    k,h=LD(str(k0)),LD(str(h0));ranges=[];weights=[]
    for axis,order in enumerate(orders):
        p,w=gauss(order)
        if axis==0:p=(p+1)/2;w=w/2
        elif axis==1:p=p/2;w=w/2
        else:p=p*h/2;w=w*h/2
        ranges.append(p);weights.append(w)
    points=np.array(list(itertools.product(*ranges)));w=np.array([a*b*c for a,b,c in itertools.product(*weights)]);x,y,z=points.T;Ji,j=inverse(k,x,z);phi,g=basis(points);grad=g@Ji;S=stress(kind,k,x,z);b=body(kind,k,x,z);measure=w*j;internal=np.einsum('p,pnb,pab->na',measure,grad,S);volume=np.einsum('p,pn,pa->na',measure,phi,b);outer=np.zeros_like(volume);reaction=np.zeros_like(volume)
    lam,mu=constants(LD);testnorm=np.einsum('p,pna->na',measure,mu*np.sum(grad*grad,axis=2)[:,:,None]+(lam+mu)*grad*grad);ref=measure@energy(strain(kind,k,x,z))
    for axis,sign in itertools.product(range(3),[-1,1]):
        others=[a for a in range(3) if a!=axis];pts=np.zeros((len(ranges[others[0]])*len(ranges[others[1]]),3),dtype=LD);pts[:,others]=np.array(list(itertools.product(*[ranges[a] for a in others])));pts[:,axis]=(0 if sign<0 else 1) if axis==0 else sign*(LD('.5') if axis==1 else h/2);fw=np.array([a*b for a,b in itertools.product(*[weights[a] for a in others])]);xx,yy,zz=pts.T;inv,det=inverse(k,xx,zz);area=sign*det[:,None]*inv[:,axis,:];fp,_=basis(pts);traction=np.einsum('pab,pb->pa',stress(kind,k,xx,zz),area);force=np.einsum('p,pn,pa->na',fw,fp,traction)
        if axis==0 and sign<0:reaction=force
        else:outer+=force
    residual=internal-volume-outer-reaction;den=np.sqrt(ref*testnorm);mask=den>0;scaled=float(np.max(abs(residual[mask])/den[mask]));relative=float(np.sqrt(np.sum(residual*residual)/np.sum(internal*internal)));allowed=np.array([a>0 for a,b,c in itertools.product(range(4),range(4),range(3))]);reaction_allowed=float(np.max(abs(reaction[allowed])))
    return dict(internal=internal,relative=relative,cauchy=scaled,allowed_reaction=reaction_allowed,reference=float(ref))
