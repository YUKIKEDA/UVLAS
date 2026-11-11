"""Linear frequency specialization. Analytic age integration is separately qualified."""
import numpy as np
from scipy.linalg import solve
from numpy.polynomial.legendre import leggauss
from uvlm_steady_reference import geometry,matrix,linear_load

def finite(p,a,b,eps):
    dl=b-a;length=np.sqrt(np.sum(dl*dl,axis=-1));e=dl/length[:,None];r=p[:,None,:]-a[None,:,:];s=np.sum(r*e[None],axis=-1);t=s-length[None];perp=r-s[:,:,None]*e[None];k2=np.sum(perp*perp,axis=-1)+eps*eps;aa=np.sqrt(s*s+k2);bb=np.sqrt(t*t+k2);same=np.real(s*t)>0
    integral=np.empty_like(s);integral[same]=np.broadcast_to(length,s.shape)[same]*(s+t)[same]/(aa*bb*(s*bb+t*aa))[same];integral[~same]=(s/aa-t/bb)[~same]/k2[~same]
    return np.cross(e[None],perp)*integral[:,:,None]/(4*np.pi)
def semi(p,a,eps):
    r=p[:,None,:]-a[None];s=r[:,:,0];perp=r.copy();perp[:,:,0]=0;k2=np.sum(perp*perp,axis=-1)+eps*eps;norm=np.sqrt(s*s+k2);factor=np.empty_like(s);minus=np.real(s)<0;factor[minus]=1/(norm*(norm-s))[minus];factor[~minus]=((1+s/norm)/k2)[~minus]
    return np.cross([1.,0,0],perp)*factor[:,:,None]/(4*np.pi)
def H(p,te,eps):return finite(p,te[:-1],te[1:],eps)+semi(p,te[1:],eps)-semi(p,te[:-1],eps)
def wake_matrix(g,omega,eps,length=256.,order=12,h=None,velocity=20.):
    T=length/velocity
    if h is None:
        a=0.;nodes=[0.];minimum=np.min(np.linalg.norm(g['P']-g['X'][-1,:-1].mean(axis=0),axis=1))/velocity
        step=min(minimum,.001)
        while a<T:
            a=min(T,a+min(step,np.pi/max(abs(omega),1e-8)));nodes.append(a);step*=1.5
    else:
        intervals=round(T/h)
        if abs(intervals*h-T)>1e-10:raise ValueError('nonintegral age window')
        nodes=np.linspace(0,T,intervals+1)
    x,w=leggauss(order);lo=np.array(nodes[:-1]);hi=np.array(nodes[1:]);ages=((lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*x/2).ravel();weights=((hi-lo)[:,None]*w/2).ravel()
    if h is None:mu=np.exp(-1j*omega*ages)
    else:
        r=np.floor(ages/h).astype(int);s=ages/h-r;mu=(1-s)*np.exp(-1j*omega*h*r)+s*np.exp(-1j*omega*h*(r+1))
    W=np.zeros((len(g['P']),g['n']),complex);te=g['X'][-1]
    for a in range(0,len(ages),8):
        tau=ages[a:a+8];p=np.tile(g['P'],(len(tau),1)).astype(complex);p[:,0]-=np.repeat(velocity*tau,len(g['P']));p[:,0]+=1e-20j;derivative=H(p,te,eps).imag/1e-20;N=np.tile(g['normal'],(len(tau),1));normal=np.einsum('pij,pj->pi',derivative,N).reshape(len(tau),len(g['P']),g['n']);W+=np.einsum('q,qpn->pn',weights[a:a+8]*mu[a:a+8]*velocity,normal)
    return W
def run(m,n,k=.2,eps=1e-5,length=256.,order=12,motion='heave',axis=.375):
    g=geometry(m,n);omega=2*k*20/1.114285714285714;A=matrix(g,eps,0.);W=wake_matrix(g,omega,eps,length,order);A=A.astype(complex);A[:,-n:]+=W
    rhs=1j*omega*g['normal'][:,2] if motion=='heave' else -(20+1j*omega*(g['P'][:,0]-axis))*g['normal'][:,2]
    gamma=solve(A,rhs);G=gamma.reshape(m,n);horse=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel();F,M,bands=linear_load(g,horse);X=g['X'];areas=[];moments=[]
    for i in range(m):
        for j in range(n):
            q=X[[i,i,i+1,i+1],[j,j+1,j+1,j]];area=np.zeros(3);moment=np.zeros(3)
            for ids in [[0,2,1],[0,3,2]]:
                tri=q[ids];v=np.cross(tri[1]-tri[0],tri[2]-tri[0])/2;area+=v;moment+=np.cross(tri.mean(axis=0)-[axis,0,0],v)
            areas.append(area);moments.append(moment)
    Fc=F.copy();Mc=M.copy();F+=1.2j*omega*(gamma[:,None]*areas).sum(axis=0);M+=1.2j*omega*(gamma[:,None]*moments).sum(axis=0)
    return dict(g=g,omega=omega,gamma=gamma,F=F,M=M,circulatory_F=Fc,circulatory_M=Mc,residual=float(np.linalg.norm(A@gamma-rhs)/np.linalg.norm(rhs)),W=W)
