"""Analytic graph reference jets, physical normal thickness, polynomial seed jets."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.interpolate import BSpline
import stable_curvature_profiles as stable

def graph(p,parameters):
    x,y=p.T;kx,ky,tau=parameters;fx=kx*x+tau*y;fy=ky*y+tau*x;g=1+fx*fx+fy*fy;v=np.column_stack((-fx,-fy,np.ones(len(p))));n=v/np.sqrt(g)[:,None]
    vx=np.broadcast_to([-kx,-tau,0.],v.shape);vy=np.broadcast_to([-tau,-ky,0.],v.shape);gx=2*(kx*fx+tau*fy);gy=2*(tau*fx+ky*fy)
    nx=vx/np.sqrt(g)[:,None]-v*(gx/(2*g**1.5))[:,None];ny=vy/np.sqrt(g)[:,None]-v*(gy/(2*g**1.5))[:,None]
    r=np.column_stack((x,y,(kx*x*x+ky*y*y)/2+tau*x*y));a1=np.column_stack((np.ones(len(p)),np.zeros(len(p)),fx));a2=np.column_stack((np.zeros(len(p)),np.ones(len(p)),fy));return r,np.stack((a1,a2),axis=2),n,np.stack((nx,ny),axis=2),g

def univariate(x,n):
    """C1 piecewise polynomial primitive space, including exact second jets."""
    coef=stable.get(n,0.,0);cell=np.minimum((np.clip(x,0,1)*2*n).astype(int),2*n-1);s=2*(x*2*n-cell)-1;out=[]
    for order in range(3):
        values=np.empty((len(x),coef.shape[-1]))
        for j in np.unique(cell):
            c=np.polynomial.chebyshev.chebder(coef[j],m=order,axis=0);mask=cell==j;values[mask]=np.polynomial.chebyshev.chebvander(s[mask],len(c)-1)@c*(4*n)**order
        out.append(values)
    return out

def trial(p,n):
    x,y=p.T;knots=np.r_[np.zeros(3),np.arange(1,n)/n,np.ones(3)];B=BSpline(knots,np.eye(n+2),2);X=[B(x),B.derivative()(x),B.derivative(2)(x)];Y=[B(y),B.derivative()(y),B.derivative(2)(y)];jets=[]
    for dx,dy in [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]:
        jets.append(np.einsum('pi,pj->pij',X[dx],Y[dy]).reshape(len(p),-1))
    for axis in [0,1]:
        U=univariate(p[:,axis],n);other=p[:,1-axis]
        for power in [0,1]:
            for j,(dx,dy) in enumerate([(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]):
                along=dx if axis==0 else dy;across=dy if axis==0 else dx
                factor=other**power if across==0 else np.ones(len(p)) if across==1 and power==1 else np.zeros(len(p))
                jets[j]=np.column_stack((jets[j],U[along]*factor[:,None]))
    return jets

def seeds(p,n):
    """Integrate piecewise Legendre P7 curvature, independently of any load."""
    columns=[]
    for axis in [0,1]:
        x=p[:,axis];v=np.zeros((len(p),16*n));dx=np.zeros_like(v);ddx=np.zeros_like(v)
        # Gauss integrates the degree <=9 primitive exactly on each seed interval.
        nodes,weights=np.polynomial.legendre.leggauss(8)
        for cell in range(2*n):
            left=cell/(2*n);width=1/(2*n);upper=np.clip(x,left,left+width);length=upper-left;q=left+length[:,None]*(nodes+1)/2;P=np.polynomial.legendre.legvander(2*(q-left)/width-1,7);I0=np.einsum('pqk,q,p->pk',P,weights,length/2);I1=np.einsum('pq,pqk,q,p->pk',q,P,weights,length/2);sl=slice(8*cell,8*(cell+1));v[:,sl]=x[:,None]*I0-I1;dx[:,sl]=I0;inside=(x>=left)&((x<left+width) if cell<2*n-1 else x<=left+width);ddx[inside,sl]=np.polynomial.legendre.legvander(2*(x[inside]-left)/width-1,7)
        zero=np.zeros_like(v)
        columns.append([v,dx if axis==0 else zero,dx if axis==1 else zero,ddx if axis==0 else zero,zero,ddx if axis==1 else zero])
    return [np.column_stack((columns[0][j],columns[1][j])) for j in range(6)]

def quadrature(n,order):
    q,w=np.polynomial.legendre.leggauss(order);x=np.concatenate([(j+(q+1)/2)/(2*n) for j in range(2*n)]);wx=np.tile(w/(4*n),2*n);X,Y=np.meshgrid(x,x,indexing='ij');W,V=np.meshgrid(wx,wx,indexing='ij');return np.column_stack((X.ravel(),Y.ravel())),(W*V).ravel()

def reference(p,zeta,h,parameters,variable=False):
    r,A,n,dn,g=graph(p,parameters);H=h*(1+.1*p[:,0]+.05*p[:,1]) if variable else np.full(len(p),h);Hz=np.tile([.1*h,.05*h],(len(p),1)) if variable else np.zeros((len(p),2));z=H*zeta/2;J=np.empty((len(p),3,3));J[:,:,:2]=A+z[:,None,None]*dn+n[:,:,None]*(zeta[:,None]*Hz/2)[:,None,:];J[:,:,2]=H[:,None]*n/2
    return r+z[:,None]*n,J,np.linalg.det(J),H,Hz
