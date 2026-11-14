"""Convex shared-position approximation with shear and H1 balls; research only."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import itertools
import numpy as np
from scipy.linalg import cho_factor,cho_solve,solve_triangular,svd
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
import check_q9_curved_bounded_lift as core

def assemble(name,n,h,observation=9,position_order=9):
    xy,cells,coords,_,_,_=core.geometry.make_mesh(name,n);_,reference,analytic=core.fields(name,xy);side=2*n+1;interior=np.array([j*side+i for j in range(1,side-1) for i in range(1,side-1)]);index={node:i for i,node in enumerate(interior)}
    points,w=core.quadrature(observation);nr=2*len(points);A=np.zeros((nr*n*n,6*len(xy)));W=np.zeros((len(A),3*len(interior)));M=np.zeros((len(xy),len(xy)));M0=np.zeros_like(M);load=np.zeros((3*len(xy),3));norm=np.zeros(3);volume=0.;area=0.;area0=0.
    hp,hw=leggauss(5);mp,mw=leggauss(position_order)
    for cell,(_,_,ids,_) in enumerate(cells):
        q=coords[ids].reshape(81);state=core.method.prepare(q,h,0.,9);J,E,B=core.method.evaluate(state,points);a=J[:,:,:2];dz=J[:,:,2]*2/h;G=a.transpose(0,2,1)@a;row=2*np.einsum('pia,kij,pj->pak',a,core.method.MB,dz);white=np.linalg.solve(np.linalg.cholesky(G),row)/np.linalg.norm(dz,axis=1)[:,None,None];weight=w*np.linalg.det(J);volume+=weight.sum();C=(np.einsum('pak,pkn->pan',white,B[:,:,:81])*np.sqrt(weight)[:,None,None]).reshape(nr,81);sl=slice(nr*cell,nr*(cell+1))
        for local,node in enumerate(ids):
            A[sl,6*node:6*node+3]+=C[:,9*local+3:9*local+6];A[sl,6*node+3:6*node+6]+=C[:,9*local+6:9*local+9]/core.LENGTH
            if node in index:W[sl,3*index[node]:3*index[node]+3]+=C[:,9*local:9*local+3]
        H=np.zeros((9,9));Q=np.zeros((9,9));f=np.zeros((9,3,3))
        for i,j in itertools.product(range(5),repeat=2):
            N,dx,dy=core.shapes.shape(hp[i],hp[j]);g=np.column_stack((dx,dy));Js=coords[ids,0,:].T@g;metric=Js.T@Js;aw=hw[i]*hw[j]*np.sqrt(np.linalg.det(metric));area+=aw;H+=aw*(np.outer(N,N)+core.LENGTH**2*g@np.linalg.inv(metric)@g.T)
        for i,j in itertools.product(range(position_order),repeat=2):
            N,dx,dy=core.shapes.shape(mp[i],mp[j]);Js=coords[ids,0,:].T@np.column_stack((dx,dy));aw=mw[i]*mw[j]*np.sqrt(np.linalg.det(Js.T@Js));area0+=aw;Q+=aw*np.outer(N,N);x,y=N@xy[ids];normal=core.geometry.surface(name,x,y)[3]
            for k,(wa,wb) in enumerate(core.WAVES):
                value=np.sin(wa*np.pi*(x+.5))*np.sin(wb*np.pi*(y+.5))*normal;f[:,:,k]+=aw*N[:,None]*value[None,:];norm[k]+=aw*(value@value)
        M[np.ix_(ids,ids)]+=H;M0[np.ix_(ids,ids)]+=Q
        dofs=(3*ids[:,None]+np.arange(3)).reshape(-1);load[dofs]+=f.reshape(27,3)
    rids=(3*interior[:,None]+np.arange(3)).reshape(-1);A/=np.sqrt(volume);W/=np.sqrt(volume);M/=area;M0/=area0;load/=area0;norm/=area0;Gr=np.kron(M0[np.ix_(interior,interior)],np.eye(3))
    return dict(A=A,W=W,M=M,G=Gr,b=load[rids],norm=norm,reference=reference,interior=interior,xy=xy,analytic=analytic)

class ReducedProblem:
    def __init__(self,a):
        self.A,self.Lp=core.whiten(a['A'],a['M']);self.Lr=np.linalg.cholesky(a['G']);self.C=solve_triangular(self.Lr,a['W'].T,lower=True,check_finite=False).T
        U,s,Vh=svd(self.A,full_matrices=False,check_finite=False);keep=s>np.finfo(float).eps*s[0];self.s=s[keep];self.Vh=Vh[keep];U=U[:,keep];self.P=U.T@self.C;perp=self.C-U@self.P;self.Q=perp.T@perp;self.z0=solve_triangular(self.Lr,a['b'],lower=True,check_finite=False);self.identity=np.eye(self.C.shape[1]);self.rank=int(keep.sum());self.spectral=float(s[0])
    def solve(self,k,bound,relative_tolerance):
        z0=self.z0[:,k];forcing=float(np.linalg.norm(self.C@z0));epsilon=relative_tolerance*forcing;s2=self.s**2
        def point(lam,alpha):
            weight=alpha/(s2+alpha);S=self.Q+self.P.T@(weight[:,None]*self.P);S=(S+S.T)/2
            z=cho_solve(cho_factor(self.identity+lam*S,lower=True,check_finite=False),z0,check_finite=False);rhs=self.P@z;coef=-self.s*rhs/(s2+alpha);vnorm=float(np.linalg.norm(coef));residual=float(np.sqrt(max(0.,z@self.Q@z+np.sum((weight*rhs)**2))))
            return z,coef,vnorm,residual
        calls=0
        def inner(lam):
            nonlocal calls
            calls+=1;p0=point(lam,0.)
            if p0[2]<=bound:return (*p0,0.)
            center=np.log(s2[0]);low=center-80;high=center
            def root(logalpha):return np.log(max(point(lam,np.exp(logalpha))[2],1e-300)/bound)
            while root(low)<0:low-=30
            while root(high)>0:high+=10
            alpha=np.exp(brentq(root,low,high,xtol=1e-10,rtol=1e-13));return (*point(lam,alpha),float(alpha))
        initial=inner(0.)
        if initial[3]<=epsilon:lam=0.;result=initial
        else:
            low=-35.;high=0.
            def root(loglam):return np.log(max(inner(np.exp(loglam))[3],1e-300)/epsilon)
            while root(high)>0:
                high+=np.log(10.)
                if high>np.log(1e24):raise RuntimeError('dual_shear_multiplier_limit')
            lam=float(np.exp(brentq(root,low,high,xtol=1e-10,rtol=1e-13)));result=inner(lam)
        z,coef,_,_,alpha=result;v=self.Vh.T@coef;r=solve_triangular(self.Lr.T,z,lower=False,check_finite=False);p=core.physical(v,self.Lp);residual=self.A@v+self.C@z;mu=lam*alpha
        station_r=z-z0+lam*(self.C.T@residual);station_p=lam*(self.A.T@residual)+mu*v
        kr=float(np.linalg.norm(station_r)/max(np.linalg.norm(z0)+lam*np.linalg.norm(self.C.T@residual),1e-30));kp=float(np.linalg.norm(station_p)/max(lam*self.spectral*np.linalg.norm(residual)+mu*np.linalg.norm(v),1e-30)) if lam else 0.
        gamma=float(np.linalg.norm(residual));pnorm=float(np.linalg.norm(v));primal=max(0.,gamma/epsilon-1,pnorm/bound-1);comp=max(abs(lam*(gamma**2-epsilon**2))/max(lam*epsilon**2,1e-30),abs(mu*(pnorm**2-bound**2))/max(mu*bound**2,1e-30))
        return dict(r=r,p=p,epsilon=epsilon,forcing=forcing,residual=gamma,relative_shear=gamma/forcing,h1=pnorm,bound=bound,lambda_shear=lam,mu_h1=mu,kkt_r=kr,kkt_p=kp,primal=float(primal),complementarity=float(comp),position_change=float(np.linalg.norm(z-z0)),dual_calls=calls)
