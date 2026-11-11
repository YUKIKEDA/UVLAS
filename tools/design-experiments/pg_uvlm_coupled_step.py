import numpy as np
from scipy.linalg import solve
import connection_shell as pg
import pg_ritz_mechanics as mechanics
import uvlm_wake_assembly_v2 as aero
import uvlm_surface_bindings as binding
class Problem:
    def __init__(self):
        self.h=.04;self.R=2;pm=pg.flat_model(1,.08,1.,1.);self.data=pg.data(pm,E=1000.,rho=2.,order=13);data=self.data;N=data['N'];X=data['X'];x,y,z=X.T;targets=np.stack([np.column_stack([0*x,-2*y*z,y*y]),np.column_stack([-y*y*z,-2*x*y*z,x*y*y])],axis=-1);B=np.empty((len(data['q0']),3,2))
        for i in range(2):B[:,:,i]=solve(data['M'],N.T@((data['weights']*2)[:,None]*targets[:,:,i]))
        self.B=B;self.flatB=B.reshape(-1,2);self.projection_error=float(np.linalg.norm(np.einsum('pn,nci->pci',N,B)-targets)/np.linalg.norm(targets));assert self.projection_error<1e-10;self.mass=sum(B[:,c].T@data['M']@B[:,c] for c in range(3));self.model=aero.Model([binding.pg_shell('pg',1,1,pm,data['q0'])],.025);self.flow=aero.Flow([2,0,.12]);self.star=np.array([.08,-.03]);self.nb=self.model.nb;self.n=1;self.free_size=self.R*2*3;self.dimension=8+2*self.nb+self.free_size+self.R;self.fixed=np.zeros(2);self.sstar=self.equilibrium_state();self.fixed=self.internal(self.star)[0]-self.static_load(self.sstar)
    def coordinates(self,r):return self.data['q0']+np.einsum('nci,i->nc',self.B,r)
    def internal(self,r,tangent=False):
        U,f,K,diag=mechanics.evaluate(self.data,np.einsum('nci,i->nc',self.B,r),tangent);assert diag['min_det']>0
        return self.flatB.T@f.ravel(),None if K is None else self.flatB.T@K@self.flatB
    def decode(self,s):
        model=self.model;info=model.surface_info[0];q=self.coordinates(s[:2]);v=np.einsum('nci,i->nc',self.B,s[2:4]);g=model.geometry(q,v);o=8;gamma=s[o:o+self.nb];previous=s[o+self.nb:o+2*self.nb];o+=2*self.nb;free=s[o:o+self.free_size].reshape(self.R,2,3);muold=s[o+self.free_size:].reshape(self.R,1);rows=(('pg','attached'),)+tuple(('pg','canonical_age',r) for r in range(1,self.R+1));ids=(tuple(model.vertex_ids[i] for i in info['te']),)+tuple(tuple((row,j) for j in range(2)) for row in rows[1:]);w=aero.Wake('pg',np.concatenate([g['X'][info['te']][None],free]),np.concatenate([(model.T@gamma).reshape(1,1),muold]),ids,rows,tuple(-self.h*r for r in range(self.R+1)));return aero.State(model.identity,100,0.,q,v,gamma,(w,),previous,self.h,'explicit_state',self.flow.identity,'prescribed',self.R)
    def equilibrium_state(self):
        q=self.coordinates(self.star);g=self.model.geometry(q);te=g['X'][self.model.surface_info[0]['te']];free=te[None]+np.arange(1,self.R+1)[:,None,None]*self.h*self.flow.U;zero=np.zeros(self.nb);s=np.r_[self.star,[0.,0.],self.star,[0.,0.],zero,zero,free.ravel(),np.zeros(self.R)];state=self.decode(s);A=[]
        for p,n in zip(g['P'],g['n']):
            Bb,W,_,_=aero.field_basis(self.model,g,state.wakes,p);A.append(n@(Bb+sum(np.split(W,self.R+1,axis=1))@self.model.T))
        gamma=solve(np.array(A),-g['n']@self.flow.U);s[8:8+self.nb]=gamma;s[8+self.nb:8+2*self.nb]=gamma;s[-self.R:]=np.tile(self.model.T@gamma,self.R);self.model.validate_state(self.decode(s));return s
    def static_load(self,s):
        st=self.decode(s);g=self.model.geometry(st.q,st.v);out=aero.loads(self.model,g,st.wakes,st.gamma,np.zeros(self.nb),self.flow,0.);return self.flatB.T@out['Q'].ravel()
    def evaluate(self,x,s,tangent=False):
        h=self.h;alpha=1.5/h;r=x[:2];u=x[2:];parent=self.decode(s);q=self.coordinates(r);v=np.einsum('nci,i->nc',self.B,u);trial=aero.begin_trial(self.model,parent,q,v,h,self.flow);ws=aero.Workspace();aero.solve_trial(self.model,parent,trial,self.flow,ws);gamma=trial['gamma'];gd=(3*gamma-4*parent.gamma+parent.gamma_previous)/(2*h);g=self.model.geometry(q,v);loads=aero.loads(self.model,g,trial['wakes'],gamma,gd,self.flow,h);force,K=self.internal(r,tangent);res=np.r_[(3*r-4*s[:2]+s[4:6])/(2*h)-u,self.mass@((3*u-4*s[2:4]+s[6:8])/(2*h))+force-self.fixed-self.flatB.T@loads['Q'].ravel()];out=np.r_[r,u,s[:4],gamma,parent.gamma,trial['wakes'][0].vertices[1:].ravel(),trial['wakes'][0].mu[1:].ravel()];J=None
        if tangent:
            J=np.block([[alpha*np.eye(2),-np.eye(2)],[K,alpha*self.mass]])
            for k in range(4):
                dq=self.B[:,:,k] if k<2 else np.zeros_like(q);dv=self.B[:,:,k-2] if k>=2 else np.zeros_like(v);gg=self.model.geometry(q,v,dq,dv);dw=[np.zeros_like(trial['wakes'][0].vertices)];dw[0][0]=gg['dX'][self.model.surface_info[0]['te']];fixed=[]
                for i,p in enumerate(gg['P']):
                    ind,dind=aero.field(self.model,gg,trial['wakes'],gamma,p,dp=gg['dP'][i],dwakes=dw);fixed.append(gg['dn'][i]@(ind+self.flow.value(p,h)-gg['CV'][i])+gg['n'][i]@(dind-gg['dCV'][i]))
                dg=solve(ws.A,-np.array(fixed));dm=np.zeros(trial['wakes'][0].mu.size);dm[:self.n]=self.model.T@dg;der=aero.loads(self.model,gg,trial['wakes'],gamma,gd,self.flow,h,dgamma=dg,dgamma_dot=alpha*dg,dwakes=dw,dmu=dm);J[2:,k]-=self.flatB.T@der['dQ'].ravel()
        return res,out,J
    def advance(self,s,guess=None):
        x=s[:4].copy() if guess is None else np.array(guess,copy=True)
        for i in range(10):
            R,out,J=self.evaluate(x,s,True)
            if np.linalg.norm(R)<1e-11:return out,i,float(np.linalg.norm(R))
            dx=solve(J,-R);x+=dx
        raise ValueError('coupled_Newton_unresolved')
