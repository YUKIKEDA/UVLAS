from dataclasses import replace
import numpy as np
from scipy.linalg import solve
from pg_uvlm_coupled_step import Problem
import uvlm_wake_assembly_v2 as aero
class FreeProblem(Problem):
    def __init__(self):
        self.parent_cache={};self.mean_records=[];super().__init__()
    def decode(self,s):return replace(super().decode(s),convection_mode='free')
    def mean(self,s,tangent=False):
        model=self.model;g=model.geometry(self.coordinates(self.star));state=self.decode(s);A=[]
        for p,n in zip(g['P'],g['n']):
            B,W,_,_=aero.field_basis(model,g,state.wakes,p);A.append(n@(B+sum(np.split(W,self.R+1,axis=1))@model.T))
        A=np.array(A);gamma=solve(A,-g['n']@self.flow.U);s=s.copy();s[8:8+self.nb]=gamma;s[8+self.nb:8+2*self.nb]=gamma;s[-self.R:]=np.tile(model.T@gamma,self.R);state=self.decode(s);V=aero.convection(model,state,self.flow,'free')[0];X=state.wakes[0].vertices;h=self.h;R=np.r_[((X[1]-X[0])/h-(V[0]+V[1])/2).ravel(),((3*X[2]-4*X[1]+X[0])/(2*h)-V[2]).ravel()];J=None
        if tangent:
            J=np.empty((12,12))
            for k in range(12):
                dw=[np.zeros_like(X)];dw[0][1:].ravel()[k]=1.;rr=[]
                for i,p in enumerate(g['P']):rr.append(g['n'][i]@aero.field(model,g,state.wakes,gamma,p,dwakes=dw)[1])
                dg=solve(A,-np.array(rr));dm=np.tile(model.T@dg,self.R+1);dV=np.zeros_like(V)
                for r in range(3):
                    for j in range(2):dV[r,j]=aero.field(model,g,state.wakes,gamma,X[r,j],vertex_id=state.wakes[0].vertex_ids[r][j],dp=dw[0][r,j],dwakes=dw,dgamma=dg,dmu=dm)[1]
                J[:,k]=np.r_[(dw[0][1]/h-(dV[0]+dV[1])/2).ravel(),((3*dw[0][2]-4*dw[0][1])/(2*h)-dV[2]).ravel()]
        return R,J,s
    def equilibrium_state(self):
        s=super().equilibrium_state();offset=8+2*self.nb
        for it in range(10):
            R,J,s=self.mean(s,True);error=float(np.linalg.norm(R));self.mean_records.append(error)
            if error<1e-11:return s
            dx=solve(J,-R);ok=False
            for k in range(8):
                trial=s.copy();trial[offset:offset+12]+=dx/2**k
                try:rr,_,ss=self.mean(trial)
                except ValueError:continue
                if np.linalg.norm(rr)<error:s=ss;ok=True;break
            if not ok:raise ValueError('free_age_equilibrium_linesearch')
        raise ValueError('free_age_equilibrium_unresolved')
    def evaluate(self,x,s,tangent=False):
        model=self.model;h=self.h;alpha=1.5/h;r=x[:2];u=x[2:4];parent=self.decode(s);q=self.coordinates(r);v=np.einsum('nci,i->nc',self.B,u);g=model.geometry(q,v);wp=parent.wakes[0];info=model.surface_info[0];release=('pg','released',parent.epoch+1,0.,h);w=aero.Wake('pg',np.concatenate([g['X'][info['te']][None],x[4:].reshape(2,2,3)]),np.concatenate([np.zeros((1,1)),wp.mu[:2]]),(wp.vertex_ids[0],tuple((release,j) for j in range(2)),wp.vertex_ids[1]),(wp.row_ids[0],release,wp.row_ids[1]),(h,0.,-h));trial=dict(identity=model.identity,flow_identity=self.flow.identity,parent_epoch=parent.epoch,parent_fingerprint=parent.fingerprint(),time=h,h=h,q=q,v=v,wakes=(w,),mode='free',retired=(wp.row_ids[-1],),tail='finite_closed_window',max_intervals=2,qualified=False);ws=aero.Workspace();aero.solve_trial(model,parent,trial,self.flow,ws);gamma=trial['gamma'];w=trial['wakes'][0];gd=(3*gamma-4*parent.gamma+parent.gamma_previous)/(2*h);loads=aero.loads(model,g,(w,),gamma,gd,self.flow,h);force,K=self.internal(r,tangent);key=parent.fingerprint()
        if key not in self.parent_cache:self.parent_cache[key]=aero.convection(model,parent,self.flow,'free')[0][0]
        V0=self.parent_cache[key];V=np.empty((2,2,3))
        for row in range(1,3):
            for j in range(2):V[row-1,j]=self.flow.value(w.vertices[row,j],h)+aero.field(model,g,(w,),gamma,w.vertices[row,j],vertex_id=w.vertex_ids[row][j])[0]
        previous_te=model.geometry(self.coordinates(s[4:6]))['X'][info['te']];R=np.r_[(3*r-4*s[:2]+s[4:6])/(2*h)-u,self.mass@((3*u-4*s[2:4]+s[6:8])/(2*h))+force-self.fixed-self.flatB.T@loads['Q'].ravel(),((w.vertices[1]-wp.vertices[0])/h-(V0+V[0])/2).ravel(),((3*w.vertices[2]-4*wp.vertices[1]+previous_te)/(2*h)-V[1]).ravel()];out=np.r_[r,u,s[:4],gamma,parent.gamma,w.vertices[1:].ravel(),w.mu[1:].ravel()];J=None
        if tangent:
            J=np.zeros((16,16));J[:4,:4]=np.block([[alpha*np.eye(2),-np.eye(2)],[K,alpha*self.mass]])
            for k in range(16):
                dq=self.B[:,:,k] if k<2 else np.zeros_like(q);dv=self.B[:,:,k-2] if 2<=k<4 else np.zeros_like(v);gg=model.geometry(q,v,dq,dv);dw=[np.zeros_like(w.vertices)];dw[0][0]=gg['dX'][info['te']]
                if k>=4:dw[0][1:].ravel()[k-4]=1.
                rr=[]
                for i,p in enumerate(gg['P']):
                    ind,der=aero.field(model,gg,(w,),gamma,p,dp=gg['dP'][i],dwakes=dw);rr.append(gg['dn'][i]@(ind+self.flow.value(p,h)-gg['CV'][i])+gg['n'][i]@(der-gg['dCV'][i]))
                dg=solve(ws.A,-np.array(rr));dm=np.zeros(w.mu.size);dm[0]=model.T@dg;der=aero.loads(model,gg,(w,),gamma,gd,self.flow,h,dgamma=dg,dgamma_dot=alpha*dg,dwakes=dw,dmu=dm);J[2:4,k]-=self.flatB.T@der['dQ'].ravel();dV=np.empty_like(V)
                for row in range(1,3):
                    for j in range(2):dV[row-1,j]=aero.field(model,gg,(w,),gamma,w.vertices[row,j],vertex_id=w.vertex_ids[row][j],dp=dw[0][row,j],dwakes=dw,dgamma=dg,dmu=dm)[1]
                J[4:,k]=np.r_[(dw[0][1]/h-dV[0]/2).ravel(),(3*dw[0][2]/(2*h)-dV[1]).ravel()]
        return R,out,J
    def guess(self,s):return np.r_[s[:4],s[8+2*self.nb:8+2*self.nb+self.free_size]]
    def advance(self,s,guess=None):
        x=self.guess(s) if guess is None else np.array(guess,copy=True)
        for i in range(10):
            R,out,J=self.evaluate(x,s,True)
            if np.linalg.norm(R)<1e-11:return out,i,float(np.linalg.norm(R))
            dx=solve(J,-R);x+=dx
        raise ValueError('free_coupled_Newton_unresolved')
