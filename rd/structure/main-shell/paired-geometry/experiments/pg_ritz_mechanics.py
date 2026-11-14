"""Same PG-Ritz S coupled to layered 3D Green/second-PK material kernels."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import solve
import pg_ritz_shell_v3 as shell
import pg_ritz_trace_v2 as trace
MB=np.zeros((6,3,3))
for i in range(3):MB[i,i,i]=1.
for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2)
S0=np.diag([.0005,.0004,.0002]);FSTAR=np.array([[1.02,.015,-.007],[.003,.99,.011],[.004,-.006,1.01]])

def material():
    H=np.diag([1/90,1/9,1/6,1/2.5,1/3.5,1/4]);H[0,1]=H[1,0]=-.27/90;H[0,2]=H[2,0]=-.23/90;H[1,2]=H[2,1]=-.29/9;D=np.diag([1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2)]);return D@solve(H,np.eye(6))@D

def reference_coefficients(model):
    p,w=shell.lift.geometry.quadrature(model['n'],13);R=shell.native(p,model['n'])[0][0];r=shell.lift.geometry.graph(p,model['parameters'])[0];q=np.zeros((model['ns'],3));q[:model['nr']]=np.linalg.lstsq(R,r,rcond=1e-13)[0];nc=model['normal_completion'];nd0=(2*model['n']+1)**2;q[model['nr']:model['nr']+nd0]=nc['coeff'];q[model['nr']+nd0:model['nr']+model['nd']]=nc['s'][:,None]*nc['V'].T;return q

def prepare(model,h=.05,split=False,order=13,zorder=5):
    plane,wp=shell.lift.geometry.quadrature(model['n'],order);layers=[(-1.,0.,.23,1.,1.2),(0.,1.,-.41,.7,.8)]
    if split:layers=[(lo,(lo+hi)/2,theta,factor,rho) for lo,hi,theta,factor,rho in layers]+[((lo+hi)/2,hi,theta,factor,rho) for lo,hi,theta,factor,rho in layers]
    nodes,w=np.polynomial.legendre.leggauss(zorder);points=[];weights=[];angles=[];factors=[];densities=[]
    for lo,hi,theta,factor,rho in layers:
        z=(lo+hi)/2+(hi-lo)*nodes/2;points.append(np.column_stack((np.repeat(plane,len(z),axis=0),np.tile(z,len(plane)))));weights.append((wp[:,None]*w*(hi-lo)/2).ravel());angles.extend([theta]*(len(plane)*len(z)));factors.extend([factor]*(len(plane)*len(z)));densities.extend([rho]*(len(plane)*len(z)))
    p=np.vstack(points);weights=np.concatenate(weights);angles=np.array(angles);factors=np.array(factors);density=np.array(densities);N,G,B,det,X=shell.kinematics(model,p,h,True);_,J,_,H,Ha=shell.lift.geometry.reference(p[:,:2],p[:,2],h,model['parameters'],True);normal=np.cross(J[:,:,0],J[:,:,1]);normal/=np.linalg.norm(normal,axis=1)[:,None];fiber=np.column_stack((np.cos(angles),np.sin(angles),np.full(len(p),.1)));a=fiber-normal*np.sum(fiber*normal,axis=1)[:,None];a/=np.linalg.norm(a,axis=1)[:,None];A=np.stack((a,np.cross(normal,a),normal),axis=2);T=np.einsum('kij,pia,lab,pjb->pkl',MB,A,MB,A,optimize=True);Cm=material()[None]*factors[:,None,None];Cg=T@Cm@T.transpose(0,2,1);estar=np.einsum('kij,ij->k',MB,(FSTAR.T@FSTAR-np.eye(3))/2);s0=np.einsum('kij,ij->k',MB,S0);e0=estar-np.linalg.solve(Cg,np.broadcast_to(s0,(len(p),6))[:,:,None])[:,:,0];M=N.T@((weights*det*density)[:,None]*N);return dict(model=model,h=h,p=p,weights=weights*det,N=N,G=G,C=Cg,e0=e0,A=A,T=T,rho=density,M=M,X=X,min_det=float(det.min()),natural_metric_min=float(np.linalg.eigvalsh(np.eye(3)+2*np.einsum('pk,kij->pij',e0,MB)).min()))

def evaluate(data,u,tangent=False):
    G=data['G'];weight=data['weights'];F=np.eye(3)+np.einsum('na,pnb->pab',u,G,optimize=True);E=np.einsum('kij,pij->pk',MB,(F.transpose(0,2,1)@F-np.eye(3))/2);e=E-data['e0'];s=np.einsum('pij,pj->pi',data['C'],e);S=np.einsum('pk,kij->pij',s,MB);P=F@S;U=np.sum(weight*np.sum(e*s,axis=1))/2;force=np.einsum('p,pnb,pab->na',weight,G,P,optimize=True);K=None
    if tangent:
        ns=len(u);K=np.zeros((3*ns,3*ns),dtype=np.result_type(u,float))
        for start in range(0,len(G),256):
            sl=slice(start,start+256);B=np.einsum('kij,pai,pnj->pkna',MB,F[sl],G[sl],optimize=True).reshape(len(G[sl]),6,3*ns);CB=np.einsum('pij,pjn->pin',data['C'][sl],B,optimize=True);K+=np.einsum('p,pki,pkj->ij',weight[sl],B,CB,optimize=True)
        geometric=np.einsum('p,pni,pij,pmj->nm',weight,G,S,G,optimize=True);K+=np.kron(geometric,np.eye(3))
    return U,force,K,dict(min_det=float(np.linalg.det(F.real).min()),P=P,F=F)

def external(model,h,order=15,zorder=8):
    P=FSTAR@S0;out=np.zeros((model['ns'],3));fixed=np.zeros_like(out);plane,wp=shell.lift.geometry.quadrature(model['n'],order);g,w=np.polynomial.legendre.leggauss(order);edge=np.concatenate([(j+(g+1)/2)/(2*model['n']) for j in range(2*model['n'])]);ew=np.tile(w/(4*model['n']),2*model['n']);z,wz=np.polynomial.legendre.leggauss(zorder)
    for axis in range(3):
        for sign in [-1,1]:
            if axis==2:p=np.column_stack((plane,np.full(len(plane),sign)));weight=wp
            else:
                p=np.zeros((len(edge)*len(z),3));p[:,axis]=(sign+1)/2;p[:,1-axis]=np.repeat(edge,len(z));p[:,2]=np.tile(z,len(edge));weight=(ew[:,None]*wz).ravel()
            N,G,B,det,X=shell.kinematics(model,p,h,True);_,J,_,H,Ha=shell.lift.geometry.reference(p[:,:2],p[:,2],h,model['parameters'],True);area=sign*det[:,None]*np.linalg.inv(J)[:,axis,:];F=np.einsum('p,pn,ab,pb->na',weight,N,P,area,optimize=True)
            if axis==0 and sign<0:fixed+=F
            else:out+=F
    return out,fixed
