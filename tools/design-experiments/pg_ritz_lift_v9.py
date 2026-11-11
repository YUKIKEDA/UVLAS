"""Reference-only membrane Ritz lift with exact first/second analytic jets."""
import numpy as np
from scipy.linalg import svd,qr,solve_triangular
import pg_ritz_geometry as geometry
import pg_ritz_particular_v2 as particular
import pg_ritz_clamp as clamp
import pg_ritz_lift_exact as parent
from pg_ritz_coarse import transform

def mandel(A):return np.stack((A[:,0,0],A[:,1,1],A[:,2,2],np.sqrt(2)*A[:,1,2],np.sqrt(2)*A[:,0,2],np.sqrt(2)*A[:,0,1]),axis=1)

def physical_map(p,parameters):
    r,A,n,dn,g=geometry.graph(p,parameters);metric=np.einsum('pai,paj->pij',A,A);T=A@np.linalg.inv(metric);t1,t2=T[:,:,0],T[:,:,1];V1=mandel(t1[:,:,None]*t1[:,None,:]);V2=mandel(t2[:,:,None]*t2[:,None,:]);V12=mandel((t1[:,:,None]*t2[:,None,:]+t2[:,:,None]*t1[:,None,:])/2)
    return V1,V2,V12,np.sqrt(g),A[:,:,0,],A[:,:,1]

def prepare(n,parameters,order=13):
    p,w=geometry.quadrature(n,order);raw=geometry.trial(p,n);A=raw[0]*np.sqrt(w)[:,None];U,s,V=svd(A,full_matrices=False,check_finite=False);ranks=[int(sum(s>tol*max(1.,np.linalg.norm(A)))) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('Ritz_trial_rank:'+str(ranks))
    T=V[:ranks[0]].T/s[:ranks[0]];R=[v@T for v in raw];W=geometry.seeds(p,n);V1,V2,V12,area,a1,a2=physical_map(p,parameters);fx=a1[:,2];fy=a2[:,2];B1=V1[:,:,None]*R[1][:,None,:]+V12[:,:,None]*R[2][:,None,:];B2=V2[:,:,None]*R[2][:,None,:]+V12[:,:,None]*R[1][:,None,:];E=V1[:,:,None]*(fx[:,None]*W[1])[:,None,:]+V2[:,:,None]*(fy[:,None]*W[2])[:,None,:]+V12[:,:,None]*(fx[:,None]*W[2]+fy[:,None]*W[1])[:,None,:];part=particular.jets(p,n,parameters);E[:,:,:16*n]=-V2[:,:,None]*parameters[1]*part[0][:,2,:16*n][:,None,:];E[:,:,16*n:]=-V1[:,:,None]*parameters[0]*part[0][:,2,16*n:][:,None,:];weight=np.sqrt(w*area);coeff=np.zeros((2,T.shape[1],32*n));metrics=[]
    for axis in [0,1]:
        edge=np.zeros((227,2));edge[:,1-axis]=np.linspace(0,1,227);jets=[v@T for v in geometry.trial(edge,n)];C=np.vstack((jets[0],jets[axis+1]));_,cs,Cv=svd(C,full_matrices=True,check_finite=False);cr=[int(sum(cs>tol*max(1.,np.linalg.norm(C)))) for tol in [1e-9,1e-10,1e-11]]
        if len(set(cr))!=1:raise ValueError('Ritz_seed_clamp_rank:'+str(cr))
        Z=Cv[cr[0]:].T;B=np.concatenate((np.einsum('pkn,nj->pkj',B1,Z,optimize=True),np.einsum('pkn,nj->pkj',B2,Z,optimize=True)),axis=2);A=(B*weight[:,None,None]).reshape(-1,2*Z.shape[1]);sl=slice(axis*16*n,(axis+1)*16*n);b=-(E[:,:,sl]*weight[:,None,None]).reshape(-1,16*n);scale=1/np.linalg.norm(A,axis=0);Q,F=qr(A*scale,mode='economic',check_finite=False);sv=svd(F,compute_uv=False,check_finite=False);fr=[int(sum(sv>tol*sv[0])) for tol in [1e-9,1e-10,1e-11]]
        if fr!=[A.shape[1]]*3:raise ValueError('Ritz_membrane_rank:'+str(fr))
        y=solve_triangular(F,Q.T@b,check_finite=False);u=scale[:,None]*y;coeff[0,:,sl]=Z@u[:Z.shape[1]];coeff[1,:,sl]=Z@u[Z.shape[1]:];res=A@u-b;stationarity=float(np.linalg.norm((A*scale).T@res)/max(1.,np.linalg.norm(b)));metrics.append(dict(axis=axis,columns=A.shape[1],rank=fr,support_rank=cr,rcond=float(sv[-1]/sv[0]),stationarity=stationarity,membrane_residual=float(np.linalg.norm(res)/max(1e-30,np.linalg.norm(b)))))
    model=dict(n=n,parameters=parameters,T=T,coeff=coeff,trial_ranks=ranks,metrics=metrics);model['clamp']=clamp.prepare(model,geometry);return model

def evaluate(model,p):
    v,d,ds=parent.evaluate(model,p)
    part=particular.jets(p,model['n'],model['parameters'])
    r,A,normal,dnormal,g=geometry.graph(p,model['parameters'])
    H=np.concatenate((A.transpose(0,2,1),normal[:,None,:]),axis=1)
    rhs=np.stack((-np.einsum('pa,pan->pn',normal,part[1]),-np.einsum('pa,pan->pn',normal,part[2]),np.zeros((len(p),32*model['n']))),axis=1)
    dp=np.linalg.solve(H,rhs);derivatives=[]
    kx,ky,tau=model['parameters']
    for axis,(i,j) in enumerate([(3,4),(4,5)]):
        dn=dnormal[:,:,axis]
        db=np.stack((-np.einsum('pa,pan->pn',dn,part[1])-np.einsum('pa,pan->pn',normal,part[i]),-np.einsum('pa,pan->pn',dn,part[2])-np.einsum('pa,pan->pn',normal,part[j]),np.zeros((len(p),32*model['n']))),axis=1)
        dH=np.zeros_like(H);dH[:,0,2]=kx if axis==0 else tau;dH[:,1,2]=tau if axis==0 else ky;dH[:,2]=dn
        derivatives.append(np.linalg.solve(H,db-dH@dp))
    t=transform(model['n']);T=np.zeros((32*model['n'],16));T[:16*model['n'],:8]=t[:16*model['n']];T[16*model['n']:,8:]=t[:16*model['n']];part=[a@T for a in part];dp=dp@T;derivatives=[a@T for a in derivatives]
    return [np.concatenate((a,b),axis=2) for a,b in zip(v,part)],np.concatenate((d,dp),axis=2),[np.concatenate((a,b),axis=2) for a,b in zip(ds,derivatives)]
