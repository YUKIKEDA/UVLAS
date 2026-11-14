"""General graph PG-Ritz scalar Cartesian position space, with unit normal thickness."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import svd
import pg_ritz_lift_v7 as lift

def native(p,n):
    x,y=p.T;knots=np.r_[np.zeros(3),np.arange(1,n)/n,np.ones(3)];B=BSpline(knots,np.eye(n+2),2);X=[B(x),B.derivative()(x)];Y=[B(y),B.derivative()(y)];R=[np.einsum('pi,pj->pij',X[dx],Y[dy]).reshape(len(p),-1) for dx,dy in [(0,0),(1,0),(0,1)]];result=[R]
    for degree in [2,1]:
        # Piecewise C0 tensor Lagrange, represented by one-dimensional B-splines.
        knots=np.r_[np.zeros(degree+1),np.repeat(np.arange(1,n)/n,degree),np.ones(degree+1)];N=degree*n+1;B=BSpline(knots,np.eye(N),degree);X=[B(x),B.derivative()(x)];Y=[B(y),B.derivative()(y)];result.append([np.einsum('pi,pj->pij',X[dx],Y[dy]).reshape(len(p),-1) for dx,dy in [(0,0),(1,0),(0,1)]])
    return result

def fields(model,p):
    n=model['n'];fs=native(p,n);r,A,normal,dn,g=lift.geometry.graph(p,model['parameters']);comp=model['normal_completion'];D=fs[1]
    normal_jets=[normal,dn[:,:,0],dn[:,:,1]]
    for j in range(3):fs[1][j]=np.column_stack((D[j],(normal_jets[j]-D[j]@comp['coeff'])@comp['V']/comp['s']))
    v,d,ds=lift.evaluate(model['lift'],p);raw=[np.stack((v[0].reshape(len(p),-1),d.reshape(len(p),-1))),np.stack((v[1].reshape(len(p),-1),ds[0].reshape(len(p),-1))),np.stack((v[2].reshape(len(p),-1),ds[1].reshape(len(p),-1)))];e=model['completion'];nr=fs[0][0].shape[1];out=[]
    for j,part in enumerate(raw):
        R,D,C=[field[j] for field in fs];rr=(part[0]-R@e['coeff'][:nr])@e['V']/e['s'];dd=(part[1]-D@e['coeff'][nr:])@e['V']/e['s'];arr=np.zeros((3,len(p),model['ns']));arr[0,:,:model['nr']]=R;arr[1,:,model['nr']:model['nr']+model['nd']]=D;arr[2,:,model['nr']+model['nd']:model['nr']+model['nd']+model['nc']]=C;arr[0,:,model['nr']+model['nd']+model['nc']:]=rr;arr[1,:,model['nr']+model['nd']+model['nc']:]=dd;out.append(arr)
    return out

def prepare(n,parameters,order=13):
    model=dict(n=n,parameters=parameters,lift=lift.prepare(n,parameters,order));p,w=lift.geometry.quadrature(n,order);fs=native(p,n);normal=lift.geometry.graph(p,parameters)[2];D=fs[1][0];weight=np.sqrt(w);coef=np.linalg.lstsq(D*weight[:,None],normal*weight[:,None],rcond=1e-13)[0];P=(normal-D@coef)*weight[:,None];U,s,V=svd(P,full_matrices=False,check_finite=False);rank=int(sum(s>1e-11));model['normal_completion']=dict(coeff=coef,V=V[:rank].T,s=s[:rank]);D=np.column_stack((D,(normal-D@coef)@V[:rank].T/s[:rank]));R=fs[0][0];v,d,ds=lift.evaluate(model['lift'],p);P=np.vstack((v[0].reshape(len(p),-1),d.reshape(len(p),-1)));A=np.zeros((2*len(p),R.shape[1]+D.shape[1]));A[:len(p),:R.shape[1]]=R;A[len(p):,R.shape[1]:]=D;weight=np.tile(weight,2);coef=np.linalg.lstsq(A*weight[:,None],P*weight[:,None],rcond=1e-13)[0];res=(P-A@coef)*weight[:,None];norm=np.linalg.norm(P*weight[:,None],axis=0);keep=norm>1e-11*max(norm);scale=1/np.maximum(norm[keep],1e-300);U,s,V=svd(res[:,keep]*scale,full_matrices=False,check_finite=False);ranks=[int(sum(s>tol*max(1.,np.linalg.norm(res[:,keep]*scale)))) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('general_PG_rank:'+str(ranks))
    rank=ranks[0];mapping=np.zeros((P.shape[1],rank));mapping[keep]=scale[:,None]*V[:rank].T;model.update(completion=dict(coeff=coef,V=mapping,s=s[:rank],ranks=ranks),nr=R.shape[1],nd=D.shape[1],nc=fs[2][0].shape[1],ns=R.shape[1]+D.shape[1]+fs[2][0].shape[1]+rank,discarded_zero_columns=int(sum(~keep)));return model

def kinematics(model,p3,h,variable=False):
    p,back=np.unique(p3[:,:2],axis=0,return_inverse=True);fs=fields(model,p);X,J,det,H,Ha=lift.geometry.reference(p3[:,:2],p3[:,2],h,model['parameters'],variable);z=H*p3[:,2]/2;za=Ha*p3[:,2,None]/2;parts=[v[:,back] for v in fs];N=parts[0][0]+z[:,None]*parts[0][1]+z[:,None]**2/2*parts[0][2];G=np.empty((len(p3),model['ns'],3))
    for a in range(2):G[:,:,a]=parts[a+1][0]+z[:,None]*parts[a+1][1]+z[:,None]**2/2*parts[a+1][2]+za[:,a,None]*(parts[0][1]+z[:,None]*parts[0][2])
    G[:,:,2]=H[:,None]/2*(parts[0][1]+z[:,None]*parts[0][2]);G=G@np.linalg.inv(J);MB=np.zeros((6,3,3));MB[0,0,0]=MB[1,1,1]=MB[2,2,2]=1
    for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2)
    B=np.einsum('kij,pnj->pkni',MB,G).reshape(len(p3),6,3*model['ns']);return N,G,B,det,X
