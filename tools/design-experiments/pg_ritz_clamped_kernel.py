"""Fixed linear support kernel with an explicit full r/d/c trace lifting.

The NEW trial fields are H=(I-L T) S Z. Pre-lift S Z is not accepted as
homogeneous. All mass, gradient and work evaluations use H itself.
"""
import numpy as np
import pg_ritz_shell_nested as parent
import pg_ritz_trace_prelift as pretrace
lift=parent.lift

def prepare(n,parameters,order=13):
    original=parent.prepare(n,parameters,order)
    Z,ranks,report=pretrace.coefficients(original,parent)
    return dict(n=n,parameters=parameters,lift=original['lift'],parent=original,Z=Z,ns=Z.shape[1],pretrace=report,compression=original['compression'])

def fields(model,p):
    original=model['parent'];Z=model['Z'];base=parent.fields(original,p)
    edge=p.copy();edge[:,0]=0;boundary=parent.fields(original,edge)
    values=[a@Z for a in base];b0=boundary[0]@Z;by=boundary[2]@Z
    x=p[:,0][None,:,None]
    return [values[0]-(1-x)*b0,values[1]+b0,values[2]-(1-x)*by]

def kinematics(model,p3,h,variable=False):
    p,back=np.unique(p3[:,:2],axis=0,return_inverse=True);fs=fields(model,p);X,J,det,H,Ha=lift.geometry.reference(p3[:,:2],p3[:,2],h,model['parameters'],variable);z=H*p3[:,2]/2;za=Ha*p3[:,2,None]/2;parts=[v[:,back] for v in fs];N=parts[0][0]+z[:,None]*parts[0][1]+z[:,None]**2/2*parts[0][2];G=np.empty((len(p3),model['ns'],3))
    for a in range(2):G[:,:,a]=parts[a+1][0]+z[:,None]*parts[a+1][1]+z[:,None]**2/2*parts[a+1][2]+za[:,a,None]*(parts[0][1]+z[:,None]*parts[0][2])
    G[:,:,2]=H[:,None]/2*(parts[0][1]+z[:,None]*parts[0][2]);G=G@np.linalg.inv(J);MB=np.zeros((6,3,3));MB[0,0,0]=MB[1,1,1]=MB[2,2,2]=1
    for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2)
    B=np.einsum('kij,pnj->pkni',MB,G).reshape(len(p3),6,3*model['ns']);return N,G,B,det,X
