"""Smooth reference-only coarse fields, retaining the full native correction space."""
import numpy as np

def transform(n):
    x,w=np.polynomial.legendre.leggauss(9);P=np.polynomial.legendre.legvander(x,7);out=np.zeros((32*n,8))
    for j in range(2*n):
        s=(j+(x+1)/2)/(2*n);Q=np.polynomial.legendre.legvander(2*s-1,7);out[8*j:8*(j+1)]=(np.arange(8)*2+1)[:,None]/2*(P.T@(w[:,None]*Q))
    return out

def fields(model,shell,p3,h,variable=False):
    p,back=np.unique(p3[:,:2],axis=0,return_inverse=True);v,d,ds=shell.lift.evaluate(model['lift'],p);sl=slice(32*model['n'],None);v=[a[:,:,sl] for a in v];d=d[:,:,sl];ds=[a[:,:,sl] for a in ds];r,A,normal,dn,g=shell.lift.geometry.graph(p,model['parameters']);T=np.eye(16)[:,:8];v=[x@T for x in v];d=d@T;ds=[x@T for x in ds];eta0=np.zeros((len(p),2,2,8));eta0[:,1,1]=-model['parameters'][1]*v[0][:,2]
    X,J,det,H,Ha=shell.lift.geometry.reference(p3[:,:2],p3[:,2],h,model['parameters'],variable);z=H*p3[:,2]/2;eta=np.zeros((len(p3),3,3,8));eta[:,:2,:2]=eta0[back]
    for alpha in range(2):
        for beta in range(2):
            e1=(np.einsum('pa,pan->pn',A[:,:,alpha],ds[beta])+np.einsum('pa,pan->pn',A[:,:,beta],ds[alpha])+np.einsum('pa,pan->pn',dn[:,:,alpha],v[beta+1])+np.einsum('pa,pan->pn',dn[:,:,beta],v[alpha+1]))/2;e2=(np.einsum('pa,pan->pn',dn[:,:,alpha],ds[beta])+np.einsum('pa,pan->pn',dn[:,:,beta],ds[alpha]))/2;eta[:,alpha,beta]+=z[:,None]*e1[back]+z[:,None]**2*e2[back]
    inv=np.linalg.inv(J);E=np.einsum('pai,pabt,pbj->pijt',inv,eta,inv);mandel=np.stack((E[:,0,0],E[:,1,1],E[:,2,2],np.sqrt(2)*E[:,1,2],np.sqrt(2)*E[:,0,2],np.sqrt(2)*E[:,0,1]),axis=1);return mandel,v[0][back]+z[:,None,None]*d[back]
