"""Independent clamped analytic displacement families on quadratic graph surfaces."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
import pg_ritz_geometry as geometry

def jets(p,parameters):
    x,y=p.T;kx,ky,tau=parameters;a=np.array([2.,4.,6.]);xx=x[:,None];i0=np.arctan(a*xx)/a;i1=np.log1p(a*a*xx*xx)/(2*a*a);i2=xx/(a*a)-np.arctan(a*xx)/(a**3);b=xx*i0-i1;C=(xx*xx*i0-i2)/2;k=1/(1+a*a*xx*xx);fx=kx*x+tau*y;fy=ky*y+tau*x;z=np.zeros_like(b)
    v=[np.stack((-kx*C-tau*y[:,None]*b,-ky*y[:,None]*b+tau*xx*b-2*tau*C,b),axis=1),np.stack((-fx[:,None]*i0,-fy[:,None]*i0+tau*b,i0),axis=1),np.stack((-tau*b,-ky*b,z),axis=1),np.stack((-kx*i0-fx[:,None]*k,-fy[:,None]*k,k),axis=1),np.stack((-tau*i0,-ky*i0,z),axis=1),np.stack((z,z,z),axis=1)]
    r,A,normal,dn,g=geometry.graph(p,parameters);H=np.concatenate((A.transpose(0,2,1),normal[:,None,:]),axis=1);rhs=np.stack((-np.einsum('pa,pan->pn',normal,v[1]),-np.einsum('pa,pan->pn',normal,v[2]),z),axis=1);d=np.linalg.solve(H,rhs);ds=[]
    for axis,(j1,j2) in enumerate([(3,4),(4,5)]):
        dH=np.zeros_like(H);dH[:,0,2]=kx if axis==0 else tau;dH[:,1,2]=tau if axis==0 else ky;dH[:,2]=dn[:,:,axis];db=np.stack((-np.einsum('pa,pan->pn',dn[:,:,axis],v[1])-np.einsum('pa,pan->pn',normal,v[j1]),-np.einsum('pa,pan->pn',dn[:,:,axis],v[2])-np.einsum('pa,pan->pn',normal,v[j2]),z),axis=1);ds.append(np.linalg.solve(H,db-dH@d))
    return v,d,ds

def fields(p3,h,parameters,variable=False):
    v,d,ds=jets(p3[:,:2],parameters);r,A,n,dn,g=geometry.graph(p3[:,:2],parameters);X,J,det,H,Ha=geometry.reference(p3[:,:2],p3[:,2],h,parameters,variable);z=H*p3[:,2]/2;eta=np.zeros((len(p3),3,3,3),dtype=np.result_type(p3,float));eta[:,1,1]=-parameters[1]*v[0][:,2]
    for alpha in range(2):
        for beta in range(2):
            e1=(np.einsum('pa,pan->pn',A[:,:,alpha],ds[beta])+np.einsum('pa,pan->pn',A[:,:,beta],ds[alpha])+np.einsum('pa,pan->pn',dn[:,:,alpha],v[beta+1])+np.einsum('pa,pan->pn',dn[:,:,beta],v[alpha+1]))/2;e2=(np.einsum('pa,pan->pn',dn[:,:,alpha],ds[beta])+np.einsum('pa,pan->pn',dn[:,:,beta],ds[alpha]))/2;eta[:,alpha,beta]+=z[:,None]*e1+z[:,None]**2*e2
    inv=np.linalg.inv(J);E=np.einsum('pai,pabt,pbj->pijt',inv,eta,inv);mandel=np.stack((E[:,0,0],E[:,1,1],E[:,2,2],np.sqrt(2)*E[:,1,2],np.sqrt(2)*E[:,0,2],np.sqrt(2)*E[:,0,1]),axis=1);return mandel,v[0]+z[:,None,None]*d
