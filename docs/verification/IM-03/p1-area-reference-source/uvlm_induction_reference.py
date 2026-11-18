"""Separate regularized doublet AREA integral and projection-fan quadrature.

No analytic boundary/sheet expressions or AD from the candidate are used.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import math
import numpy as np

def integrate(point,x,eps,order=64,fan=False):
    x=np.asarray(x);point=np.asarray(point);e=x[1]-x[0];f=x[2]-x[0];cross=np.cross(e,f);jac=np.sqrt(np.sum(cross*cross));normal=cross/jac;q,w=np.polynomial.legendre.leggauss(order);q=(q+1)/2;w=w/2;chunks=[]
    if not fan:
        u,v=np.meshgrid(q,q,indexing='ij');bary=np.stack(((1-u)*(1-v),u,(1-u)*v),axis=-1).reshape(-1,3);yy=bary@x;weights=(w[:,None]*w*(1-u)*jac).ravel();chunks=[(point-yy,weights,bary)]
    else:
        # Fan construction is used for values only. Derivatives use fixed
        # reference barycentric quadrature above, including geometry weights.
        d=float(np.dot(point-x[0],normal));projection=point-d*normal;E=np.column_stack((e,f));uv=np.linalg.solve(E.T@E,E.T@(projection-x[0]));beta=np.r_[1-sum(uv),uv];length=max(np.linalg.norm(a-projection) for a in x);delta=np.sqrt(d*d+eps*eps)/length
        if delta==0:raise ValueError('singular_reference_fan')
        cuts=[0.];r=min(delta,1.)
        while r<1.:cuts.append(r);r*=2
        cuts.append(1.)
        for i,j in [(0,1),(1,2),(2,0)]:
            signed=np.dot(np.cross(x[i]-projection,x[j]-projection),normal)
            if abs(signed)<1e-16*float(jac):continue
            angular=(1-q[:,None])*np.eye(3)[i]+q[:,None]*np.eye(3)[j]
            for lo,hi in zip(cuts[:-1],cuts[1:]):
                r=lo+(hi-lo)*q;bary=(1-r[:,None,None])*beta[None,None,:]+r[:,None,None]*angular[None,:,:];weights=r[:,None]*(hi-lo)*w[:,None]*w[None,:]*signed;edge=(1-q[:,None])*(x[i]-projection)+q[:,None]*(x[j]-projection);R=d*normal[None,None,:]-r[:,None,None]*edge[None,:,:];chunks.append((R.reshape(-1,3),weights.ravel(),bary.reshape(-1,3)))
    parts=[]
    for R,weights,bary in chunks:
        r2=np.sum(R*R,axis=1);rn=R@normal;den=r2+eps*eps;kernel=((r2-2*eps*eps)[:,None]*normal-3*R*rn[:,None])/den[:,None]**2.5;weighted=np.einsum('p,pk,pa->pka',weights,bary,kernel)/(4*np.pi)
        if np.iscomplexobj(weighted):parts.append(np.sum(weighted,axis=0))
        else:parts.append(np.array([[math.fsum(weighted[:,k,a].tolist()) for a in range(3)] for k in range(3)]))
    if any(np.iscomplexobj(a) for a in parts):return sum(parts)
    return np.array([[math.fsum(a[k,c] for a in parts) for c in range(3)] for k in range(3)])

def jacobian(point,x,eps,order=192):
    values=np.r_[point,np.asarray(x).ravel(),eps];out=np.zeros((3,3,13))
    for j in range(13):
        z=values.astype(complex);z[j]+=1e-25j;out[:,:,j]=integrate(z[:3],z[3:12].reshape(3,3),z[12],order,False).imag/1e-25
    return out
