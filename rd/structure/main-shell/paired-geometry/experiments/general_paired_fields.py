"""Independent continuum lifts on nonquadratic graphs; rational and polynomial curvature.

The reference fields use continuous one-dimensional moment formulas, not the
FE basis, its coefficient inverse, or any assembled operator.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from functools import lru_cache
from fractions import Fraction as F
import numpy as np
import polynomial_reference_surface as geometry
import paired_exact_polynomials as exact


def rational_moments(x,rates,maximum):
    x=x[:,None];a=np.asarray(rates)[None];z=a*x
    result=[np.arctan(z)/a,np.log1p(z*z)/(2*a*a)]
    for m in range(2,maximum+1):result.append(x**(m-1)/(a*a*(m-1))-result[m-2]/(a*a))
    small=abs(z)<.5
    for m in range(maximum+1):
        term=np.ones_like(z);series=term/(m+1)
        for j in range(1,45):term=term*(-z*z);series=series+term/(m+2*j+1)
        stable=x**(m+1)*series
        result[m]=np.where(small,stable,result[m])
    return result,1/(1+z*z)


@lru_cache(None)
def polynomial_moment_coefficients(maximum):
    from paired_coarse_polynomials import compose_univariate
    output=[]
    for k in range(8):
        P=compose_univariate(exact.legendre(k),F(-1),F(2))
        columns=[]
        for m in range(maximum+1):
            coef=np.zeros(k+m+2)
            for (i,j),v in P.items():coef[i+m+1]=float(v/F(i+m+1))
            columns.append(coef)
        output.append(columns)
    return output


def polynomial_moments(x,maximum):
    coefficients=polynomial_moment_coefficients(maximum)
    I=[np.column_stack([np.polynomial.polynomial.polyval(x,coefficients[k][m]) for k in range(8)]) for m in range(maximum+1)]
    k=np.polynomial.legendre.legvander(2*x-1,7)
    return I,k


def jets(p,parameters,kind='rational',rates=(2.,4.,6.)):
    x,y=p.T;coefficients=geometry.coefficients(parameters);maximum=max((i for i,j in coefficients),default=0)+1
    if kind=='rational':I,k=rational_moments(x,rates,maximum)
    elif kind=='polynomial':I,k=polynomial_moments(x,maximum)
    else:raise ValueError('curvature_kind')
    columns=k.shape[1];dtype=np.result_type(p,float)
    J=[np.zeros((len(p),columns),dtype=dtype) for _ in range(4)]
    Js=[np.zeros_like(J[0]) for _ in range(4)]
    for (i,j),a in coefficients.items():
        for derivative in range(min(j,3)+1):
            factor=float(a)
            for d in range(derivative):factor*=j-d
            factor=factor*y**(j-derivative)
            J[derivative]+=factor[:,None]*I[i]
            Js[derivative]+=factor[:,None]*I[i+1]
    f=geometry.values(p,parameters);fx=geometry.values(p,parameters,1,0);fy=geometry.values(p,parameters,0,1)
    fxx=geometry.values(p,parameters,2,0);fxy=geometry.values(p,parameters,1,1);fyy=geometry.values(p,parameters,0,2)
    width=float(geometry.chart(parameters)[1]);z=np.zeros_like(I[0]);b=x[:,None]*I[0]-I[1]
    v=[np.stack(((-f[:,None]*I[0]+J[0])/width,-x[:,None]*J[1]+Js[1],b),axis=1),
       np.stack((-fx[:,None]*I[0]/width,-J[1],I[0]),axis=1),
       np.stack(((-fy[:,None]*I[0]+J[1])/width,-x[:,None]*J[2]+Js[2],z),axis=1),
       np.stack(((-fxx[:,None]*I[0]-fx[:,None]*k)/width,-fy[:,None]*k,k),axis=1),
       np.stack((-fxy[:,None]*I[0]/width,-J[2],z),axis=1),
       np.stack(((-fyy[:,None]*I[0]+J[2])/width,-x[:,None]*J[3]+Js[3],z),axis=1)]
    for a in v:a[:,1]-=float(geometry.shear(parameters))*a[:,0]
    _,A,normal,dn,g=geometry.graph(p,parameters)
    H=np.concatenate((A.transpose(0,2,1),normal[:,None,:]),axis=1)
    rhs=np.stack((-np.einsum('pa,pan->pn',normal,v[1]),-np.einsum('pa,pan->pn',normal,v[2]),z),axis=1)
    d=np.linalg.solve(H,rhs);ds=[]
    for axis,(a,bj) in enumerate([(3,4),(4,5)]):
        dH=np.zeros_like(H);dH[:,0,2]=fxx if axis==0 else fxy;dH[:,1,2]=fxy if axis==0 else fyy;dH[:,2]=dn[:,:,axis]
        db=np.stack((-np.einsum('pa,pan->pn',dn[:,:,axis],v[1])-np.einsum('pa,pan->pn',normal,v[a]),
                     -np.einsum('pa,pan->pn',dn[:,:,axis],v[2])-np.einsum('pa,pan->pn',normal,v[bj]),z),axis=1)
        ds.append(np.linalg.solve(H,db-dH@d))
    return v,d,ds


def fields(p3,h,parameters,variable=False,kind='rational',rates=(2.,4.,6.)):
    p,back=np.unique(p3[:,:2],axis=0,return_inverse=True)
    v,d,ds=jets(p,parameters,kind,rates);r,A,n,dn,g=geometry.graph(p,parameters)
    X,J,det,H,Ha=geometry.reference(p3[:,:2],p3[:,2],h,parameters,variable);z=H*p3[:,2]/2
    eta=np.zeros((len(p3),3,3,d.shape[-1]),dtype=np.result_type(p3,float))
    # eta_ss and eta_sy vanish analytically; eta_yy is retained.
    eta[:,1,1]=v[2][back,1]+float(geometry.shear(parameters))*v[2][back,0]
    for a in range(2):
        for b in range(2):
            e1=(np.einsum('pi,pit->pt',A[:,:,a],ds[b])+np.einsum('pi,pit->pt',A[:,:,b],ds[a])
                +np.einsum('pi,pit->pt',dn[:,:,a],v[b+1])+np.einsum('pi,pit->pt',dn[:,:,b],v[a+1]))/2
            e2=(np.einsum('pi,pit->pt',dn[:,:,a],ds[b])+np.einsum('pi,pit->pt',dn[:,:,b],ds[a]))/2
            eta[:,a,b]+=z[:,None]*e1[back]+z[:,None]**2*e2[back]
    inv=np.linalg.inv(J);E=np.einsum('pai,pabt,pbj->pijt',inv,eta,inv)
    mandel=np.stack((E[:,0,0],E[:,1,1],E[:,2,2],np.sqrt(2)*E[:,1,2],np.sqrt(2)*E[:,0,2],np.sqrt(2)*E[:,0,1]),axis=1)
    return mandel,v[0][back]+z[:,None,None]*d[back]


class Coarse:
    @staticmethod
    def fields(model,shell,p3,h,variable=False):return fields(p3,h,model['parameters'],variable,'polynomial')
