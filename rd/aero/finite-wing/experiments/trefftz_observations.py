"""Distinct far-field observations; none silently changes the bound load law."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
def midpoint(te,mu,rho=1.2,eps=0.):
    te=np.asarray(te,float);mu=np.asarray(mu,float)
    if te.shape!=(len(mu)+1,3) or eps<0:raise ValueError('invalid_trefftz_input')
    yz=te[:,1:];delta=np.diff(yz,axis=0)
    if np.any(np.linalg.norm(delta,axis=1)==0):raise ValueError('zero_transverse_strip')
    kappa=np.r_[0.,mu]-np.r_[mu,0.];p=(yz[:-1]+yz[1:])/2;total=0.
    for j in range(0,len(mu),128):
        r=p[j:j+128,None]-yz[None];den=np.sum(r*r,axis=-1)+eps*eps
        if np.any(den==0):raise ValueError('nonincident_trefftz_contact')
        velocity_dot=np.sum(kappa[None]*np.einsum('ijk,ik->ij',r,delta[j:j+128])/den,axis=1)/(2*np.pi);total+=np.dot(mu[j:j+128],velocity_dot)
    return -rho*total/2
