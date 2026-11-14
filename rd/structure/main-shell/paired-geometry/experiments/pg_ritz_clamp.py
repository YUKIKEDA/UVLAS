"""Algebraic reference clamp projector; retain boundary-zero identities explicitly."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.interpolate import BSpline

def spline(y,n,order=0):
    knots=np.r_[np.zeros(11),np.repeat(np.arange(1,2*n)/(2*n),9),np.ones(11)];B=BSpline(knots,np.eye(18*n+2),10);return B(y) if order==0 else B.derivative(order)(y)

def prepare(model,geometry):
    n=model['n'];x,w=np.polynomial.legendre.leggauss(13);y=np.concatenate([(j+(x+1)/2)/(2*n) for j in range(2*n)]);B=spline(y,n);result=[]
    for axis in [0,1]:
        p=np.zeros((len(y),2));p[:,1-axis]=y;R=[v@model['T'] for v in geometry.trial(p,n)];sl=slice(axis*16*n,(axis+1)*16*n);coef=model['coeff'][:,:,sl];v=np.stack((R[0]@coef[0],R[0]@coef[1]),axis=1);d=np.stack((R[axis+1]@coef[0],R[axis+1]@coef[1]),axis=1);a=np.linalg.lstsq(B,v.reshape(len(y),-1),rcond=1e-13)[0];b=np.linalg.lstsq(B,d.reshape(len(y),-1),rcond=1e-13)[0];result.append((a.reshape(len(a),2,16*n),b.reshape(len(b),2,16*n)))
    return result

def project(model,p,delta):
    n=model['n'];out=[v.copy() for v in delta]
    for axis in [0,1]:
        along=p[:,axis];other=p[:,1-axis];a,b=model['clamp'][axis];edge=[];edge_d=[]
        for order in range(3):
            B=spline(other,n,order);edge.append(np.einsum('pn,nat->pat',B,a));edge_d.append(np.einsum('pn,nat->pat',B,b))
        sl=slice(axis*16*n,(axis+1)*16*n);indices=[0,1,2,3,4,5] if axis==0 else [0,2,1,5,4,3];correction=[edge[0]+along[:,None,None]*edge_d[0],edge_d[0],edge[1]+along[:,None,None]*edge_d[1],np.zeros_like(edge[0]),edge_d[1],edge[2]+along[:,None,None]*edge_d[2]]
        for local,index in enumerate(indices):out[index][:,:,sl]-=correction[local]
        on=along==0
        for index in [0,axis+1,2-axis]:out[index][on,:,sl]=0.
    return out
