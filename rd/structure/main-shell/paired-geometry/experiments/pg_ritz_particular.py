"""Exact reference-only particular lift; Ritz corrects its remaining membrane strain."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
import curvature_work_profiles as curvature

def jets(p,n,parameters):
    kx,ky,tau=parameters;out=[[] for _ in range(6)]
    for axis in [0,1]:
        x=p[:,axis];y=p[:,1-axis];ka,kb=(kx,ky) if axis==0 else (ky,kx);k,I=curvature.evaluate(x,n,0,[0.,0.]);I0,I1,I2=I;b=x[:,None]*I0-I1;C=(x[:,None]**2*I0-I2)/2;fa=ka*x+tau*y;fb=kb*y+tau*x;z=np.zeros_like(b)
        local=[np.stack((-ka*C-tau*y[:,None]*b,-kb*y[:,None]*b+tau*x[:,None]*b-2*tau*C,b),axis=1),np.stack((-fa[:,None]*I0,-fb[:,None]*I0+tau*b,I0),axis=1),np.stack((-tau*b,-kb*b,z),axis=1),np.stack((-ka*I0-fa[:,None]*k,-fb[:,None]*k,k),axis=1),np.stack((-tau*I0,-kb*I0,z),axis=1),np.stack((z,z,z),axis=1)]
        for j in range(6):
            value=local[j if axis==0 else [0,2,1,5,4,3][j]]
            if axis==1:value=value[:,[1,0,2]]
            out[j].append(value)
    return [np.concatenate(v,axis=2) for v in out]
