"""Independent Decimal80 SE(3) energy, matrix exponential/logarithm route."""
from decimal import Decimal,localcontext
import math
import numpy as np
D=Decimal
def dec(x):return D(str(float(x)))
def sincos(x):
    s=x;c=D(1);ts=x;tc=D(1)
    for k in range(1,100):
        ts*=-x*x/D((2*k)*(2*k+1));tc*=-x*x/D((2*k-1)*(2*k));s+=ts;c+=tc
        if abs(ts)+abs(tc)<D('1e-90'):break
    return s,c
def skew(w):
    x,y,z=w;return np.array([[D(0),-z,y],[z,D(0),-x],[-y,x,D(0)]],dtype=object)
def identity():return np.array([[D(int(i==j)) for j in range(3)] for i in range(3)],dtype=object)
def series_exponential(W):
    result=identity();term=identity()
    for k in range(1,100):
        term=(term@W)/D(k);result+=term
        if max(abs(a) for a in term.ravel())<D('1e-85'):break
    return result
def energy(r,R,L,K,ref,e0,delta):
    # Re-orthogonalize input float matrices by using their rotation-vector log
    # from a skew generator supplied through scipy Rotation, then Decimal exp.
    from scipy.spatial.transform import Rotation
    Rs=[];rs=[]
    for i in range(2):
        base=np.array([dec(a) for a in Rotation.from_matrix(R[i]).as_rotvec()]);rot=series_exponential(skew(base))
        Rs.append(rot@series_exponential(skew(delta[6*i+3:6*i+6])))
        rs.append(np.array([dec(a) for a in r[i]])+delta[6*i:6*i+3])
    Q=Rs[0].T@Rs[1];c=(np.trace(Q)-1)/2;theta=D(str(math.acos(float(c))))
    for _ in range(8):
        s,co=sincos(theta);theta+=(co-c)/s
    s,co=sincos(theta);factor=theta/(2*s);W=(Q-Q.T)*factor
    # Solve the exponential translation Jacobian by convergent matrix series.
    J=identity();term=identity()
    for k in range(1,100):
        term=term@W/D(k+1);J+=term
        if max(abs(a) for a in term.ravel())<D('1e-85'):break
    rhs=Rs[0].T@(rs[1]-rs[0]);A=np.column_stack((J,rhs))
    for i in range(3):
        pivot=max(range(i,3),key=lambda j:abs(A[j,i]));A[[i,pivot]]=A[[pivot,i]];A[i]/=A[i,i]
        for j in range(3):
            if i!=j:A[j]-=A[j,i]*A[i]
    d=np.r_[A[:,3],np.array([W[2,1],W[0,2],W[1,0]])]
    e=d/dec(L)-np.array([dec(a) for a in ref])-np.array([dec(a) for a in e0]);C=np.array([[dec(a) for a in row] for row in K])
    return dec(L)*(e@(C@e))/2
def directional(r,R,L,K,ref,e0,direction):
    with localcontext() as ctx:
        ctx.prec=80;h=D('1e-16');d=np.array([dec(a) for a in direction]);zero=np.array([D(0)]*12)
        u0=energy(r,R,L,K,ref,e0,zero);up=energy(r,R,L,K,ref,e0,h*d);um=energy(r,R,L,K,ref,e0,-h*d)
        return float((up-um)/(2*h)),float((up-2*u0+um)/(h*h))
