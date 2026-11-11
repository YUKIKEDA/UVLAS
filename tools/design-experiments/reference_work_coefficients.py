"""Coefficient algebra for RW-2, without subtracting nearly equal point fields.

Decimal arithmetic is confined to geometry polynomial construction. FE operators,
SVD, QR and all nonlinear solves remain float64. No coefficients are thresholded.
"""
import itertools,math
from decimal import Decimal,localcontext
from fractions import Fraction
from functools import lru_cache
import numpy as np
from scipy.linalg import svd
import reference_work_ancf_v2 as parent
IDS=list(itertools.product(range(9),range(9),range(5)));LOOK={v:i for i,v in enumerate(IDS)}
NORM=np.sqrt([math.prod(2*k+1 for k in v) for v in IDS])
RTROWS=np.array([3*LOOK[v]+a for a,v in parent.RTID])

def polyadd(a,b,scale=1):
    out=a.copy()
    for k,v in b.items():out[k]=out.get(k,Decimal(0))+v*scale
    return out
def polymul(a,b):
    out={}
    for i,v in a.items():
        for j,w in b.items():
            k=tuple(x+y for x,y in zip(i,j));out[k]=out.get(k,Decimal(0))+v*w
    return out

@lru_cache(None)
def powers(n):
    # Exact expansion of x**n into ordinary (not orthonormal) Legendre polynomials.
    c=[Fraction(1)]
    for _ in range(n):
        d=[Fraction(0)]*(len(c)+1)
        for j,a in enumerate(c):
            d[j+1]+=a*Fraction(j+1,2*j+1)
            if j:d[j-1]+=a*Fraction(j,2*j+1)
        c=d
    return [(i,a) for i,a in enumerate(c) if a]

def toleg(poly):
    out={}
    for key,value in poly.items():
        if not value:continue
        for choices in itertools.product(*(powers(n) for n in key)):
            target=tuple(c[0] for c in choices);rat=math.prod(c[1] for c in choices)
            out[target]=out.get(target,Decimal(0))+value*Decimal(rat.numerator)/Decimal(rat.denominator)
    return out

@lru_cache(maxsize=192)
def geometry_coefficients(qbytes,independent=False):
    q=np.frombuffer(qbytes,dtype=np.float64).reshape(3,5,5,3)
    with localcontext() as ctx:
        ctx.prec=80 if independent else 60
        nodes=[Decimal.from_float(float(x)) for x in parent.m.LOBATTO];lag=[]
        for i,a in enumerate(nodes):
            c={(0,):Decimal(1)};den=Decimal(1)
            for j,b in enumerate(nodes):
                if i!=j:c=polymul(c,{(1,):Decimal(1),(0,):-b});den*=a-b
            lag.append({k:v/den for k,v in c.items()})
        fields={};modal=np.zeros((3,5,5,3))
        for field,cart in itertools.product(range(3),repeat=2):
            P={}
            for t,s in itertools.product(range(5),repeat=2):
                value=Decimal.from_float(float(q[field,t,s,cart]))
                if value:
                    for (i,),a in lag[s].items():
                        for (j,),b in lag[t].items():P[(i,j)]=P.get((i,j),Decimal(0))+value*a*b
            fields[field,cart]=P
            for (i,j),v in toleg(P).items():modal[field,i,j,cart]=float(v)
        # J columns as polynomials in (xi,eta,zeta,h).
        J={}
        for axis,cart in itertools.product(range(3),repeat=2):
            P={}
            for field in range(3):
                factor=[Decimal(1),Decimal(1)/2,Decimal(1)/8][field]
                for (i,j),value in fields[field,cart].items():
                    if axis==0 and i:P[(i-1,j,field,field)]=value*factor*i
                    if axis==1 and j:P[(i,j-1,field,field)]=value*factor*j
                    if axis==2 and field:P[(i,j,field-1,field)]=value*factor*field
            J[axis,cart]=P
        coeff=np.zeros((len(IDS),3,15))
        for axis,cart in itertools.product(range(3),repeat=2):
            u,v=(axis+1)%3,(axis+2)%3;a,b=(cart+1)%3,(cart+2)%3
            P=polyadd(polymul(J[u,a],J[v,b]),polymul(J[u,b],J[v,a]),-1)
            for hpower in range(5):
                part={key[:3]:value for key,value in P.items() if key[3]==hpower}
                for key,value in toleg(part).items():coeff[LOOK[key],axis,5*cart+hpower]=float(value)/NORM[LOOK[key]]
    coeff.setflags(write=False);modal.setflags(write=False);return coeff,modal

def scalar_basis(p,independent=False):
    v,g=parent.basis(p,IDS,independent);return v*NORM,g*NORM[None,None,:]

def derivative(coeff):
    # coeff: scalar_mode, spatial_row, field; divergence from the same coefficients.
    out=np.zeros((len(IDS),coeff.shape[-1]))
    for j,key in enumerate(IDS):
        for axis,n in enumerate(key):
            for k in range(n-1,-1,-2):
                target=list(key);target[axis]=k
                out[LOOK[tuple(target)]]+=math.sqrt((2*n+1)*(2*k+1))*coeff[j,axis]
    return out

def restrict(coeff,axis,sign):
    out=np.zeros_like(coeff)
    for i,key in enumerate(IDS):
        target=list(key);target[axis]=0;out[LOOK[tuple(target)]]+=math.sqrt(2*key[axis]+1)*sign**key[axis]*coeff[i]
    return out

def completion(coeff,rows,independent=False):
    flat=coeff.reshape(-1,coeff.shape[-1]);scale=max(float(np.max(np.linalg.norm(flat,axis=0))),1e-30)
    residual=flat.copy();residual[rows]=0.;residual/=scale
    U,s,V=svd(residual,full_matrices=False,lapack_driver='gesvd' if independent else 'gesdd',check_finite=False);rr=parent.classify(s)
    if len(set(rr))!=1:raise ValueError('coefficient_geometric_rank_ambiguous:'+str(rr))
    r=rr[0];E=(residual@V[:r].T)/s[:r] if r else np.zeros((len(flat),0))
    # Orthonormal complement to the known coordinate rows, evaluated as a polynomial.
    result=np.zeros((len(flat),len(rows)+r));result[rows,np.arange(len(rows))]=1;result[:,len(rows):]=E
    return result.reshape(*coeff.shape[:-1],-1),dict(ranks=rr,singular=s,scale=scale)

def mixed_coefficients(q,independent=False):
    C,modal=geometry_coefficients(np.asarray(q,dtype=np.float64).tobytes(),independent);S,meta=completion(C,RTROWS,independent)
    return S,derivative(S),meta,modal

def face_coefficients(q,axis,sign,independent=False):
    C,_=geometry_coefficients(np.asarray(q,dtype=np.float64).tobytes(),independent);C=restrict(C[:,axis,:],axis,sign)
    rows=np.array([LOOK[v] for v in parent.VID if v[axis]==0]);return completion(C,rows,independent)

def reference_position(modal,h):
    out=np.zeros((75,3))
    for j,(i,k,z) in enumerate(parent.VID):
        if z==0:out[j]=(modal[0,i,k]+h*h*modal[2,i,k]/24)/math.sqrt((2*i+1)*(2*k+1))
        if z==1:out[j]=h*modal[1,i,k]/(2*math.sqrt((2*i+1)*(2*k+1)*3))
        if z==2:out[j]=h*h*modal[2,i,k]/(12*math.sqrt((2*i+1)*(2*k+1)*5))
    return out

def potential_coefficients():
    W=np.zeros((len(IDS),104))
    for i,j in enumerate(parent.WID):key=parent.PID[j];W[LOOK[key],i]=1/NORM[LOOK[key]]
    G=np.zeros((len(IDS),3,104))
    for axis in range(3):
        for j,key in enumerate(IDS):
            for k in range(key[axis]-1,-1,-2):
                target=list(key);target[axis]=k;G[LOOK[tuple(target)],axis]+=math.sqrt((2*key[axis]+1)*(2*k+1))*W[j]
    return W,G
