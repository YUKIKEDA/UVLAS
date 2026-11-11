"""Exact all-component midsurface isometries on the intended mapped cylinder.

This geometry is the analytic cylinder through the binary-rational x/y corner
coordinates. It is not asserted identical to a rounded nodal z interpolant at h=0.
"""
from fractions import Fraction as F
from functools import lru_cache
import itertools
import hashlib
import numpy as np
import reference_work_thin_limit as thin
MID=[(i,j) for i,j in itertools.product(range(6),repeat=2) if (i,j)!=(5,5)]

def add(a,b,factor=F(1)):
    out=a.copy()
    for k,v in b.items():out[k]=out.get(k,F(0))+factor*v
    return {k:v for k,v in out.items() if v}
def mul(a,b):
    out={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():out[i+k,j+l]=out.get((i+k,j+l),F(0))+v*w
    return {k:v for k,v in out.items() if v}
def deriv(a,axis):
    out={}
    for key,v in a.items():
        if key[axis]:target=list(key);target[axis]-=1;out[tuple(target)]=v*key[axis]
    return out
def bilinear(values):
    out={}
    for value,(sx,sy) in zip(values,[(-1,-1),(1,-1),(-1,1),(1,1)]):
        v=F.from_float(float(value))/4
        for key,factor in [((0,0),1),((1,0),sx),((0,1),sy),((1,1),sx*sy)]:out[key]=out.get(key,F(0))+factor*v
    return {k:v for k,v in out.items() if v}
def exact_null(rows,n):
    rows=[r.copy() for r in rows if r];pivots=[];active=0
    for col in range(n):
        pivot=next((i for i in range(active,len(rows)) if rows[i].get(col,F(0))),None)
        if pivot is None:continue
        rows[active],rows[pivot]=rows[pivot],rows[active];factor=rows[active][col];rows[active]={k:v/factor for k,v in rows[active].items()};p=rows[active]
        for i in range(len(rows)):
            if i==active:continue
            factor=rows[i].get(col,F(0))
            if factor:
                for k,v in p.items():
                    value=rows[i].get(k,F(0))-factor*v
                    if value:rows[i][k]=value
                    else:rows[i].pop(k,None)
        pivots.append(col);active+=1
    free=[i for i in range(n) if i not in pivots];Z=[]
    for j in free:
        vector=[F(0)]*n;vector[j]=F(1)
        for row,pivot in zip(rows[:active],pivots):vector[pivot]=-row.get(j,F(0))
        Z.append(vector)
    return Z,pivots

@lru_cache(maxsize=96)
def certificate(qbytes,kappa):
    q=np.frombuffer(qbytes,dtype=np.float64).reshape(75,3);corners=q[[0,4,20,24]];x=bilinear(corners[:,0]);y=bilinear(corners[:,1]);z={k:F(str(kappa))*v/2 for k,v in mul(x,x).items()};r=[x,y,z]
    tangent=[[deriv(v,axis) for v in r] for axis in range(2)];equations={}
    for i,(a,b) in enumerate(MID):
        mon={(a,b):F(1)};dx=deriv(mon,0);dy=deriv(mon,1)
        for cart in range(3):
            fields=[mul(tangent[0][cart],dx),mul(tangent[1][cart],dy),add(mul(tangent[0][cart],dy),mul(tangent[1][cart],dx))]
            for eq,poly in enumerate(fields):
                for key,value in poly.items():equations.setdefault((eq,key),{})[3*i+cart]=value
    rows=list(equations.values());vectors,pivots=exact_null(rows,105)
    # Exact rational substitution certifies both rank and the complete nullspace.
    residual=all(sum(value*vector[col] for col,value in row.items())==0 for row in rows for vector in vectors)
    membrane_vectors=vectors;membrane_pivots=pivots
    # On J_top/bottom=0, the normal reference-gradient image is exactly
    # d_zeta W = (Q55 minus its corner) x P1(zeta). A1 q0 is zeta-independent;
    # hence solve the remaining director compatibility as a polynomial identity.
    director=[{key:-F(str(kappa))*value for key,value in x.items()},{},{(0,0):F(1)}]
    def cross(a,b):return [add(mul(a[(j+1)%3],b[(j+2)%3]),mul(a[(j+2)%3],b[(j+1)%3]),F(-1)) for j in range(3)]
    normal=cross(tangent[0],tangent[1]);Cxi=[{key:value/2 for key,value in p.items()} for p in cross(tangent[1],director)];Ceta=[{key:value/2 for key,value in p.items()} for p in cross(director,tangent[0])]
    pairs=[(i,j) for i in range(3) for j in range(i,3)];coupled={};nm=len(membrane_vectors)
    for col,vector in enumerate(membrane_vectors):
        u=[{key:vector[3*i+cart] for i,key in enumerate(MID) if vector[3*i+cart]} for cart in range(3)]
        B=[[add(mul(deriv(u[a],0),Cxi[b]),mul(deriv(u[a],1),Ceta[b])) for b in range(3)] for a in range(3)]
        for eq,(a,b) in enumerate(pairs):
            P=B[a][b] if a==b else add(B[a][b],B[b][a])
            for key,value in P.items():coupled.setdefault((eq,key),{})[col]=value
    for i,key in enumerate(MID):
        mon={key:F(1)}
        for cart in range(3):
            for eq,(a,b) in enumerate(pairs):
                P={}
                if cart==a:P=add(P,mul(mon,normal[b]))
                if a!=b and cart==b:P=add(P,mul(mon,normal[a]))
                for power,value in P.items():coupled.setdefault((eq,power),{})[nm+3*i+cart]=value
    joint,joint_pivots=exact_null(list(coupled.values()),nm+105)
    joint_zero=all(sum(value*vector[col] for col,value in row.items())==0 for row in coupled.values() for vector in joint)
    vectors=[[sum(v[i]*membrane_vectors[i][j] for i in range(nm)) for j in range(105)] for v in joint]
    coeff=np.zeros((len(thin.c.IDS),3,len(vectors)));acoef=np.zeros_like(coeff)
    for j,vector in enumerate(vectors):
        for i,(a,b) in enumerate(MID):
            for k,ka in thin.c.powers(a):
                for l,lb in thin.c.powers(b):
                    row=thin.c.LOOK[k,l,0]
                    for cart in range(3):
                        coeff[row,cart,j]+=float(vector[3*i+cart]*ka*lb)/thin.c.NORM[row]
                        acoef[row,cart,j]+=float(joint[j][nm+3*i+cart]*ka*lb)/thin.c.NORM[row]
    points=np.array([[s,t] for t,s in itertools.product(thin.m.LOBATTO,repeat=2)]);xx=sum(float(v)*points[:,0]**i*points[:,1]**j for (i,j),v in x.items());yy=sum(float(v)*points[:,0]**i*points[:,1]**j for (i,j),v in y.items())
    expected=np.vstack((np.column_stack((xx,yy,kappa*xx*xx/2)),np.column_stack((-kappa*xx,np.zeros(25),np.ones(25))),np.zeros((25,3))))
    encoded='\n'.join(','.join(str(x) for x in vector) for vector in vectors).encode()
    return coeff,acoef,dict(exact_membrane_rank=len(membrane_pivots),exact_membrane_nullity=nm,exact_director_columns=nm+105,exact_director_rank=len(joint_pivots),compatible_q0_nullity=len(vectors),pivot_columns=joint_pivots,rational_nullspace_sha256=hashlib.sha256(encoded).hexdigest(),all_rational_substitutions_zero=residual and joint_zero,all_three_components=True,
        scalar_polynomial_space_dimension=35,geometry='analytic_parabolic_cylinder_of_binary_rational_xy_corners',max_rounded_reference_difference=float(np.max(np.abs(expected-q))))

def local_basis(op,kappa):
    s=op['state'];ns=s['nscalar'];nc=3*ns;q=s['q'];poly,acoef,cert=certificate(np.asarray(q,dtype=np.float64).tobytes(),kappa);n=poly.shape[-1]
    if s['S'].shape[-1]!=280 or [x['n'] for x in s['faces']]!=[25,25,15,15,15,15]:raise ValueError('director_image_proof_requires_basic_RT_benchmark')
    p1=np.zeros_like(acoef)
    for row,(i,j,k) in enumerate(thin.c.IDS):
        if k==0:p1[thin.c.LOOK[i,j,1]]=acoef[row]/np.sqrt(3.)
    P0=np.zeros((nc,n));P0[:225]=poly[thin.f.VI].reshape(225,n)
    for face in s['faces']:
        cf=thin.c.restrict(poly,face['axis'],face['sign']);v=np.einsum('mn,maj->naj',face['Mcoeff'],cf,optimize=True);P0[3*face['offset']:3*(face['offset']+face['n'])]=v.reshape(3*face['n'],n)
    particular=np.zeros((nc,n));particular[:225]=p1[thin.f.VI].reshape(225,n)
    for face in s['faces']:
        cf=thin.c.restrict(p1,face['axis'],face['sign']);v=np.einsum('mn,maj->naj',face['Mcoeff'],cf,optimize=True);particular[3*face['offset']:3*(face['offset']+face['n'])]=v.reshape(3*face['n'],n)
    J=np.vstack([face['hat'].T@((face['w']/4)[:,None]*face['jump']) for face in s['faces']]);ZT,info=thin.nullspace(np.vstack((J[:sum(x['n'] for x in s['faces'][:2])],s['Gwref'][:,2,:])));free=np.kron(ZT,np.eye(3))
    Z=np.zeros((2*nc,n+free.shape[1]));Z[:nc,:n]=P0;Z[nc:,:n]=particular;Z[nc:,n:]=free;Z=np.linalg.qr(Z,mode='reduced')[0]
    C=thin.constraints(op);residual=float(np.linalg.norm(C@Z)/max(1.,np.linalg.norm(C)))
    if residual>1e-9:raise ValueError('analytic_isometry_lift_residual:'+str(residual))
    return Z,dict(certificate=cert,q1_homogeneous_kernel=info,nullity=Z.shape[1],residual=residual,geometry_conditional=True)
