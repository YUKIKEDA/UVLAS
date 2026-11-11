"""Exact rational L2 orthogonalization before conversion to float64.

Only polynomial nullspace preprocessing is exact. FE operators/solves stay float64.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from fractions import Fraction as F
import hashlib,math
import numpy as np

def orthogonal_modal(q0vectors,joint,nm,MID,DIR,c):
    vectors=[]
    for j,u in enumerate(q0vectors):
        v={}
        for field,ids,source in [(0,MID,u),(1,DIR,joint[j][nm:])]:
            for i,(a,b) in enumerate(ids):
                for k,ka in c.powers(a):
                    for l,lb in c.powers(b):
                        for cart in range(3):
                            value=source[3*i+cart]*ka*lb
                            if value:key=(field,cart,k,l);v[key]=v.get(key,F(0))+value
        vectors.append({k:v for k,v in v.items() if v})
    def dot(a,b):return sum((v*b.get(k,F(0))/((2*k[2]+1)*(2*k[3]+1)) for k,v in a.items()),F(0))
    orth=[];squares=[]
    for v in vectors:
        v=v.copy()
        for q,norm in zip(orth,squares):
            ratio=dot(v,q)/norm
            if ratio:
                for key,value in q.items():
                    new=v.get(key,F(0))-ratio*value
                    if new:v[key]=new
                    else:v.pop(key,None)
        square=dot(v,v)
        if square<=0:raise ValueError('exact_nullspace_basis_lost_independence')
        orth.append(v);squares.append(square)
    exact_orthogonal=all(dot(a,b)==0 for j,a in enumerate(orth) for b in orth[:j])
    if not exact_orthogonal:raise ValueError('rational_orthogonalization_failed')
    out=[np.zeros((len(c.IDS),3,len(orth))) for _ in range(2)]
    for j,(v,norm) in enumerate(zip(orth,squares)):
        for (field,cart,k,l),value in v.items():
            # Ratio is bounded by the inverse polynomial L2 weight. Avoid
            # converting a potentially enormous Fraction numerator by itself.
            normalized=math.copysign(math.sqrt(float(value*value/norm)),1 if value>0 else -1)
            row=c.LOOK[k,l,0];out[field][row,cart,j]=normalized/c.NORM[row]
    encoded='\n'.join(';'.join(str(k)+':'+str(v) for k,v in sorted(q.items())) for q in orth).encode()
    return *out,dict(method='exact_rational_L2_orthogonalization_of_complete_joint_nullspace',columns=len(orth),all_exact_orthogonality_checks=exact_orthogonal,basis_sha256=hashlib.sha256(encoded).hexdigest())
