"""Canonical polynomial quotient, avoiding geometry-dependent generator rank."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from decimal import Decimal,localcontext
from fractions import Fraction as F
from functools import lru_cache
import itertools
import numpy as np
import qem_el2_q4 as base
from gsp_q3_basis_v2 import D,ONE,add,mul,scale,derivative,total
SIZE=486
def index(c,i,j):return c*81+i*9+j
def base_rows():
    rows=[];pivots=[]
    def insert(row):
        row={k:F(v) for k,v in row.items() if v}
        for pivot,b in zip(pivots,rows):
            value=row.get(pivot,0)
            if value:
                for key,v in b.items():row[key]=row.get(key,0)-value*v
                row={k:v for k,v in row.items() if v}
        if row:
            pivot=min(row);value=row[pivot];rows.append({k:v/value for k,v in row.items()});pivots.append(pivot)
    for component in [0,1]:
        for i,j in itertools.product(range(5),repeat=2):
            row={}
            if component==0:
                if i:row[index(0,i-1,j)]=i
                if j:row[index(5,i,j-1)]=F(j,2)
            else:
                if j:row[index(1,i,j-1)]=j
                if i:row[index(5,i-1,j)]=F(i,2)
            insert(row)
    for i,j in itertools.product(range(5),repeat=2):
        row={}
        if i:row[index(4,i-1,j)]=F(i,2)
        if j:row[index(3,i,j-1)]=F(j,2)
        insert(row)
    for i,j in itertools.product(range(9),repeat=2):insert({index(2,i,j):1})
    assert len(rows)==152
    return rows,pivots
BASE_ROWS,BASE_PIVOTS=base_rows()

def geometry_polynomials(q):
    nodes=base.LOBATTO;cardinal=[]
    for i,x in enumerate(nodes):
        p=ONE.copy();den=Decimal(1)
        for j,y in enumerate(nodes):
            if i!=j:p=mul(p,{(1,0,0):Decimal(1),(0,0,0):-D(y)});den*=D(x)-D(y)
        cardinal.append(scale(p,1/den))
    N=[]
    for j,i in itertools.product(range(5),repeat=2):N.append(mul(cardinal[i],{(0,k[0],0):v for k,v in cardinal[j].items()}))
    coords=q.reshape(75,3)
    r=[total(scale(N[i],D(coords[i,a])) for i in range(25)) for a in range(3)]
    d=[total(scale(N[i],D(coords[25+i,a])) for i in range(25)) for a in range(3)]
    J=[[derivative(r[a],0),derivative(r[a],1),d[a]] for a in range(3)]
    return r,J

@lru_cache(maxsize=128)
def coefficients(qbytes,digits=80):
    q=np.frombuffer(qbytes,dtype=np.float64)
    with localcontext() as ctx:
        ctx.prec=digits;r,J=geometry_polynomials(q)
        def vector(fields):
            out=[Decimal(0)]*SIZE
            for c,p in enumerate(fields):
                for (i,j,k),value in p.items():assert k==0 and max(i,j)<=8;out[index(c,i,j)]=value
            return out
        generators=[]
        rcenter=[p.get((0,0,0),Decimal(0)) for p in r]
        rc=[add(p,scale(ONE,-center)) for p,center in zip(r,rcenter)]
        for axis in np.eye(3):
            omega=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
            v=[total(scale(rc[b],D(omega[a,b])) for b in range(3)) for a in range(3)]
            normal=total(mul(J[a][2],v[a]) for a in range(3))
            sz=scale(total(mul(J[a][0],v[a]) for a in range(3)),Decimal('.5'))
            tz=scale(total(mul(J[a][1],v[a]) for a in range(3)),Decimal('.5'))
            generators.append(vector([{}, {},normal,tz,sz,{}]))
        for mb in base.MB:
            fields=[]
            for i,j in [(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]:
                fields.append(total(scale(mul(J[a][i],J[b][j]),D(mb[a,b])) for a,b in itertools.product(range(3),repeat=2) if mb[a,b]))
            generators.append(vector(fields))
        complements=[];residual_ratios=[]
        for original in generators:
            v=original.copy();original_norm=sum(a*a for a in original).sqrt()
            if not original_norm:residual_ratios.append(0.);continue
            for pivot,row in zip(BASE_PIVOTS,BASE_ROWS):
                value=v[pivot]
                if value:
                    for i,frac in row.items():v[i]-=value*Decimal(frac.numerator)/Decimal(frac.denominator)
                    v[pivot]=Decimal(0)
            for _ in range(2):
                for b in complements:
                    dot=sum(a*c for a,c in zip(v,b));v=[a-dot*c for a,c in zip(v,b)]
            norm=sum(a*a for a in v).sqrt();residual_ratios.append(float(norm/original_norm))
            if norm>Decimal('1e-12')*original_norm:complements.append([a/norm for a in v])
        arrays=[]
        for row in BASE_ROWS:
            v=np.zeros(SIZE)
            for i,value in row.items():v[i]=float(value)
            arrays.append(v)
        arrays.extend(np.array([float(x) for x in b]) for b in complements)
        result=np.array(arrays).T.reshape(6,9,9,len(arrays))
        return result,dict(base_dimension=152,enrichment_dimension=len(complements),digits=digits,residual_ratios=residual_ratios)

def evaluate(q,points,J,independent=False):
    coeff,meta=coefficients(np.asarray(q,dtype=np.float64).tobytes(),120 if independent else 80)
    if independent:
        hat=np.array([np.polynomial.polynomial.polyval2d(points[:,0],points[:,1],coeff[c]) for c in range(6)]).transpose(2,0,1)
    else:
        v=np.einsum('pi,pj->pij',points[:,0,None]**np.arange(9),points[:,1,None]**np.arange(9)).reshape(len(points),81)
        hat=(v@coeff.transpose(1,2,0,3).reshape(81,-1)).reshape(len(points),6,-1)
    tensor=np.zeros((len(points),3,3,hat.shape[-1]));pairs=[(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]
    for c,(i,j) in enumerate(pairs):tensor[:,i,j]=hat[:,c];tensor[:,j,i]=hat[:,c]
    inv=np.linalg.inv(J);physical=np.einsum('pia,pijn,pjb->pabn',inv,tensor,inv,optimize=True)
    return np.einsum('kab,pabn->pkn',base.MB,physical),meta
