"""Research-only Decimal coefficient preparation of the unchanged GSP space."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from decimal import Decimal,localcontext
import itertools,time
import numpy as np
from scipy.linalg import cholesky,solve_triangular
import gsp_q3 as old

D=lambda x:Decimal.from_float(float(x))
ONE={(0,0,0):Decimal(1)}
def add(a,b):
    c=a.copy()
    for k,v in b.items():c[k]=c.get(k,Decimal(0))+v
    return {k:v for k,v in c.items() if v}
def scale(a,s):return {k:v*s for k,v in a.items() if v*s}
def mul(a,b):
    c={}
    for i,x in a.items():
        for j,y in b.items():
            k=tuple(a+b for a,b in zip(i,j));c[k]=c.get(k,Decimal(0))+x*y
    return {k:v for k,v in c.items() if v}
def derivative(a,axis):
    out={}
    for k,v in a.items():
        if k[axis]:
            new=list(k);new[axis]-=1;out[tuple(new)]=v*k[axis]
    return out
def total(items):
    out={}
    for item in items:out=add(out,item)
    return out
def shift(a,axis):return {tuple(n+(i==axis) for i,n in enumerate(k)):v for k,v in a.items()}

def raw_polynomials(state):
    # Preserve every stored float64 Vandermonde coefficient of the primary
    # geometry, interpreted exactly; no refit or geometry rounding.
    cardinal=[{(power,0,0):D(old.BASE.geometry.COEFFICIENT[power,i]) for power in range(4)
        if old.BASE.geometry.COEFFICIENT[power,i]} for i in range(4)]
    N=[]
    for t,s in itertools.product(range(4),range(4)):
        pt={(0,k[0],0):v for k,v in cardinal[t].items()};N.append(mul(cardinal[s],pt))
    q=state['q'].reshape(16,3,3)
    fields=[[total(scale(N[n],D(q[n,k,a])) for n in range(16)) for a in range(3)] for k in range(3)]
    A=[[scale(derivative(fields[0][a],b) if b<2 else fields[1][a],1/D(state['frame_scale'][b])) for b in range(3)] for a in range(3)]
    det={}
    for perm in itertools.permutations(range(3)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        det=add(det,scale(mul(mul(A[0][perm[0]],A[1][perm[1]]),A[2][perm[2]]),Decimal(sign)))
    mapped=[]
    for k in range(6):
        components=[]
        for out in range(6):
            parts=[]
            for a,b,i,j in itertools.product(range(3),repeat=4):
                constant=D(old.MB[out,a,b])*D(old.MB[k,i,j])
                if constant:parts.append(scale(mul(A[a][i],A[b][j]),constant))
            components.append(total(parts))
        mapped.append(components)
    leg=[]
    rawlegs=[[Decimal(1)],[Decimal(0),Decimal(1)],[-Decimal(1)/2,Decimal(0),Decimal(3)/2],
        [Decimal(0),-Decimal(3)/2,Decimal(0),Decimal(5)/2]]
    for axis in range(3):
        axisleg=[]
        for degree,coeff in enumerate(rawlegs):
            p={}
            for power,value in enumerate(coeff):
                key=[0,0,0];key[axis]=power
                if value:p[tuple(key)]=value*D(np.sqrt((2*degree+1)/2))
            axisleg.append(p)
        leg.append(axisleg)
    columns=[]
    for component,imax,jmax,kmax in [(0,2,3,2),(1,3,2,2),(2,3,3,1),(3,3,2,2),(4,2,3,2),(5,2,2,2)]:
        for i,j,k in itertools.product(range(imax+1),range(jmax+1),range(kmax+1)):
            p=mul(mul(leg[0][i],leg[1][j]),leg[2][k]);columns.append([mul(x,p) for x in mapped[component]])
    phi=[ONE.copy()]+[scale(add(fields[0][a],scale(ONE,-D(state['center'].reshape(3)[a]))),1/D(state['length'])) for a in range(3)]
    phi += [shift(fields[1][a],2) for a in range(3)]
    phi += [scale(shift(shift(fields[2][a],2),2),D(state['length'])) for a in range(3)]
    for p in phi:
        numerator=mul(det,p)
        for out in range(6):columns.append([numerator if k==out else {} for k in range(6)])
    assert len(columns)==263
    return columns,det

def powers_to_legendre(n):
    table=[[Decimal(1)]]
    for degree in range(1,n):
        new=[Decimal(0)]*(degree+1)
        for k,v in enumerate(table[-1]):
            new[k+1]+=v*Decimal(k+1)/Decimal(2*k+1)
            if k:new[k-1]+=v*Decimal(k)/Decimal(2*k+1)
        table.append(new)
    return table

def prepare_coefficients(state,digits):
    start=time.monotonic()
    with localcontext() as ctx:
        ctx.prec=digits
        columns,det=raw_polynomials(state)
        selected=np.flatnonzero(np.linalg.norm(state['T'],axis=1)>0)
        rows={}
        for col,index in enumerate(selected):
            for component,p in enumerate(columns[index]):
                for powers,value in p.items():
                    assert powers[0]<=11 and powers[1]<=11 and powers[2]<=2
                    key=(*powers,component)
                    if key not in rows:rows[key]=[Decimal(0)]*len(selected)
                    rows[key][col]=value
        keys=list(rows)
        transforms=np.array([[D(x) for x in row] for row in state['T'][selected]],dtype=object)
        combined=np.array([rows[k] for k in keys],dtype=object)@transforms
        # Convert each polynomial axis before rounding coefficients to float64.
        coeff=np.full((12,12,3,6,len(selected)),Decimal(0),dtype=object)
        for key,row in zip(keys,combined):coeff[key]=row
        tables=[powers_to_legendre(12),powers_to_legendre(12),powers_to_legendre(3)]
        for axis,table in enumerate(tables):
            matrix=np.full((len(table),len(table)),Decimal(0),dtype=object)
            for power,row in enumerate(table):matrix[power,:len(row)]=row
            coeff=np.moveaxis(coeff,axis,-1)@matrix
            coeff=np.moveaxis(coeff,-1,axis)
        result=np.array(coeff,dtype=float)
        denominator=np.zeros((9,9))
        for (i,j,k),v in det.items():assert k==0;denominator[i,j]=float(v)
    return dict(coeff=result,denominator=denominator),dict(digits=digits,seconds=time.monotonic()-start,
        coefficient_bytes=result.nbytes+denominator.nbytes,coefficient_max=float(np.max(abs(result))))

def evaluate(poly,points):
    v=[np.polynomial.legendre.legvander(points[:,i],n) for i,n in enumerate([11,11,2])]
    phi=np.einsum('pi,pj,pk->pijk',*v).reshape(len(points),-1)
    numerator=(phi@poly['coeff'].reshape(432,-1)).reshape(len(points),6,-1)
    den=np.polynomial.polynomial.polyval2d(points[:,0],points[:,1],poly['denominator'])
    return numerator/den[:,None,None]

def whiten(state,poly):
    p,w=old.BASE.quadrature(9);mu,_=old.measure(state,0.,p,w);S=evaluate(poly,p)
    A=(np.einsum('ij,pjn->pin',old.CIH,S)*np.sqrt(mu)[:,None,None]).reshape(-1,S.shape[2])
    H=A.T@A;L=cholesky(H,lower=True,check_finite=False)
    # Well-conditioned correction of the evaluated polynomial basis, no new space.
    coeff=solve_triangular(L,poly['coeff'].reshape(-1,len(L)).T,lower=True,check_finite=False).T.reshape(poly['coeff'].shape)
    return dict(coeff=coeff,denominator=poly['denominator'])
