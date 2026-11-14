"""Exact arithmetic comparison in the specified y-independent bending sector."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib,json
import mortar_q4 as numerical
ROOT=rd_paths.ROOT
K=F(2,5)
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def add(a,b):
    result=a.copy()
    for k,v in b.items():result[k]=result.get(k,F(0))+v
    return {k:v for k,v in result.items() if v}
def mul(a,b):
    result={}
    for i,x in a.items():
        for j,y in b.items():result[i+j]=result.get(i+j,F(0))+x*y
    return {k:v for k,v in result.items() if v}
def scale(a,c):return {k:v*c for k,v in a.items() if v*c}
def compose(poly,a,b):
    result={}
    for n,c in poly.items():
        result=add(result,{j:c*comb(n,j)*a**(n-j)*b**j for j in range(n+1)})
    return result
def integral(poly):return sum(2*v/F(k+1) for k,v in poly.items() if k%2==0)
def rank(rows):
    A=[row.copy() for row in rows];r=0
    for col in range(20):
        pivot=next((i for i in range(r,len(A)) if A[i][col]),None)
        if pivot is None:continue
        A[r],A[pivot]=A[pivot],A[r];v=A[r][col];A[r]=[x/v for x in A[r]]
        for i in range(r+1,len(A)):
            if A[i][col]:
                v=A[i][col];A[i]=[x-v*y for x,y in zip(A[i],A[r])]
        r+=1
        if r==len(A):break
    return r

FIELDS=[{0:{0:F(1)},2:{2:-K/2},3:{3:-2*K/3},4:{4:-3*K/4}},
        {j+1:{j:F(1)} for j in range(4)},
        {j+1:{j-1:-F(j)} for j in range(1,4)},
        {j+1:{j:-K*j} for j in range(1,4)}]

def build(amplitude,degree=None):
    mesh=numerical.broken_mesh(.4,2,amplitude);rows=[];edges=[]
    for left,right,axis in mesh['faces']:
        local=[4,24] if axis==0 else [20,24];nodes=mesh['cells'][left][local]
        xa,xb=[F.from_float(float(mesh['ref'][node,0])) for node in nodes];ya,yb=[F.from_float(float(mesh['ref'][node,1])) for node in nodes]
        a=(xa+xb)/2;b=(xb-xa)/2;d=(yb-ya)/2;x={0:a,1:b};sign=F(1 if axis==0 else -1)
        area=[{0:sign*d},scale(add({0:F(1)},scale(mul(x,x),K*K)),-sign*b),scale(x,sign*K*d)]
        edges.append(dict(left=left,right=right,axis=axis,x_center=str(a),x_slope=str(b),y_slope=str(d)))
        def row_for(field,weight=None,power=None):
            row=[F(0)]*20
            for local,poly in field.items():
                value=compose(poly,a,b)
                value=value.get(power,F(0)) if power is not None else integral(mul(value,weight))
                row[5*left+local]+=value;row[5*right+local]-=value
            if any(row):rows.append(row)
        if degree is None:
            for field in FIELDS:
                for power in range(5):row_for(field,power=power)
        else:
            for field in FIELDS:
                for power in range(degree+1):row_for(field,{power:F(1)})
            for field in FIELDS[:2]:
                for a0 in area:row_for(field,a0)
    for cell,ids in enumerate(mesh['cells']):
        if sum(mesh['ref'][ids[:25],0]==0.)!=5:continue
        for local in [0,1,2]:
            row=[F(0)]*20;row[5*cell+local]=F(1);rows.append(row)
    return rows,edges

def main():
    dest=ROOT/'rd/structure/main-shell/mortar/records/mortar-q4-sector-audit.json';assert not dest.exists();records=[]
    for amp in [0.,.2]:
        strong,_=build(amp);rs=rank(strong)
        for degree in [0,1]:
            rows,edges=build(amp,degree);r=rank(rows);union=rank(rows+strong)
            records.append(dict(amplitude=amp,degree=degree,unknowns=20,rows=len(rows),exact_rank=r,sector_dimension=20-r,
                strong_rank=rs,union_rank=union,same_sector_as_strong=bool(r==rs==union),edges=edges))
    sources=[Path(__file__),Path(numerical.__file__),ROOT/'rd/structure/main-shell/mortar/records/mortar-q4-protocol.md']
    dest.write_text(json.dumps(dict(records=records,arithmetic='fractions.Fraction',scope='ideal_bilinear_geometry_y_independent_cubic_W_quartic_U_c_zero_sector',
        full_discrete_kernel_exactly_certified=False,finite_h_lower_bound=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources}),indent=2),encoding='utf-8')
    print([(r['amplitude'],r['degree'],r['sector_dimension'],r['same_sector_as_strong']) for r in records])
if __name__=='__main__':main()
