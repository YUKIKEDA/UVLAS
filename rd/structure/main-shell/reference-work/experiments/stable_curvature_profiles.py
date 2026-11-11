"""Decimal reference-only basis conversion; float64 Chebyshev FE evaluation."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from decimal import Decimal as D,localcontext
from pathlib import Path
import hashlib,json,time
import numpy as np
ROOT=rd_paths.ROOT
PI=D('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148')
def cosine(x):
    x=x%(2*PI)
    if x>PI:x-=2*PI
    term=D(1);total=term
    for j in range(1,101):term*=-(x*x)/D((2*j-1)*(2*j));total+=term
    return total
def atan(x):
    term=x;total=x
    for j in range(1,180):term*=-(x*x);total+=term/D(2*j+1)
    return total
def integrals(x,k):
    if not k:return [x**(j+1)/D(j+1) for j in range(10)]
    k2=k*k;g=1+k2*x*x;h=[atan(k*x)/k,g.ln()/(2*k2)]
    for j in range(2,8):h.append(x**(j-1)/(k2*D(j-1))-h[j-2]/k2)
    out=[x/(2*g)+atan(k*x)/(2*k),-1/(2*k2*g)]
    for j in range(2,10):out.append((h[j-2]-out[j-2])/k2)
    return out
def polys(left,width):
    b=-1-2*left/width;a=2/width;P=[[D(1)],[b,a]]
    for j in range(1,7):
        v=[D(0)]*(j+2)
        for i,c in enumerate(P[-1]):v[i]+=D(2*j+1)*b*c;v[i+1]+=D(2*j+1)*a*c
        for i,c in enumerate(P[-2]):v[i]-=D(j)*c
        P.append([c/D(j+1) for c in v])
    return P
def build(n,kappa,kind,degree=40):
    tick=time.monotonic();segments=2*n;nb=8*segments;columns=2*nb;rows=segments*(degree+1)
    with localcontext() as ctx:
        ctx.prec=100;k=D(str(kappa));width=D(1)/segments
        nodes=[cosine(PI*(D(j)+D('.5'))/D(degree+1)) for j in range(degree+1)]
        transform=[[cosine(PI*D(i)*(D(j)+D('.5'))/D(degree+1))*D(2)/D(degree+1)/(2 if i==0 else 1) for j in range(degree+1)] for i in range(degree+1)]
        pol=[polys(D(cell)*width,width) for cell in range(segments)];ends=[];starts=[]
        for cell in range(segments):
            a=integrals(D(cell)*width,k);b=integrals(D(cell+1)*width,k)
            starts.append([[sum(c*a[j+power] for j,c in enumerate(P)) for power in range(3)] for P in pol[cell]])
            ends.append([[sum(c*(b[j+power]-a[j+power]) for j,c in enumerate(P)) for power in range(3)] for P in pol[cell]])
        coeff=[[D(0)]*rows for _ in range(columns)]
        for cell in range(segments):
            vals=[[D(0)]*(degree+1) for _ in range(columns)]
            for pt,node in enumerate(nodes):
                x=(D(cell)+(node+1)/2)*width;J=integrals(x,k)
                for origin in range(cell+1):
                    for mode,P in enumerate(pol[origin]):
                        I=ends[origin][mode] if origin<cell else [sum(c*J[j+power] for j,c in enumerate(P))-starts[cell][mode][power] for power in range(3)]
                        index=8*origin+mode;i0,i1,i2=I
                        vals[index][pt]=x*i0-i1 if kind==0 else i0
                        vals[nb+index][pt]=(x*x*i0-i2)/2 if kind==0 else x*i0
            for col in range(columns):
                for i,row in enumerate(transform):coeff[col][cell*(degree+1)+i]=sum(a*b for a,b in zip(row,vals[col]))
        basis=[];remaining=coeff;norms=[sum(v*v for v in col) for col in remaining];largest=max(norms);smallest=None;discard=None
        while remaining:
            pivot=max(range(len(norms)),key=norms.__getitem__);norm=norms[pivot]
            if norm<=largest*D('1e-100'):discard=norm.sqrt();break
            length=norm.sqrt();q=[v/length for v in remaining.pop(pivot)];norms.pop(pivot);basis.append(q);smallest=length
            for j,col in enumerate(remaining):
                inner=sum(a*b for a,b in zip(q,col));remaining[j]=[a-inner*b for a,b in zip(col,q)];norms[j]=sum(v*v for v in remaining[j])
        encoded=np.array([[float(v) for v in q] for q in basis]).T.reshape(segments,degree+1,len(basis))
        return encoded,dict(n=n,kappa=kappa,kind=kind,degree=degree,decimal_digits=100,rank=len(basis),raw_columns=columns,last_kept=str(smallest),first_discarded=str(discard),seconds=time.monotonic()-tick)
_cache={}
def get(n,k,kind):
    key=(n,float(k),kind)
    if key not in _cache:
        directory=ROOT/'build/analysis/stable-curvature-profiles';directory.mkdir(exist_ok=True);stem=f'n{n}-k{k:g}-f{kind}-d40-p100';path=directory/(stem+'.npz');meta=directory/(stem+'.json')
        if path.exists():
            with np.load(path) as data:coef=data['coef'].copy()
        else:
            coef,record=build(n,k,kind);np.savez_compressed(path,coef=coef);record['artifact_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();meta.write_text(json.dumps(record,indent=2),encoding='utf-8');print('STABLE_PROFILE',record,flush=True)
        _cache[key]=coef
    return _cache[key]
def profiles(x,n,axis,kind,metric):
    coef=get(n,metric[axis],kind);lo=0. if axis==0 else -.5;u=np.clip(np.asarray(x)-lo,0.,1.);cell=np.minimum((u*(2*n)).astype(int),2*n-1);local=2*(u*(2*n)-cell)-1;value=np.empty((len(u),coef.shape[-1]));deriv=np.empty_like(value)
    for j in np.unique(cell):
        mask=cell==j;value[mask]=np.polynomial.chebyshev.chebvander(local[mask],coef.shape[1]-1)@coef[j]
        dc=np.polynomial.chebyshev.chebder(coef[j],axis=0);deriv[mask]=np.polynomial.chebyshev.chebvander(local[mask],dc.shape[0]-1)@dc*(4*n)
    return value,deriv
