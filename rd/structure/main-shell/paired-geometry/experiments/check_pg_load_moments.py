"""80-digit direct body+six-face thickness integrals, independent of series moments."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from decimal import Decimal as D,localcontext
from pathlib import Path
import hashlib,json
import numpy as np
import pg_continuum_loads as loads
from developable_general_reference import decimal_gauss
ROOT=rd_paths.ROOT

def fields(x,z,k,mode,target):
    nu=D('.3');mu=1/(2*(1+nu));lam=nu/((1+nu)*(1-2*nu));g=1+k*k*x*x;j=g-k*z
    if mode=='rational':a=D([2,4,6][target]);v=1/(1+a*a*x*x);vp=-2*a*a*x/(1+a*a*x*x)**2
    elif target==0:v=x*x;vp=2*x
    elif target==1:v=20*x**3;vp=60*x*x
    else:v=x.exp();vp=v
    f=-z*v/j;fx=-z*(vp*j-v*2*k*k*x)/j**2;fz=-v*g/j**2
    T=[lam*g+2*mu,lam*g,lam*g+2*mu*k*k*x*x,2*mu*k*x];Tx=[2*lam*k*k*x,2*lam*k*k*x,2*lam*k*k*x+4*mu*k*k*x,2*mu*k]
    S=[f*t for t in T];Sx=[fx*t+f*dt for t,dt in zip(T,Tx)];Sz=[fz*t for t in T]
    b=[-(Sx[0]+k*x*Sx[3]-k*x*Sz[0]+(1-k*z)*Sz[3])/j,D(0),-(Sx[3]+k*x*Sx[2]-k*x*Sz[3]+(1-k*z)*Sz[2])/j]
    return S,b,j

def direct(x,k,h,mode):
    out=np.zeros((3,3,3));side_x=np.zeros_like(out);side_y=np.zeros_like(out)
    with localcontext() as ctx:
        ctx.prec=80;x=D(str(x));k=D(str(k));h=D(str(h));nodes,weights=decimal_gauss(20)
        for target in range(3):
            for l in range(3):
                volume=[D(0)]*3;surface=[D(0)]*3;xx=[D(0)]*3;yy=[D(0)]*3;fac=D(2 if l==2 else 1)
                for t,w in zip(nodes,weights):
                    z=h*t/2;weight=w*h/2;S,b,j=fields(x,z,k,mode,target);phi=z**l/fac
                    for i in range(3):volume[i]+=weight*phi*j*b[i]
                    tx=[S[0]+k*x*S[3],D(0),S[3]+k*x*S[2]];ty=[D(0),j*S[1],D(0)]
                    for i in range(3):xx[i]+=weight*phi*tx[i];yy[i]+=weight*phi*ty[i]
                for sign in [-1,1]:
                    z=sign*h/2;S,_,_=fields(x,z,k,mode,target);phi=z**l/fac;traction=[-k*x*S[0]+(1-k*z)*S[3],D(0),-k*x*S[3]+(1-k*z)*S[2]]
                    for i in range(3):surface[i]+=sign*phi*traction[i]
                out[l,:,target]=[float(a+b) for a,b in zip(volume,surface)];side_x[l,:,target]=[float(a) for a in xx];side_y[l,:,target]=[float(a) for a in yy]
    return out,side_x,side_y

def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-load-moment-audit.json';assert not dest.exists();records=[]
    for mode in ['standard','rational']:
        for h in [.01,.0001]:
            x=np.array([0.,.17,.43,.79,1.]);v,vp=loads.curvature(x,mode);actual=loads.moments(x,.4,h,v,vp);expected=[np.stack([direct(xx,.4,h,mode)[j] for xx in x],axis=1) for j in range(3)]
            differences=[float(np.linalg.norm(a-b)/np.linalg.norm(b)) for a,b in zip(actual,expected)];records.append(dict(mode=mode,h=h,relative_differences=differences,passed=max(differences)<=1e-10));print(records[-1],flush=True)
    paths=[Path(__file__),Path(loads.__file__),ROOT/'rd/structure/main-shell/q3/experiments/developable_general_reference.py',ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-physical-pde-protocol.md']
    result=dict(passed=all(r['passed'] for r in records),records=records,precision=80,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest.write_text(json.dumps(result,indent=2),encoding='utf-8')

if __name__=='__main__':main()
