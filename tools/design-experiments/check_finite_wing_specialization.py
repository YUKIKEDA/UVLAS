"""Specialization proof against the adopted complete W-P1 implementation."""
from pathlib import Path
import numpy as np,json,hashlib,time
import uvlm_wake_assembly_v2 as full
from uvlm_steady_reference import geometry,matrix,run
ROOT=Path(__file__).resolve().parents[2]
def binding(m,n,d):
    g=geometry(m,n,d);xi=(1-np.cos(np.linspace(0,np.pi,m+1)))/2;xv=np.r_[.75*xi[:-1]+.25*xi[1:],1.];xc=.25*xi[:-1]+.75*xi[1:];y=g['y']
    def weights(x,y):
        x,y=np.broadcast_arrays(x,y);x=x.ravel();y=y.ravel();N=np.zeros((len(x),6));Dx=N.copy();Dy=N.copy()
        for k,(xx,yy) in enumerate(zip(x,y)):
            cell=0 if yy<0 else 1;s=(yy+4)/4 if cell==0 else yy/4
            for node,wx,dx in [(0,1-xx,-1),(1,xx,1)]:
                N[k,2*cell+node]=(1-s)*wx;N[k,2*(cell+1)+node]=s*wx;Dx[k,2*cell+node]=(1-s)*dx;Dx[k,2*(cell+1)+node]=s*dx;Dy[k,2*cell+node]=-wx/4;Dy[k,2*(cell+1)+node]=wx/4
        return N,Dx,Dy
    Nv=weights(xv[:,None],y[None,:])[0];Nc,Dx,Dy=weights(xc[:,None],((y[:-1]+y[1:])/2)[None,:]);q=[]
    for yy in [-4.,0.,4.]:
        for xx in [0.,1.]:q.append([abs(yy)*np.tan(np.deg2rad(35))+xx*(1.5-.9*abs(yy)/4),yy,abs(yy)*np.tan(np.deg2rad(d))])
    return dict(sid='wing',m=m,n=n,Nv=Nv,Nc=Nc,Dx=Dx,Dy=Dy,q0=np.array(q)),g
def main():
    start=time.monotonic();records=[]
    for d in [0.,5.]:
        m,n=2,4;b,g=binding(m,n,d);model=full.Model([b],.001);geom=model.geometry(b['q0']);assert np.max(abs(geom['X']-g['X'].reshape(-1,3)))<1e-12;assert np.max(abs(geom['P']-g['P']))<1e-12
        te=geom['X'][model.surface_info[0]['te']];tail=te+[32.,0,0];ids0=tuple(model.vertex_ids[i] for i in model.surface_info[0]['te']);ids1=tuple(('tail',j) for j in range(n+1));w=full.Wake('wing',np.array([te,tail]),np.zeros((2,n)),(ids0,ids1),(('wing','attached'),('wing','released',1)),(0.,-1.6));A,W=full.Workspace().assemble(model,geom,(w,));closed=A+W[:,n:2*n]@model.T;fast=matrix(g,.001,32.);err=float(np.linalg.norm(closed-fast)/np.linalg.norm(closed));assert err<1e-10
        records.append(dict(case='full_WP1_constant_wake_specialization',dihedral=d,matrix_error=err))
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'tools/design-experiments/uvlm_steady_reference.py',ROOT/'tools/design-experiments/uvlm_wake_assembly_v2.py',ROOT/'docs/verification/uvlm-independent-reference-protocol.md']}
    result=dict(passed=True,records=records,failures=[],sources=sources,seconds=time.monotonic()-start)
    (ROOT/'docs/verification/finite-wing-specialization-result.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
