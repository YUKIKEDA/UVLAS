import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np,json,hashlib,time
import uvlm_wake_assembly_v2 as full
from check_finite_wing_specialization import binding
from uvlm_drag_force_kernel import matrix
ROOT=rd_paths.ROOT
def main():
    start=time.monotonic();rows=[]
    for eps in [0.,.02]:
        m,n=2,4;b,g=binding(m,n,0.);model=full.Model([b],eps);geom=model.geometry(b['q0']);te=geom['X'][model.surface_info[0]['te']];tail=te+[32.,0,0];gamma=np.random.default_rng(71).normal(size=m*n);mu=np.tile(gamma[-n:],(2,1));ids=tuple(model.vertex_ids[i] for i in model.surface_info[0]['te']);w=full.Wake('wing',np.array([te,tail]),mu,(ids,tuple(('tail',j) for j in range(n+1))),(('wing','attached'),('wing','released',1)),(0.,-1.6));X=g['X'].copy();X[-1,:,0]+=32.;p=(g['left']+g['right'])/2;dl=g['right']-g['left'];N=np.column_stack([np.zeros(len(dl)),-dl[:,2],dl[:,1]]);C=matrix(p,N,X,eps,np.arange(m*n));reference=[]
        for k in range(m*n):
            a,b=model.faces[k,:2];edge=model.edges.index(tuple(sorted((int(a),int(b)))));v,_=full.field(model,geom,(w,),gamma,p[k],edge_id=edge);reference.append(N[k]@v)
        reference=np.array(reference);err=float(np.linalg.norm(C@gamma-reference)/np.linalg.norm(reference));assert err<1e-8;rows.append(dict(eps=eps,relative_error=err))
    bad=p.copy();bad[0,0]+=.01
    try:matrix(bad,N,X,.02,np.arange(m*n))
    except ValueError as e:assert str(e)=='false_material_self_line'
    else:raise AssertionError('false self line accepted')
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'rd/aero/finite-wing/experiments/uvlm_drag_force_kernel.py',ROOT/'rd/aero/finite-wing/experiments/check_finite_wing_specialization.py',ROOT/'rd/aero/wake-loads/experiments/uvlm_wake_assembly_v2.py']};r=dict(passed=True,records=rows,false_self_rejected=True,sources=sources,seconds=time.monotonic()-start);dest=ROOT/'rd/aero/finite-wing/records/drag-material-self-result.json';assert not dest.exists();dest.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
