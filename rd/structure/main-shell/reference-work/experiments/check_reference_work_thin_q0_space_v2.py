"""Post-derived global q0-space rank; no pure-bending restriction in assembly."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.linalg import svd,qr
import check_reference_work_thin_global_v4 as solver
t=solver.t;ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def field(q,p,k):
    N,_=t.m.shape(p,0);x,y=(N@q.reshape(75,3)[:25])[:,:2].T
    return np.column_stack((-k*x**3/3,np.zeros(len(x)),x*x/2,-k*x**4/8,np.zeros(len(x)),x**3/6,-k*x**3*y/3,-k*x**4/24,x*x*y/2,-k*x*x/2,np.zeros(len(x)),x,-k*x*x*y/2,-k*x**3/6,x*y))
def main():
    dest=ROOT/'rd/structure/main-shell/reference-work/records/reference-work-thin-q0-space-v2-result.json';assert not dest.exists();records=[]
    for n in [2,4]:
        g=t.m.mesh(.4,n,.2);cells=[t.local(g['ref'][ids].reshape(-1),.4) for ids in g['cells']];Z,info,_,_=solver.assemble_constraints(cells,n)
        image=[];known=[]
        for cell in cells:
            ys=slice(cell['yoffset'],cell['yoffset']+cell['ny']);image.append(cell['Z'][:cell['nc']]@Z[ys]);s=t.operators(cell['q'])['state'];v=t.f.reduce(s,lambda p:field(cell['q'],p,.4))
            known.append(v.reshape(s['nscalar'],5,3).transpose(0,2,1).reshape(cell['nc'],5))
        image=np.vstack(image);known=np.vstack(known);U,sv,_=svd(image,full_matrices=False,check_finite=False);ranks=t.rank_info(sv);Q=qr(known,mode='economic',check_finite=False)[0]
        difference=float(np.linalg.norm(U[:,:ranks[0]]-Q@(Q.T@U[:,:ranks[0]])))
        checks=dict(q0_dimension=ranks==[5]*3,span=difference<=1e-8)
        records.append(dict(n=n,q0_ranks=ranks,q01_dimensions=info['nullity'],known_fields=['cylindrical quadratic','cylindrical cubic','quadratic y-dependent twisting','linear rotation-like field','linear y-dependent twisting'],span_difference=difference,
            all_three_positions_and_y_dependence_retained_before_rank=True,checks=checks))
        print(records[-1],flush=True)
    ps=[Path(__file__),Path(solver.__file__),ROOT/'rd/structure/main-shell/reference-work/records/reference-work-thin-global-v4-result.json'];dest.write_text(json.dumps(dict(records=records,geometry_conditional=True,all_checks_pass=all(all(x['checks'].values()) for x in records),sources={p.relative_to(ROOT).as_posix():digest(p) for p in ps}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
