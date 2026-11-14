"""Reobserve saved fields and independently verify the new mass-root rejection."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.linalg import svd
import check_pg_ritz_finite_v3 as engine
import pg_ritz_shell_v3 as protected_shell
import pg_ritz_clamped_kernel as extended

ROOT=rd_paths.ROOT


def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-r4-followup-result.json';assert not dest.exists()
    sources=[];records=[]
    model=protected_shell.prepare(4,(.4,.2,0.));engine.shell=protected_shell
    for h in [.01,.0001]:
        path=ROOT/f'build/analysis/pg-ritz-protected-elliptic-n4-h{h}.npz';sources.append(path)
        a=np.load(path);Z=np.kron(a['Zs'],np.eye(3))
        sol=dict(coef=a['coef'],beta=a['beta'],Z=Z,scale=np.ones(Z.shape[1]))
        # Different observation quadrature; raw residual scale is not an acceptance claim.
        o=engine.observe(model,h,sol,order=17)
        r=dict(case='protected_reference_lift',parameters=[.4,.2,0.],n=4,h=h,
               error=o['error'],position=o['position'])
        records.append(r);print(r,flush=True)
    ratio=np.array(records[1]['position'])/records[0]['position']
    assert np.all(ratio>2.)
    model=extended.prepare(2,(.4,.2,0.));plane,w=extended.lift.geometry.quadrature(2,15)
    z,wz=np.polynomial.legendre.leggauss(9)
    p=np.column_stack((np.repeat(plane,9,axis=0),np.tile(z,len(plane))))
    weight=(w[:,None]*wz).ravel()
    for h in [.01,.0001]:
        N,G,B,det,X=extended.kinematics(model,p,h)
        A=N*np.sqrt(weight*det)[:,None];A/=np.linalg.norm(A,axis=0)
        ss=svd(A,compute_uv=False,check_finite=False)
        ranks=[int(sum(ss>tol*ss[0])) for tol in [1e-9,1e-10,1e-11]]
        boundary=p[:137].copy();boundary[:,0]=0
        trace=extended.kinematics(model,boundary,h)[0]
        defect=float(np.max(abs(trace)))
        assert defect<=1e-9
        records.append(dict(case='independent_full_trace_lift_mass_root',n=2,h=h,
            dofs=model['ns'],rank=ranks,rcond=float(ss[-1]/ss[0]),trace=defect))
        print(records[-1],flush=True)
    assert len(set(records[-1]['rank']))>1
    # Freeze the actual imported local module closure, not just the top-level script.
    import ast
    queue=[Path(__file__)];seen=set()
    while queue:
        path=queue.pop()
        if path in seen:continue
        seen.add(path);sources.append(path)
        tree=ast.parse(path.read_text(encoding='utf-8-sig'))
        for node in ast.walk(tree):
            names=[x.name.split('.')[0] for x in node.names] if isinstance(node,ast.Import) else [node.module.split('.')[0]] if isinstance(node,ast.ImportFrom) and node.module else []
            for name in names:
                local=path.parent/(name+'.py')
                if local.exists():queue.append(local)
    result=dict(rejection_reproduced=True,accepted=False,R4_complete=False,records=records,
        position_thin_ratio=ratio.tolist(),sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8')


if __name__=='__main__':main()
