"""Position-gradient, side-selected normal and panel-normal differential audit."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
import moment_position_recovery as r
import moment_recovery_normals as normals
import check_moment_position_recovery as primary
f=r.f;ROOT=rd_paths.ROOT
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))
def main():
    out=ROOT/'rd/structure/main-shell/moment-enriched/records/moment-recovery-geometry-audit.json';assert not out.exists();start=time.monotonic()
    grid=f.m.mesh(.4,2,.2);states=[f.prepare(grid['ref'][ids].reshape(-1),.001) for ids in grid['cells']];R=r.Recovery(states,2)
    X=np.zeros((R.size,3))
    for cell,s in enumerate(states):X[np.argmax(R.native[cell][:s['ni']],axis=1)]=s['X0'][:s['ni']]
    for item in R.facet.values():
        face=item['face'];s=states[item['cell']];o=item['offset'];X[o:o+face['n']]=s['X0'][face['offset']:face['offset']+face['n']]
    rng=np.random.default_rng(71104);v=1e-5*rng.normal(size=X.shape);current=X+1e-7*rng.normal(size=X.shape);records=[]
    p=np.array([[-.7,-.6,1.],[.4,-.1,1.],[.2,.7,1.]])
    for cell,s in enumerate(states):
        G=R.gradient(cell,p);physical=f.geometry(s['q'],s['h'],p,True)[1].transpose(0,2,1)
        gradient=rel(G@X-physical,physical)
        tang=G@current;dt=G@v;n,dn=normals.normal(tang[:,0,:],tang[:,1,:],dt[:,0,:],dt[:,1,:]);eps=1e-3
        plus=G@(current+eps*v);minus=G@(current-eps*v);fd=(normals.normal(plus[:,0,:],plus[:,1,:])-normals.normal(minus[:,0,:],minus[:,1,:]))/(2*eps)
        normal_error=rel(dn-fd,dn)
        pp=np.array([[-.7,-.6,1.],[.7,-.6,1.],[.7,.6,1.],[-.7,.6,1.]])
        A=R.evaluate(cell,pp);verts=A@current;dv=A@v;pn,pdn=normals.panel(verts,dv);pfd=(normals.panel(verts+eps*dv)-normals.panel(verts-eps*dv))/(2*eps);panel_error=rel(pdn-pfd,pdn)
        reversed_normal=rel(normals.panel(verts[[0,3,2,1]])+pn,pn)
        records.append(dict(cell=cell,reference_gradient=gradient,surface_normal_differential=normal_error,panel_normal_differential=panel_error,reversed_connectivity=reversed_normal))
    rejected=False
    try:normals.panel(np.array([[0.,0,0],[1.,0,0],[2.,0,0],[3.,0,0]]))
    except ValueError:rejected=True
    native=np.vstack(R.native);sv=np.linalg.svd(native,compute_uv=False);injective=sv.min()>1e-9*sv.max()
    # Internal volume/TB coordinates cannot influence another cell's position.
    local_leakage=0.
    for cell,s in enumerate(states):
        cols=np.argmax(R.native[cell][:s['ni']],axis=1)
        for other in range(4):
            if cell!=other:local_leakage=max(local_leakage,float(np.linalg.norm(R.evaluate(other,p)[:,cols])))
    passed=max([x[k] for x in records for k in x if k!='cell']+[local_leakage])<=1e-8 and rejected and injective
    paths=[Path(__file__),Path(r.__file__),Path(normals.__file__),ROOT/'rd/structure/main-shell/moment-enriched/records/moment-position-recovery-result.json']
    out.write_text(json.dumps(dict(records=records,degenerate_rejected=rejected,native_map_injective=bool(injective),native_singular_bounds=[float(sv.min()),float(sv.max())],internal_neighbor_leakage=local_leakage,passed=bool(passed),physical_adoption=False,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8')
    print('geometry',passed,records,'locality',local_leakage,flush=True)
if __name__=='__main__':main()
