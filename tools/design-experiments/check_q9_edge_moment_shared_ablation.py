"""Ablation of work bubble following the failed finite-h shared-space gate."""
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from scipy.linalg import lstsq
import check_q9_edge_moment_shared_space_v2 as screen
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    target=F/'q9-edge-moment-shared-ablation-result.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-edge-moment-shared-space-v2-result.json';r=json.loads(parent.read_text());assert len(r['cases'])==45 and not r['diagnostic_checks_pass']
    assert r['source_sha256']==digest(screen.__file__)
    for name,sha in r['dependencies'].items():assert sha==digest(ROOT/'tools/design-experiments'/name)
    original=screen.method.evaluate;start=time.monotonic();cases=[]
    def without_bubble(state,points):
        J,E,R=screen.method.r1(state,points)
        return J,screen.method.mandel(E),screen.method.mandel(R)
    try:
        for h,n in itertools.product([.1,.01,.001],[2,4]):
            record=dict(h=h,n=n)
            for stage in ['R1','R']:
                screen.method.evaluate=without_bubble if stage=='R1' else original
                a=screen.assemble('TRAP',n,h);D,W=a['D'],a['W'];U,s,_=np.linalg.svd(D,full_matrices=False);sc=np.linalg.svd(np.column_stack((D,W)),compute_uv=False);scale=sc[0];rd=int(np.count_nonzero(s>1e-10*scale));rc=int(np.count_nonzero(sc>1e-10*scale));sz=np.linalg.svd(W-U[:,:rd]@(U[:,:rd].T@W),compute_uv=False)
                xy=a['xy'][a['interior']];w=n*np.sin(np.pi*(xy[:,0]+.5))*np.sin(np.pi*(xy[:,1]+.5));d=lstsq(D,-W@w,cond=1e-10)[0]
                record[stage]=dict(w_dofs=W.shape[1],compatible_w=W.shape[1]-(rc-rd),projected_w=W.shape[1]-int(np.count_nonzero(sz>1e-10*scale)),best_fit_norm=float(np.linalg.norm(D@d+W@w)))
            cases.append(record);print(json.dumps(record),flush=True)
            assert time.monotonic()-start<300 and screen.memory.memory_bytes()<1024**3
    finally:screen.method.evaluate=original
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),protocol_sha256=digest(F/'q9-edge-moment-shared-ablation-protocol.md'),dependencies={Path(m.__file__).name:digest(m.__file__) for m in [screen,screen.method]},cases=cases,seconds=time.monotonic()-start,peak_bytes=screen.memory.memory_bytes(),completed=len(cases)==6,product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
