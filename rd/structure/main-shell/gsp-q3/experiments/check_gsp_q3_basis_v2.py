"""Frozen two-state gate before any GSP local v2 continuation."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import qr,svdvals
import gsp_q3_basis_v2 as v2
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))

def main():
    dest=ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-basis-v2-result.json'
    assert not dest.exists()
    start=time.monotonic();out=ROOT/'build/analysis'/('gsp-q3-basis-v2-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    parentpath=ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-local-result.json';parent=json.loads(parentpath.read_text())
    records=[];stop=None
    try:
        for index,item in enumerate(parent['states']):
            source=ROOT/item['artifact'];assert digest(source)==item['artifact_sha256']
            data=np.load(source);state={k:data[k] for k in ['q','frame_scale','length','center','T','active','column_norms']}
            bases={};preparations=[]
            for digits in [80,120]:
                raw,prep=v2.prepare_coefficients(state,digits);poly=v2.whiten(state,raw)
                bases[digits]=(raw,poly);preparations.append(prep)
                print('prepared',index,prep,flush=True)
            rng=np.random.default_rng(20261103);points=rng.uniform(-.97,.97,(17,3))
            agreement=rel(v2.evaluate(bases[80][1],points)-v2.evaluate(bases[120][1],points),v2.evaluate(bases[120][1],points))
            gates=[]
            for order in [9,11,13]:
                p,w=v2.old.BASE.quadrature(order);mu,_=v2.old.measure(state,0.,p,w)
                S=v2.evaluate(bases[120][1],p)
                A=(np.einsum('ij,pjn->pin',v2.old.CIH,S)*np.sqrt(mu)[:,None,None]).reshape(-1,S.shape[2])
                raw=v2.old.stress_generator(state,p)
                G=(np.einsum('ij,pjn->pin',v2.old.CIH,raw)*np.sqrt(mu)[:,None,None]).reshape(-1,263)
                Q,_=qr(A,mode='economic',check_finite=False)
                rank=v2.old.ranks(svdvals(G[:,state['active']]/state['column_norms'][state['active']],check_finite=False))
                gates.append(dict(order=order,ranks=rank,span=rel(G-Q@(Q.T@G),G),
                    gram_identity=rel(A.T@A-np.eye(A.shape[1]),np.eye(A.shape[1])),
                    gram=(A.T@A)))
            path=out/f'state-{index}.npz'
            np.savez_compressed(path,**state,coeff80=bases[80][1]['coeff'],coeff120=bases[120][1]['coeff'],
                raw_coeff120=bases[120][0]['coeff'],denominator=bases[120][1]['denominator'],
                arbitrary_points=points,**{f'gram{r["order"]}':r['gram'] for r in gates})
            for gate in gates:gate.pop('gram')
            checks=dict(precision_agreement=agreement<=1e-10,
                original_gates=all(g['ranks']==[state['T'].shape[1]]*3 and max(g['span'],g['gram_identity'])<=1e-10 for g in gates[:2]),
                resources=time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3)
            records.append(dict(kappa=item['kappa'],n=item['n'],amplitude=item['amplitude'],cell=item['cell'],
                preparations=preparations,precision_agreement=agreement,gates=gates,checks={k:bool(v) for k,v in checks.items()},
                source_artifact=item['artifact'],source_sha256=item['artifact_sha256'],
                artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path)))
            print('gate',index,gates,checks,flush=True)
            if not all(checks.values()):raise RuntimeError('fixed_GSP_basis_v2_gate_failed')
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(v2.__file__),Path(v2.old.__file__),parentpath,
        ROOT/'rd/structure/main-shell/gsp-q3/gsp-q3-basis-v2-contract.md',ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-basis-v2-protocol.md',
        ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-basis-v2-preparation.json']
    dest.write_text(json.dumps(dict(records=records,all_checks_pass=len(records)==2 and stop is None,
        stop_reason=stop,local_30_executed=False,physical_adoption=False,seconds=time.monotonic()-start,
        peak_bytes=memory.memory_bytes(),sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources}),indent=2),encoding='utf-8')

if __name__=='__main__':main()
