"""Independent quadrature/derivatives with change of geometric face coordinates."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,importlib,json,time
import numpy as np
import check_reference_work_local as checker
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a-b)/max(1e-30,np.linalg.norm(b)))
def main():
    for module,prefix in [('reference_work_ancf','reference-work'),('reference_work_ancf_v2','reference-work-v2')]:
        f=importlib.import_module(module);checker.f=f;source=ROOT/rd_paths.record(f'{prefix}-local-result.json');parent=json.loads(source.read_text());dest=ROOT/rd_paths.record(f'{prefix}-local-audit.json');assert not dest.exists()
        records=[];stops=[];start=time.monotonic()
        for r in parent['records']:
            try:
                path=ROOT/r['artifact'];assert digest(path)==r['artifact_sha256'];data=np.load(path);s=f.prepare(data['q'],r['h']);si=f.prepare(data['q'],r['h'],order=15,zorder=9,independent=True)
                assert s['nscalar']==si['nscalar'];T=np.eye(s['nscalar'])
                for fm,fi in zip(s['faces'],si['faces']):
                    sl=slice(fi['offset'],fi['offset']+fi['n']);M=f.face_values(s,fm,fi['points'],True);T[sl,sl]=fi['hat'].T@((fi['w']/4)[:,None]*M)
                Xi=T@data['current'];a=(T@data['v1'].reshape(-1,3)).reshape(-1);b=(T@data['v2'].reshape(-1,3)).reshape(-1)
                U,force,K,info=f.nonlinear(si,Xi,a,b);back=(T.T@force.reshape(-1,3)).reshape(-1)
                metrics=dict(energy=rel(U,data['energy']),force=rel(back,data['force']),tangent=rel(K,data['tangent']),face_transform=float(np.linalg.norm(T.T@T-np.eye(len(T)))),
                    gradient=rel(si['Gw']@T,f.evaluate_gradient(s,si['p'],True)))
                wk=checker.work(si);poly=checker.patch(si) if r['geometry']=='flat' else None
                records.append(dict(geometry=r['geometry'],h=r['h'],metrics=metrics,work=wk,polynomial=poly,reference_gradient=si['reference_gradient_error'],gram=si['gram_error'],
                    passed=bool(max(metrics.values())<=1e-8 and max(v for k,v in wk.items() if k!='fields')<=1e-9 and (poly is None or poly['passed']))))
            except Exception as exc:stops.append(dict(geometry=r['geometry'],h=r['h'],reason=repr(exc)))
        reproduced_stops=[]
        for r in parent['stops']:
            q=checker.prior.generic()
            try:f.prepare(q,r['h'],order=15,zorder=9,independent=True);reproduced_stops.append(dict(geometry=r['geometry'],h=r['h'],failed_again=False))
            except Exception as exc:reproduced_stops.append(dict(geometry=r['geometry'],h=r['h'],failed_again=True,reason=repr(exc)))
        paths=[Path(__file__),Path(f.__file__),source,ROOT/'rd/structure/main-shell/reference-work/records/reference-work-protocol.md']
        dest.write_text(json.dumps(dict(records=records,stops=stops,parent_stops_rechecked=reproduced_stops,all_reproduction_checks=len(records)==len(parent['records']) and not stops and all(x['passed'] for x in records),
            independent_global_resolve=False,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
        print(prefix,'records',len(records),'failed',[(r['geometry'],r['h'],r['metrics']) for r in records if not r['passed']],'stops',stops,flush=True)
if __name__=='__main__':main()
