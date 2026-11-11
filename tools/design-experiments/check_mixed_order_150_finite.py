from pathlib import Path
import hashlib,json,time
import numpy as np
import mixed_order_150_finite as f
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    dest=ROOT/'docs/verification/mixed-order-150-finite-result.json';assert not dest.exists()
    parent=json.loads((ROOT/'docs/verification/mixed-order-150-result.json').read_text());audit=json.loads((ROOT/'docs/verification/mixed-order-150-audit.json').read_text())
    assert parent['all_numerical_checks'] and parent['both_curved_improved'] and audit['all_reproduction_checks']
    start=time.monotonic();out=ROOT/'build/analysis'/('mixed-order-150-finite-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();records=[];stop=None
    try:
        for n,h in [(4,.01),(4,.0001),(8,.001),(8,.0001)]:
            tick=time.monotonic();a=f.assemble(n,h);u,fit=f.solve(a);observed=f.observe(a['mesh'],h,u,True)
            rh=f.m.reference_h1('quartic',.4,40);bound=2*rh
            checks=dict(local_condition=fit['local_rcond']>=1e-10,global_condition=fit['global_rcond']>=1e-10,
                linear=fit['linear_residual']<=1e-8,kkt=fit['kkt']<=1e-8,quadrature=abs(observed['error']-fit['error'])<=1e-8,
                H1=abs(observed['h1']/fit['h1']-1)<=1e-8,sphere=observed['h1']<=bound*(1+1e-8),
                resources=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3)
            path=out/f'n{n}-h{h:g}.npz';mesh=a['mesh'];np.savez_compressed(path,ref=mesh['ref'],cells=mesh['cells'],u=u,
                xy0=mesh['xy'][0],xy1=mesh['xy'][1],xy2=mesh['xy'][2],free0=mesh['free'][0],free1=mesh['free'][1],free2=mesh['free'][2])
            records.append(dict(n=n,h=h,fit=fit,observed=observed,h1_bound=bound,full_free_dofs=3*len(a['free']),boundary_dofs=3*len(a['boundary']),
                checks={k:bool(v) for k,v in checks.items()},seconds=time.monotonic()-tick,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path)))
            print(n,h,fit,observed,checks,flush=True)
            del a,u
            if not all(checks.values()):raise RuntimeError('finite_thickness_numerical_or_H1_gate_failed')
    except Exception as exc:stop=repr(exc)
    complete=len(records)==4 and stop is None;thin_ratio=records[1]['observed']['error']/records[0]['observed']['error'] if len(records)>1 else None
    accuracy=complete and all(r['observed']['error']<=.01 for r in records if r['n']==8) and thin_ratio<=2
    sources=[Path(__file__),Path(f.__file__),Path(f.m.__file__),Path(f.reference.__file__),Path(f.qrtools.__file__),ROOT/'docs/verification/mixed-order-150-finite-protocol.md',ROOT/'docs/verification/mixed-order-150-finite-preparation.json',ROOT/'docs/verification/mixed-order-150-result.json',ROOT/'docs/verification/mixed-order-150-audit.json']
    dest.write_text(json.dumps(dict(records=records,all_numerical_checks=complete,accuracy_gate_pass=bool(accuracy),thin_amplification=thin_ratio,
        stop_reason=stop,physical_adoption=False,finite_thickness_lower_bound=False,unconstrained_solution_admissible=complete,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
