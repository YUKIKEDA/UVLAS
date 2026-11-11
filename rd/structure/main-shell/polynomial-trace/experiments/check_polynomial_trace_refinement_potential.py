"""Original four finite-h quality cases after the global comparison gate."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
import polynomial_trace_fit_v2 as fit
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    dest=ROOT/'rd/structure/main-shell/polynomial-trace/records/polynomial-trace-potential-refinement-result.json';assert not dest.exists()
    parent=json.loads((ROOT/'rd/structure/main-shell/polynomial-trace/records/polynomial-trace-shared-result.json').read_text());decision=next(r for r in parent['dispositions'] if r['kind']=='potential')
    assert decision['all_global_checks'] and decision['all_curved_improved']
    out=ROOT/'build/analysis'/('polynomial-trace-refinement-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();start=time.monotonic();records=[];stop=None
    try:
        for n,h in [(4,.01),(4,.0001),(8,.001),(8,.0001)]:
            tick=time.monotonic();a=fit.assemble(.4,n,.2,h,'potential');solution=fit.solve(a);observed=fit.observe(a,solution)
            reference=[fit.m.reference_h1(kind,.4,40) for kind in fit.KIND_TARGETS]
            checks=dict(local_condition=solution['fit']['local_rcond']>=1e-10,global_condition=solution['fit']['global_rcond']>=1e-10,
                local_ranks=all(c['local_ranks']==[399]*3 for c in a['local']),linear=max(solution['fit']['linear_residual'])<=1e-8,kkt=max(solution['fit']['kkt'])<=1e-8,
                quadrature=max(abs(x-y) for x,y in zip(observed['error'],solution['fit']['error']))<=1e-8,H1=all(x<=2*y*(1+1e-8) for x,y in zip(observed['h1'],reference)),
                resources=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3)
            path=out/f'potential-k0.4-a0.2-n{n}-h{h:g}.npz';np.savez_compressed(path,local=solution['local'],trace=solution['trace'],ref=a['original']['ref'],cells=a['original']['cells'])
            records.append(dict(kind='potential',kappa=.4,amplitude=.2,n=n,h=h,fit=solution['fit'],observed=observed,checks={k:bool(v) for k,v in checks.items()},
                seconds=time.monotonic()-tick,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path)))
            (out/'checkpoint.json').write_text(json.dumps(dict(records=records),indent=2),encoding='utf-8')
            print(n,h,'error',observed['error'],'H1',observed['h1'],'fit',solution['fit'],'checks',checks,'seconds',time.monotonic()-tick,flush=True)
            del a,solution
            if not all(checks.values()):raise RuntimeError('refined_numerical_or_H1_gate_failed')
    except Exception as exc:stop=repr(exc)
    complete=len(records)==4 and stop is None;ratios=[records[1]['observed']['error'][j]/records[0]['observed']['error'][j] for j in range(2)] if len(records)>=2 else None
    accuracy=complete and all(r['observed']['error'][0]<=.01 for r in records if r['n']==8) and ratios[0]<=2
    sources=[Path(__file__),Path(fit.__file__),ROOT/'rd/structure/main-shell/hybrid-gradient/experiments/hybrid_gradient_fit.py',ROOT/'rd/structure/main-shell/hybrid-gradient/experiments/hybrid_gradient_ancf.py',ROOT/'rd/structure/main-shell/polynomial-trace/records/polynomial-trace-protocol.md',ROOT/'rd/structure/main-shell/polynomial-trace/records/polynomial-trace-potential-refinement-preparation.json',ROOT/'rd/structure/main-shell/polynomial-trace/records/polynomial-trace-shared-result.json']
    dest.write_text(json.dumps(dict(records=records,stop_reason=stop,all_numerical_checks=complete,quartic_accuracy_gate_pass=bool(accuracy),thin_amplification=ratios,
        exponential_additional_observation_only=True,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
