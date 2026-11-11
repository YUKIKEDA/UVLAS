"""Shared all-component interior position approximation: predeclared 152 QCQPs."""
from pathlib import Path
import itertools,json,time
import numpy as np
from scipy.linalg import solve_triangular
import q9_curved_best_approximation as method
core=method.core;ROOT=core.ROOT;F=core.F;EPS=[1e-3,1e-4,1e-5,1e-6]
def error(a,r,k):return float(np.sqrt(max(0.,r@a['G']@r-2*r@a['b'][:,k]+a['norm'][k])/a['norm'][k]))
def main():
    target=F/'q9-curved-best-approximation-result.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-curved-bounded-lift-fine-audit.json';assert json.loads(parent.read_text())['all_criteria_pass'];start=time.monotonic();out=ROOT/'build/analysis'/('q9-curved-best-approximation-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();cases=[];stop=None;original=core.method.prepare
    configs=[(name,2,h) for name,h in itertools.product(['SADDLE','NORMAL-C'],[.01,.001])]+[(name,2,.001) for name in ['RECT','TRAP']]+[(name,4,.001) for name in ['SADDLE','NORMAL-C']]
    try:
        for name,n,h in configs:
            cache={}
            def prepare(q,width,offset,order):
                key=(q.tobytes(),width,offset,order)
                if key not in cache:cache[key]=original(q,width,offset,order)
                return cache[key]
            core.method.prepare=prepare;a=method.assemble(name,n,h);problem=method.ReducedProblem(a);records=[];saved_r=[];saved_p=[]
            for k in ([0] if n==4 else range(3)):
                baseline_r=solve_triangular(problem.Lr.T,problem.z0[:,k],lower=False);baseline=error(a,baseline_r,k);ref=core.pnorm(a['reference'][:,:,k],a['M'])
                for factor in ([2.] if n==4 else [2.,4.]):
                    for eps in EPS:
                        fit=problem.solve(k,factor*ref,eps);r=fit.pop('r');p=fit.pop('p');saved_r.append(r);saved_p.append(p);gamma=a['A']@p+a['W']@r;Hp=(a['M']@p.reshape(len(a['M']),6)).reshape(-1);grad_r=a['G']@r-a['b'][:,k]+fit['lambda_shear']*(a['W'].T@gamma);grad_p=fit['lambda_shear']*(a['A'].T@gamma)+fit['mu_h1']*Hp
                        physical_r=float(np.linalg.norm(grad_r)/max(np.linalg.norm(a['b'][:,k])+fit['lambda_shear']*np.linalg.norm(a['W'].T@gamma),1e-30));physical_p=float(np.linalg.norm(grad_p)/max(fit['lambda_shear']*np.linalg.norm(a['A'].T@gamma)+fit['mu_h1']*np.linalg.norm(Hp),1e-30)) if fit['lambda_shear'] else 0.
                        fit.update(wave=list(core.WAVES[k]),wave_index=k,bound_factor=factor,relative_tolerance=eps,position_error=error(a,r,k),baseline_error=baseline,physical_kkt_r=physical_r,physical_kkt_p=physical_p);records.append(fit)
                        print(json.dumps(dict(geometry=name,n=n,h=h,wave=fit['wave'],bound=factor,epsilon=eps,error=fit['position_error'],primal=fit['primal'],kkt=max(fit['kkt_r'],fit['kkt_p'],physical_r,physical_p))),flush=True)
                        if time.monotonic()-start>600:raise RuntimeError('resource_time_limit')
            del problem
            independent=method.assemble(name,n,h,observation=11,position_order=11)
            for record,r,p in zip(records,saved_r,saved_p):
                gamma=float(np.linalg.norm(independent['A']@p+independent['W']@r));record['observation_difference']=abs(gamma-record['residual'])/max(gamma,record['residual'],1e-8*record['forcing']);record['residual11']=gamma;record['relative_shear11']=gamma/record['forcing'];record['independent_constraint_error']=max(0.,gamma/record['epsilon']-1);record['position_error11']=error(independent,r,record['wave_index']);record['position_observation_difference']=abs(record['position_error11']-record['position_error'])
            monotone=True;bounds=True
            for k in ([0] if n==4 else range(3)):
                for factor in ([2.] if n==4 else [2.,4.]):
                    group=[r for r in records if r['wave_index']==k and r['bound_factor']==factor]
                    monotone=monotone and all(b['position_error']>=a0['position_error']*(1-1e-6) for a0,b in zip(group,group[1:]))
                if n!=4:
                    for eps in EPS:
                        group=[r for r in records if r['wave_index']==k and r['relative_tolerance']==eps];bounds=bounds and group[1]['position_error']<=group[0]['position_error']*(1+1e-6)
            checks=dict(analytic=a['analytic']<=1e-11,primal=all(r['primal']<=1e-7 for r in records),kkt=all(max(r['kkt_r'],r['kkt_p'],r['physical_kkt_r'],r['physical_kkt_p'],r['complementarity'])<=1e-7 for r in records),observation=all(r['observation_difference']<=1e-4 for r in records),position_observation=all(r['position_observation_difference']<=1e-7 for r in records),independent_feasible=all(r['independent_constraint_error']<=1e-4 for r in records),monotone=bool(monotone and bounds))
            path=out/f'{name}-n{n}-h{h:g}.npz';np.savez_compressed(path,A=a['A'],W=a['W'],M=a['M'],G=a['G'],b=a['b'],norm=a['norm'],interior=a['interior'],xy=a['xy'],r=np.array(saved_r),p=np.array(saved_p),A11=independent['A'],W11=independent['W'],G11=independent['G'],b11=independent['b'],norm11=independent['norm'])
            cases.append(dict(geometry=name,n=n,h=h,position_unknowns=a['W'].shape[1],direction_unknowns=a['A'].shape[1],records=records,checks={k:bool(v) for k,v in checks.items()},artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=core.digest(path)));print(json.dumps(dict(case=[name,n,h],checks=cases[-1]['checks'])),flush=True);core.method.prepare=original
            if time.monotonic()-start>600 or core.memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
            del a,independent,cache
        assert sum(p.stat().st_size for p in out.iterdir())<=128*1024**2
    except Exception as exc:stop=repr(exc)
    finally:core.method.prepare=original
    result=dict(source_sha256=core.digest(__file__),parent_sha256=core.digest(parent),contract_sha256=core.digest(ROOT/'docs/q9-curved-best-approximation-contract.md'),protocol_sha256=core.digest(F/'q9-curved-best-approximation-experiment.md'),dependencies={Path(m.__file__).name:core.digest(m.__file__) for m in [method,core,core.method,core.method.base,core.geometry]},cases=cases,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=core.memory.memory_bytes(),diagnostic_checks_pass=bool(stop is None and len(cases)==8 and sum(len(c['records']) for c in cases)==152 and all(all(c['checks'].values()) for c in cases)),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['cases','dependencies']}))
    if not result['diagnostic_checks_pass']:raise SystemExit(1)
if __name__=='__main__':main()
