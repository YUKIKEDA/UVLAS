"""Membrane/bending/fixed energy budgets on two meshes and three thicknesses."""
from pathlib import Path
import itertools,json,time
import numpy as np
from scipy.linalg import solve_triangular
import q9_curved_energy_scaling as method
import check_q9_curved_best_approximation_bounds as bounds
core=method.core;ROOT=core.ROOT;F=core.F
def position_error(a,r):return float(np.sqrt(max(0.,r@a['G']@r-2*r@a['b'][:,0]+a['norm'][0])/a['norm'][0]))
def main():
    target=F/'q9-curved-energy-scaling-result.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-curved-best-approximation-decimal-v2-audit.json';assert json.loads(parent.read_text())['all_criteria_pass'];start=time.monotonic();references={};cases=[];stop=None;original=core.method.prepare;out=ROOT/'build/analysis'/('q9-curved-energy-scaling-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    try:
        for name in ['RECT','TRAP','SADDLE','NORMAL-C']:
            r12=method.reference(name,12,[.01,.001,.0001]);r16=method.reference(name,16,[.01,.001,.0001]);error=abs(r12['bending']/r16['bending']-1)
            if name in ['SADDLE','NORMAL-C']:error=max(error,abs(r12['membrane']/r16['membrane']-1))
            else:assert r16['membrane']<=1e-24*r16['bending']
            if name=='RECT':assert abs(r16['bending']/((core.shapes.lam+2*method.MU)*np.pi**4)-1)<=1e-11
            assert error<=1e-10 and r16['derivative_error']<=1e-11;references[name]=dict(order12=r12,order16=r16,relative_difference=error)
        tests=bounds.synthetic();assert all(t['contains'] and t['feasible'] and t['objective_order'] for t in tests)
        configs=[(name,n,h) for name,n,h in itertools.product(['SADDLE','NORMAL-C'],[2,4],[.01,.001,.0001])]+[(name,n,.001) for name,n in itertools.product(['RECT','TRAP'],[2,4])]
        for name,n,h in configs:
            cache={}
            def prepare(q,width,offset,order):
                key=(q.tobytes(),width,offset,order)
                if key not in cache:cache[key]=original(q,width,offset,order)
                return cache[key]
            core.method.prepare=prepare;a=method.assemble(name,n,h);problem=method.ReducedProblem(a);reference=references[name]['order16'];Im=reference['membrane'];Ib=reference['bending'];Um=h*Im/2;Ub=h**3*Ib/24;ref_h=.01;fixed=np.sqrt((1e-2*ref_h**3*Ib/12)/(method.MU*reference['analytic_volume'][str(ref_h)]));bound=2*core.pnorm(a['reference'][:,:,0],a['M']);r0=solve_triangular(problem.Lr.T,problem.z0[:,0],lower=False);base_error=position_error(a,r0);forcing=float(np.linalg.norm(problem.C@problem.z0[:,0]));settings=[]
            if name in ['SADDLE','NORMAL-C']:
                for eta in [1e-2,1e-4]:settings.append(('membrane',eta,np.sqrt(2*eta*Um/(method.MU*a['volume']))))
            for eta in [1e-2,1e-4]:settings.append(('bending',eta,np.sqrt(2*eta*Ub/(method.MU*a['volume']))))
            settings.append(('fixed',1e-2,fixed));records=[];saved_r=[];saved_p=[]
            for family,eta,epsilon in settings:
                fit=problem.solve(0,bound,float(np.nextafter(epsilon/forcing,0.)*(1-8*np.finfo(float).eps)));r=fit.pop('r');p=fit.pop('p');saved_r.append(r);saved_p.append(p);g=a['A']@p+a['W']@r;Hp=(a['M']@p.reshape(len(a['M']),6)).reshape(-1);grad_r=a['G']@r-a['b'][:,0]+fit['lambda_shear']*a['W'].T@g;grad_p=fit['lambda_shear']*a['A'].T@g+fit['mu_h1']*Hp;pr=float(np.linalg.norm(grad_r)/max(np.linalg.norm(a['b'][:,0])+fit['lambda_shear']*np.linalg.norm(a['W'].T@g),1e-30));pp=float(np.linalg.norm(grad_p)/max(fit['lambda_shear']*np.linalg.norm(a['A'])*np.linalg.norm(g)+fit['mu_h1']*np.linalg.norm(Hp),1e-30)) if fit['lambda_shear'] else 0.
                fit.update(family=family,eta=eta,requested_epsilon=float(epsilon),relative_epsilon=float(epsilon/forcing),position_error=position_error(a,r),baseline_error=base_error,physical_kkt_r=pr,physical_kkt_p=pp,component_energy9=float(.5*method.MU*a['volume']*fit['residual']**2));fit['kkt_pass']=bool(max(fit['kkt_r'],fit['kkt_p'],pr,pp,fit['complementarity'])<=1e-7);records.append(fit)
            del problem
            observation=method.assemble(name,n,h,observation=11,position_order=11);combined=dict(a,A11=observation['A'],W11=observation['W'],G11=observation['G'],b11=observation['b'],norm11=observation['norm']);residuals11=[]
            for record,r,p in zip(records,saved_r,saved_p):
                v=observation['A']@p+observation['W']@r;residuals11.append(v);norm=float(np.linalg.norm(v));record['residual11']=norm;record['observation_difference']=abs(norm-record['residual'])/max(norm,record['residual'],1e-8*forcing);record['position_error11']=position_error(observation,r);record['position_observation_difference']=abs(record['position_error11']-record['position_error']);record['bounds']=bounds.certify(combined,r,p,record['lambda_shear'],record['epsilon'],record['bound'],0);record['certified_component_energy_upper9']=float(.5*method.MU*a['volume']*record['bounds']['shear_upper9']**2);record['certified_component_energy_upper11']=float(.5*method.MU*observation['volume']*record['bounds']['shear_upper11']**2)
                print(json.dumps(dict(geometry=name,n=n,h=h,family=record['family'],eta=record['eta'],epsilon=record['requested_epsilon'],relative_epsilon=record['relative_epsilon'],error_interval=[record['bounds']['lower_relative_error'],record['bounds']['upper_relative_error']],kkt=record['kkt_pass'])),flush=True)
            checks=dict(reference=references[name]['relative_difference']<=1e-10,energy_identity=max(a['energy_identity_error'],observation['energy_identity_error'])<=1e-11,edge=a['edge_valid'] and observation['edge_valid'],metric=min(a['metric_min'],observation['metric_min'])>0,primal=all(r['primal']<=1e-7 for r in records),observation=all(r['observation_difference']<=1e-4 for r in records),position_observation=all(r['position_observation_difference']<=1e-7 for r in records),bounds=all(r['bounds']['feasible'] and r['bounds']['objective_order'] for r in records))
            path=out/f'{name}-n{n}-h{h:g}.npz';np.savez_compressed(path,A=a['A'],W=a['W'],M=a['M'],G=a['G'],b=a['b'],norm=a['norm'],reference=a['reference'],xy=a['xy'],interior=a['interior'],r=np.array(saved_r),p=np.array(saved_p),residuals11=np.array(residuals11),G11=observation['G'],b11=observation['b'],norm11=observation['norm']);cases.append(dict(geometry=name,n=n,h=h,volume9=a['volume'],volume11=observation['volume'],membrane_scale=Um,bending_scale=Ub,fixed_epsilon=float(fixed),baseline_error=base_error,energy_identity_error=max(a['energy_identity_error'],observation['energy_identity_error']),metric_eigenvalues=[min(a['metric_min'],observation['metric_min']),max(a['metric_max'],observation['metric_max'])],edge_condition=max(a['edge_condition'],observation['edge_condition']),records=records,checks={k:bool(v) for k,v in checks.items()},all_kkt_pass=bool(all(r['kkt_pass'] for r in records)),artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=core.digest(path)));core.method.prepare=original
            print(json.dumps(dict(case=[name,n,h],checks=cases[-1]['checks'],seconds=time.monotonic()-start)),flush=True)
            del a,observation,combined,cache
            if time.monotonic()-start>900 or core.memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
        assert sum(p.stat().st_size for p in out.iterdir())<=256*1024**2
    except Exception as exc:stop=repr(exc)
    finally:core.method.prepare=original
    complete=bool(stop is None and len(cases)==16 and sum(len(c['records']) for c in cases)==72);result=dict(source_sha256=core.digest(__file__),parent_sha256=core.digest(parent),contract_sha256=core.digest(ROOT/'docs/q9-curved-energy-scaling-contract.md'),protocol_sha256=core.digest(F/'q9-curved-energy-scaling-experiment.md'),dependencies={Path(m.__file__).name:core.digest(m.__file__) for m in [method,core,core.method,core.method.base,core.geometry,bounds]},material=dict(C=method.CM.tolist(),mu_ref=method.MU),references=references,cases=cases,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=core.memory.memory_bytes(),bounds_checks_pass=bool(complete and all(all(c['checks'].values()) for c in cases)),all_kkt_pass=bool(complete and all(c['all_kkt_pass'] for c in cases)),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['cases','references','dependencies','material']}))
    if not result['bounds_checks_pass']:raise SystemExit(1)
if __name__=='__main__':main()
