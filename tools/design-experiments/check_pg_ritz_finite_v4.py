"""Complete general-surface cases with coefficient-exact, conditioned full traces."""
from pathlib import Path
import gc,hashlib,json,time
import numpy as np
import check_pg_ritz_finite_v3 as engine
import pg_ritz_shell_v3 as shell
import pg_ritz_trace_v2 as trace
engine.shell=shell;engine.trace=trace
ROOT=Path(__file__).resolve().parents[2]

def main():
    dest=ROOT/'docs/verification/pg-ritz-finite-v4-result.json';assert not dest.exists();parent=ROOT/'docs/verification/pg-ritz-finite-v3-result.json';old=json.loads(parent.read_text());records=[dict(r,execution_route='cylinder_reference_branch_v2') for r in old['records']];assert len(records)==4;stops=[];out=ROOT/'build/analysis'/('pg-ritz-general-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();start=time.monotonic()
    for parameters,n,variable in [((.4,.2,0.),2,False),((.4,.2,0.),4,False),((.4,0.,.07),2,False),((.4,0.,.07),4,False),((.4,-.2,.07),2,True),((.4,-.2,.07),4,True)]:
        try:
            model=shell.prepare(n,parameters);print('PREPARED',parameters,n,model['ns'],flush=True)
            for h in [.01,.0001]:
                tick=time.monotonic();sol=engine.solve(model,h,variable);o=engine.observe(model,h,sol,variable);p=np.array([[0.,y,z] for y in np.linspace(0,1,113) for z in [-1.,-.31,0.,.47,1.]]);N,G,B,det,X=shell.kinematics(model,p,h,variable);_,Ug=engine.coarse.fields(model,shell,p,h,variable);u=(N@sol['coef'].reshape(model['ns'],9)).reshape(len(p),3,3)+Ug@sol['beta'];fixed=float(np.max(abs(u)));trace_defect=float(np.linalg.norm(N@sol['Zs'])/max(1.,np.linalg.norm(N)));checks=dict(stationarity=max(o['stationarity'])<=1e-8,condition=sol['rcond']>=1e-10,mass=sol['mass_min']>1e-10,fixed=fixed<=1e-9,trace=trace_defect<=1e-9)
                artifact=out/f'p{parameters}-n{n}-h{h:g}.npz';np.savez_compressed(artifact,coef=sol['coef'],beta=sol['beta'],Zs=sol['Zs'],ritz_T=model['lift']['T'],ritz_coeff=model['lift']['coeff'],normal_coeff=model['normal_completion']['coeff'],normal_V=model['normal_completion']['V'],normal_s=model['normal_completion']['s'],extra_coeff=model['completion']['coeff'],extra_V=model['completion']['V'],extra_s=model['completion']['s']);record=dict(execution_route='general_reference_clamped_v3',parameters=parameters,n=n,h=h,variable=variable,observation=o,rcond=sol['rcond'],mass_min=sol['mass_min'],fixed=fixed,trace_defect=trace_defect,checks=checks,seconds=time.monotonic()-tick,artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(artifact.read_bytes()).hexdigest());records.append(record);(out/'checkpoint.json').write_text(json.dumps(dict(records=records,stops=stops),indent=2),encoding='utf-8');print(record,flush=True);del sol;gc.collect()
                if not all(checks.values()):raise RuntimeError('general_finite_numerical_gate')
        except Exception as exc:stops.append(dict(parameters=parameters,n=n,reason=repr(exc)));print('STOP',stops[-1],flush=True);break
    paths=[parent,Path(__file__),Path(shell.__file__),Path(shell.lift.__file__),Path(shell.lift.clamp.__file__),Path(trace.__file__),ROOT/'tools/design-experiments/pg_trace_nullspace.py',Path(engine.__file__)];result=dict(all_cases=len(records)==16 and not stops,records=records,stops=stops,R4_complete=False,main_shell_adopted=False,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest.write_text(json.dumps(result,indent=2),encoding='utf-8')

if __name__=='__main__':main()
