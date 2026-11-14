"""Re-solve selected H1 balls at volume order 9 and observe unchanged at 11."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,itertools,time
import numpy as np
from scipy.linalg import svd
import check_q9_curved_bounded_lift as core
ROOT=core.ROOT;F=core.F
def main():
    target=F/'q9-curved-bounded-lift-fine-result.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-curved-bounded-lift-result.json';r=json.loads(parent.read_text());assert len(r['cases'])==22 and r['stop_reason'] is None and not r['holdout']['checks_pass'];start=time.monotonic();cases=[];stop=None;original=core.method.prepare;out=ROOT/'build/analysis'/('q9-curved-bounded-lift-fine-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    try:
        for name,h in itertools.product(['SADDLE','NORMAL-C'],[.1,.01,.001]):
            cache={}
            def prepare(q,width,offset,order):
                key=(q.tobytes(),width,offset,order)
                if key not in cache:cache[key]=original(q,width,offset,order)
                return cache[key]
            core.method.prepare=prepare
            a=core.assemble(name,2,h,observation=9);Aw,chol=core.whiten(a['A'],a['M']);U,s,Vh=svd(Aw,full_matrices=False,check_finite=False);rows=[];solutions=[];coarse=next(c for c in r['cases'] if c['geometry']==name and c['n']==2 and c['h']==h);old=np.load(ROOT/coarse['artifact'])['bounded']
            for k,wave in enumerate(core.WAVES):
                f=a['F'][:,k];forcing=float(np.linalg.norm(f));reference=core.pnorm(a['reference'][:,:,k],a['M']);fits=[];free=[];rhs=U.T@f
                for threshold in [1e-9,1e-10,1e-11]:
                    keep=s>threshold*s[0];coef=np.zeros_like(s);coef[keep]=-rhs[keep]/s[keep];v=Vh.T@coef;free.append(dict(threshold=threshold,rank=int(keep.sum()),relative=float(np.linalg.norm(Aw@v+f)/forcing),h1=float(np.linalg.norm(v)),direction_h1=float(np.linalg.norm(v.reshape(len(chol),6)[:,:3])),curvature_h1=float(np.linalg.norm(v.reshape(len(chol),6)[:,3:]))))
                for factor in [2.,4.]:
                    fit=core.bounded(U,s,Vh,f,factor*reference,Aw,chol);solutions.append(fit.pop('p'));fit.update(bound_factor=factor,bound=factor*reference,relative=fit['residual']/forcing,thin_scaled=fit['residual']/forcing/h);fits.append(fit)
                rows.append(dict(wave=list(wave),forcing=forcing,reference_h1=reference,unbounded=free,bounded=fits))
            observation=core.assemble(name,2,h,observation=11);differences=[]
            for index,p in enumerate(solutions):
                k=index//2;fit=rows[k]['bounded'][index%2];norm11=float(np.linalg.norm(observation['A']@p+observation['F'][:,k]));coarse11=float(np.linalg.norm(observation['A']@old[index]+observation['F'][:,k]));forcing11=float(np.linalg.norm(observation['F'][:,k]));error=abs(norm11-fit['residual'])/max(norm11,fit['residual'],1e-30);differences.append(error);fit.update(residual11=norm11,relative11=norm11/forcing11,coarse_solution_residual11=coarse11,coarse_solution_relative11=coarse11/forcing11,observation_difference=error)
            checks=dict(analytic=a['analytic']<=1e-11,constraint=all(fit['constraint_error']<=1e-10 for row in rows for fit in row['bounded']),kkt=all(fit['kkt']<=1e-9 for row in rows for fit in row['bounded']),observation=max(differences)<=1e-4)
            path=out/f'{name}-h{h:g}.npz';np.savez_compressed(path,A=a['A'],F=a['F'],M=a['M'],reference=a['reference'],bounded=np.array(solutions),A11=observation['A'],F11=observation['F'],singular_values=s)
            cases.append(dict(geometry=name,n=2,h=h,rows=rows,maximum_observation_difference=max(differences),checks={k:bool(v) for k,v in checks.items()},artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=core.digest(path)));print(json.dumps(dict(geometry=name,h=h,checks=cases[-1]['checks'],relative11=[row['bounded'][0]['relative11'] for row in rows],observation=max(differences))),flush=True)
            core.method.prepare=original
            if time.monotonic()-start>300 or core.memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
    except Exception as exc:stop=repr(exc)
    finally:core.method.prepare=original
    result=dict(source_sha256=core.digest(__file__),parent_sha256=core.digest(parent),protocol_sha256=core.digest(F/'q9-curved-bounded-lift-fine-experiment.md'),dependencies={Path(m.__file__).name:core.digest(m.__file__) for m in [core,core.method,core.method.base,core.geometry]},cases=cases,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=core.memory.memory_bytes(),diagnostic_checks_pass=bool(stop is None and len(cases)==6 and all(all(c['checks'].values()) for c in cases)),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['cases','dependencies']}))
    if not result['diagnostic_checks_pass']:raise SystemExit(1)
if __name__=='__main__':main()
