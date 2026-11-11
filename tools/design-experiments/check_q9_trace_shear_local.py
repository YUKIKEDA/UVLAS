"""Frozen local gates for T-S36; no global solve until all gates pass."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import q9_trace_shear as candidate
import check_q9_shear_projection_gate as baseline
import check_q9_curved_load as geometry
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=F/'q9-trace-shear-local-result.json'
    if output.exists():raise FileExistsError(output)
    out=ROOT/'build/analysis'/('q9-trace-shear-local-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();start=time.monotonic();records=[];stop=None
    try:
        for name in ['RECT','SKEW','TRAP']:
            q=geometry.geometry(name);data={o:candidate.matrices(q,.1,o) for o in [7,9]};d=data[9];scale=1/np.sqrt(d['Kfull'].diagonal());scaled=lambda k:scale[:,None]*k*scale[None,:];K=scaled(d['Ktie']);norm=np.linalg.norm(K,2);singular=np.linalg.svd(K,compute_uv=False);ranks={str(t):int(np.count_nonzero(singular>t*singular[0])) for t in [1e-11,1e-10,1e-9]}
            fields=baseline.patch_fields(q);patches=[float(np.linalg.norm(scale*((d['Ktie']-d['Kfull'])@u))/np.linalg.norm(scale*(d['Kfull']@u))) for u in fields]
            coords=q.reshape(9,3,3);rigid=[]
            for axis in np.eye(3):
                u=np.zeros_like(coords);u[:,0,:]=axis;rigid.append(u.reshape(81));rigid.append(np.cross(np.broadcast_to(axis,coords.shape),coords).reshape(81))
            re=max(float(np.linalg.norm(K@(u/scale))/(norm*np.linalg.norm(u/scale))) for u in rigid);integration=float(np.linalg.norm(scaled(data[7]['Ktie']-d['Ktie']))/np.linalg.norm(scaled(d['Kfull'])))
            symmetry=float(np.linalg.norm(K-K.T)/np.linalg.norm(K));minimum=float(np.linalg.eigvalsh((K+K.T)/2)[0]);condition=float(np.linalg.cond(d['factor'][1]));u=scale*np.random.default_rng(831).normal(size=81);beta=d['T']@u
            energy=abs(float((u@d['Kn']@u+beta@d['M']@beta-u@d['Ktie']@u)/2))/max(1.,abs(float(u@d['Ktie']@u/2)));force=float(np.linalg.norm(scale*(d['Kn']@u+d['T'].T@(d['M']@beta)-d['Ktie']@u))/np.linalg.norm(scale*(d['Ktie']@u)))
            duplicate=d['M'].copy();duplicate[-1,:]=duplicate[0,:];duplicate[:,-1]=duplicate[:,0];rejected=False
            try:baseline.factor(duplicate)
            except ValueError:rejected=True
            # Independent constant physical shear values, after reconstruction at volume points.
            p=np.array([[-.73,-.41,-.2],[.22,.37,.8],[.71,-.63,.3]]);J,_=candidate.kinematics(q,.1,p);V=np.linalg.solve(J[:,:2,:2].transpose(0,2,1),candidate.interpolation(p));constant=0.
            for axis in range(2):
                b=d['T']@fields[4 if axis==0 else 3];value=np.einsum('pij,j->pi',V,b);expected=np.eye(2)[axis]
                constant=max(constant,float(np.max(abs(value-expected))))
            checks=dict(gram=condition<=1e8,symmetry=symmetry<=1e-12,integration=integration<=1e-8,rank=all(v==75 for v in ranks.values()),rigid=re<=1e-11,
                positive=minimum>=-1e-10*max(1.,norm),patch=max(patches)<=1e-10,constant_shear=constant<=1e-12,energy=energy<=1e-10,force=force<=1e-10,duplicate_rejected=rejected,leading_metric=d['leading_metric_difference']<=1e-12)
            checks={k:bool(v) for k,v in checks.items()}
            file=out/(name+'.npz');np.savez_compressed(file,**{k:v for k,v in d.items() if k!='factor'},q=q,rigid=np.array(rigid),patch_fields=np.array(fields))
            records.append(dict(case=name,ranks=ranks,singular_values=singular.tolist(),gram_condition=condition,patch_force_differences=patches,constant_shear_difference=constant,rigid_action=re,integration_difference=integration,
                symmetry=symmetry,minimum_eigenvalue=minimum,substitution_energy_difference=energy,substitution_force_difference=force,checks=checks,artifact=file.relative_to(ROOT).as_posix(),artifact_sha256=digest(file)))
            print(json.dumps(dict(case=name,ranks=ranks,checks=checks,patches=patches)),flush=True)
            if time.monotonic()-start>120 or memory.memory_bytes()>512*1024**2:raise RuntimeError('resource_limit')
        assert sum(p.stat().st_size for p in out.iterdir())<=16*1024**2
    except Exception as exc:stop=repr(exc)
    result=dict(source_sha256=digest(__file__),candidate_sha256=digest(ROOT/'docs/q9-trace-shear-candidate.md'),protocol_sha256=digest(F/'q9-trace-shear-gate-experiment.md'),
        dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [candidate,baseline,geometry,memory]},cases=records,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=bool(stop is None and len(records)==3 and all(all(r['checks'].values()) for r in records)),product_shell_adopted=False)
    output.write_text(json.dumps(result,indent=2)+'\n')
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
