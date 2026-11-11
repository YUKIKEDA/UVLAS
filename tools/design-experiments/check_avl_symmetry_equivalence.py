"""Equivalence gate: all native per-panel data are checked before fine symmetry runs."""
from pathlib import Path
import hashlib,json,math,time
import numpy as np
import avl_symmetry_reconstruction as expansion
import check_avl_conservative_observation as observation
import check_q9_plate_convergence as memory

ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):observation.write(p,x)
def guards(folder):
    a=np.loadtxt(folder/'uh-guard-tests.txt',skiprows=1)
    return bool(a.shape==(21,3) and np.all((a[:,1]!=0)==(a[:,2]!=0)))
def main():
    begin=time.monotonic();out=FOLDER/'avl-symmetry-equivalence-result.json'
    if out.exists():raise FileExistsError(out)
    raw_path=FOLDER/'avl-symmetry-regression-runs.json';raw=json.loads(raw_path.read_text());assert raw['completed']
    prep_path=FOLDER/'avl-symmetry-preparation.json';prep=json.loads(prep_path.read_text())
    records=[];matrix=[];stop=None;checks={};synthetic=expansion.synthetic()
    try:
        for run in raw['runs']:
            if run['kind']!='half_measured':continue
            folder=ROOT/run['directory']/'out';full_folder=ROOT/run['comparison_directory']
            derived=ROOT/run['directory']/'physical-fullwing'
            a,forces,info=expansion.expand(folder,run['m'],run['physical_n'],derived)
            b=np.loadtxt(full_folder/'uv-rows.txt',skiprows=2);reference=np.loadtxt(full_folder/'uv-forces.txt',skiprows=2)
            f=expansion.canonical_forces(forces);g=expansion.canonical_forces(reference)
            assert np.array_equal(f[:,:3],g[:,:3]) and a.shape==b.shape
            gamma_diff=abs(a[:,13]-b[:,13]);target=np.maximum(1e-11,1e-10*abs(b[:,13]))
            run_physical=dict(case=run['case'],m=run['m'],n=run['physical_n'],directory=derived.relative_to(ROOT).as_posix())
            model,obs,meta=observation.native_model(run_physical,'symmetric')
            old_model,old_obs,old_meta=observation.native_model(dict(run_physical,directory=run['comparison_directory']),'full')
            bands=np.array([[float(v[k]) for k in ['lift','moment_point']] for v in obs]);old_bands=np.array([[float(v[k]) for k in ['lift','moment_point']] for v in old_obs])
            coefficients=observation.prior.old.coefficients(observation.prior.old.fields(folder/'total.txt'))
            old_coeff=observation.prior.old.coefficients(observation.prior.old.fields(full_folder/'total.txt'))
            global_difference=max(np.max(abs(np.array(coefficients['force'])-old_coeff['force'])),np.max(abs(np.array(coefficients['moment'])-old_coeff['moment'])))
            native=np.loadtxt(folder/'uv-rows.txt',skiprows=2)
            residual=float(np.max(abs(native[:,20]))/np.max(np.linalg.norm(native[:,14:17],axis=1)))
            rec=dict(case=run['case'],m=run['m'],physical_n=run['physical_n'],native_unknowns=len(native),physical_panels=len(a),
                raw_directory=(folder.relative_to(ROOT).as_posix()),derived_directory=derived.relative_to(ROOT).as_posix(),reference_directory=run['comparison_directory'],
                geometry_difference=float(np.max(abs(a[:,1:13]-b[:,1:13]))),point_difference=float(np.max(abs(f[:,3:6]-g[:,3:6]))),
                gamma_difference=float(np.max(gamma_diff)),gamma_eligible=bool(np.all(gamma_diff<=target)),
                velocity_difference=float(max(np.max(abs(a[:,14:17]-b[:,14:17])),np.max(abs(a[:,17:20]-b[:,17:20])))),
                force_difference=float(np.max(abs(f[:,6:9]-g[:,6:9]))),global_coefficient_difference=float(global_difference),
                conservative_band_difference=float(np.max(abs(bands-old_bands))),residual_relative=residual,guard_cases_pass=guards(folder),
                raw_total_semantics='AERO explicitly doubles/zeros totals; equivalence also checks individual forces',
                derivation=info,files={p.name:dict(sha256=digest(p),bytes=p.stat().st_size) for p in derived.iterdir()})
            rec['eligible']=rec['geometry_difference']<=1e-11 and rec['point_difference']<=1e-11 and rec['gamma_eligible'] and rec['velocity_difference']<=1e-10 and max(rec[k] for k in ['force_difference','global_coefficient_difference','conservative_band_difference'])<=1e-10 and residual<=1e-11 and rec['guard_cases_pass']
            records.append(rec)
            if run['m']==6:
                full_run=next(r for r in raw['runs'] if r['kind']=='full_control' and r['case']==run['case'])
                full=ROOT/full_run['directory']/'out';nr=run['physical_n']//2;m=run['m'];h=nr*m
                af=np.fromfile(full/'uh-matrix.bin',dtype='<f8').reshape((2*h,2*h),order='F')
                ah=np.fromfile(folder/'uh-matrix.bin',dtype='<f8').reshape((h,h),order='F')
                positive=np.arange(h,2*h);negative=np.array([(nr-1-s)*m+j for s in range(nr) for j in range(m)])
                bpp=af[np.ix_(positive,positive)];cpn=af[np.ix_(positive,negative)];expected=bpp+cpn
                fr=np.loadtxt(full/'uv-rows.txt',skiprows=2);hr=native
                rhsf=-np.sum(fr[positive,4:7]*fr[positive,14:17],axis=1);rhsh=-np.sum(hr[:,4:7]*hr[:,14:17],axis=1)
                comparison=dict(case=run['case'],full_unknowns=2*h,half_unknowns=h,
                    maximum_matrix_difference=float(np.max(abs(expected-ah))),matrix_eligible=bool(np.all(abs(expected-ah)<=1e-10+1e-10*abs(expected))),
                    rhs_difference=float(np.max(abs(rhsf-rhsh))),block_symmetry_difference=float(max(np.max(abs(bpp-af[np.ix_(negative,negative)])),np.max(abs(cpn-af[np.ix_(negative,positive)])))))
                comparison['eligible']=comparison['matrix_eligible'] and comparison['rhs_difference']<=1e-12
                matrix.append(comparison)
            print(json.dumps({k:rec[k] for k in ['case','m','physical_n','geometry_difference','gamma_difference','force_difference','eligible']}),flush=True)
        checks=dict(native_half_output_invariance=all(all(r.get('parent_half_outputs_equal',{}).values()) and all(r['before_after_equal'].values()) for r in raw['runs']),
            full_controls=len([r for r in raw['runs'] if r['kind']=='full_control'])==2 and all(all(r['previous_outputs_equal'].values()) for r in raw['runs'] if r['kind']=='full_control'),
            six_full_equivalences=len(records)==6 and all(r['eligible'] for r in records),two_matrix_blocks=len(matrix)==2 and all(r['eligible'] for r in matrix),
            vector_reflection=max(synthetic.values())<=1e-12,resource_budget=time.monotonic()-begin<120 and memory.memory_bytes()<512*1024**2)
    except Exception as exc:stop=repr(exc);raise
    finally:
        result=dict(protocol='docs/verification/avl-symmetry-reduction-experiment.md',source_sha256=digest(Path(__file__)),
            dependency_sha256={Path(m.__file__).name:digest(Path(m.__file__)) for m in [expansion,observation,memory]},
            preparation_sha256=digest(prep_path),regression_runs_sha256=digest(raw_path),records=records,matrix=matrix,synthetic=synthetic,
            checks=checks,all_criteria_pass=bool(checks) and all(checks.values()),completed=len(records)==6 and len(matrix)==2 and stop is None,
            stop_reason=stop,reference_convergence_accepted=False,scope='even static two-wing reference only; antisymmetric/dynamic modes unverified',
            seconds=time.monotonic()-begin,peak_process_working_bytes=memory.memory_bytes())
        write(out,result)
    print(json.dumps(dict(checks=checks,seconds=result['seconds']),indent=2))

if __name__=='__main__':main()
