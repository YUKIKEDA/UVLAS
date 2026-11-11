"""Same-capacity N1280 span/chord observations with unchanged asymptotic criteria."""
from pathlib import Path
import hashlib,json,math,time
import numpy as np
import avl_symmetry_reconstruction as expansion
import check_avl_conservative_observation as observation
import check_avl_symmetry_equivalence as equivalence
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=F/'avl-span1280-observation-result.json'
    if output.exists():raise FileExistsError(output)
    begin=time.monotonic();raw=json.loads((F/'avl-span1280-fine-ks-runs.json').read_text());assert raw['completed']
    coarse=json.loads((F/'avl-span1280-capacity-ks-runs.json').read_text());assert coarse['completed']
    gate=json.loads((F/'avl-span1280-capacity-ks-artifact-audit.json').read_text());assert gate['all_criteria_pass']
    runs=[r for r in coarse['runs'] if r['m']==24]+raw['runs'];assert len(runs)==14
    models={};observed={};meshes=[];series=[];pairs=[];budgets=[];checks={};stop=None
    def budget():
        if time.monotonic()-begin>120:raise TimeoutError('postprocess_time_limit')
        if memory.memory_bytes()>1024**3:raise MemoryError('postprocess_memory_limit')
    try:
        for run in runs:
            folder=ROOT/run['directory']/'out';derived=ROOT/run['directory']/'physical-fullwing'
            rows,forces,info=expansion.expand(folder,run['m'],run['physical_n'],derived)
            n=run['m']*run['physical_n'];native=np.loadtxt(folder/'uv-rows.txt',skiprows=2)
            assert rows.shape==(n,21) and forces.shape==(3*n,9) and len(native)*2==n
            assert np.all(np.isfinite(rows)) and np.all(np.isfinite(forces))
            key=(run['case'],run['m'],run['physical_n'])
            models[key],observed[key],record=observation.native_model(dict(case=key[0],m=key[1],n=key[2],directory=derived.relative_to(ROOT).as_posix()),'spacing-pi64-ks1000')
            coeff=observation.prior.old.coefficients(observation.prior.old.fields(folder/'total.txt'))
            sums=observation.prior.old.aggregate(forces[:,3:6],forces[:,6:9])
            cf=np.array([math.fsum(sums['force'][:,j]) for j in range(3)])/observation.S
            cm=np.array([math.fsum(sums['moment'][:,j]) for j in range(3)])/(observation.S*np.array([8.,observation.C,8.]))
            meta={s.split()[0]:int(s.split()[1]) for s in (folder/'uh-meta.txt').read_text().splitlines()}
            assert meta['native_unknowns']==len(native)==run['native_unknowns'] and meta['physical_panels']==n and meta['real_bits']==64
            record.update(native_unknowns=len(native),physical_panels=n,raw_directory=folder.relative_to(ROOT).as_posix(),derived_directory=derived.relative_to(ROOT).as_posix(),
                origin='fine' if run['physical_n']==1280 else 'capacity-regression',derivation=info,guard_pass=equivalence.guards(folder),
                residual_relative=float(np.max(abs(native[:,20]))/np.max(np.linalg.norm(native[:,14:17],axis=1))),normal_x_max=float(np.max(abs(native[:,4]))),
                force_sum_difference=float(np.max(abs(cf-coeff['force']))),moment_sum_difference=float(np.max(abs(cm-coeff['moment']))),
                files={p.name:dict(sha256=digest(p),bytes=p.stat().st_size) for p in derived.iterdir()})
            record['eligible']=record['guard_pass'] and record['residual_relative']<=1e-11 and record['normal_x_max']==0 and max(record[k] for k in ['force_sum_difference','moment_sum_difference','load_identity_difference','independent_band_difference','p1_moment_preservation'])<=1e-10
            meshes.append(record)
            budget()
        for case in ['FW-V01','FW-V02']:
            for label,indices in [('span1280',[(24,n) for n in [160,320,640,1280]]),('chord1280',[(m,1280) for m in [3,6,12,24]])]:
                keys=[(case,m,n) for m,n in indices];vectors=[observation.channels(observed[k]) for k in keys]
                series.append(dict(case=case,label=label,meshes=indices,channels={c:observation.norm_diagnostic([v[c] for v in vectors]) for c in vectors[0]}))
                for j in range(3):
                    pair=observation.pair_bounds(models[keys[j]],models[keys[j+1]],observed[keys[j]],observed[keys[j+1]])
                    pair.update(case=case,label=label,pair=j);pairs.append(pair)
                    budget()
            span=next(s for s in series if s['case']==case and s['label']=='span1280');chord=next(s for s in series if s['case']==case and s['label']=='chord1280');channels={}
            for c,s in span['channels'].items():
                ce=chord['channels'][c]['conditional_extrapolation'];se=s['conditional_extrapolation']
                combined=ce+se if ce is not None and se is not None else None;scale=chord['channels'][c]['scale']
                channels[c]=dict(chord_status=chord['channels'][c]['status'],span_status=s['status'],chord_conditional=ce,span_conditional=se,
                    combined_conditional=combined,relative=combined/scale if combined is not None and scale>0 else None,
                    state='conditional_estimate_only' if combined is not None else 'unassessed_missing_axis')
            baseline=next(r for r in meshes if (r['case'],r['m'],r['n'])==(case,24,1280));tau=channels['tau_p1']['combined_conditional'];moment=channels['j_point']['combined_conditional']
            lift=2*math.sqrt(8)/observation.S*tau/abs(baseline['lift']) if tau is not None else None
            cm=2*math.sqrt(8)/(observation.S*observation.C)*moment/max(abs(baseline['moment_point']),.1) if moment is not None else None
            distribution=channels['tau0']['relative']
            budgets.append(dict(case=case,base='M24N1280',same_base_without_transfer=True,channels=channels,
                conditional_lift_relative=lift,conditional_moment_relative=cm,conditional_distribution_relative=distribution,
                conditional_targets_met=all(v is not None and v<=t for v,t in [(lift,.003),(cm,.003),(distribution,.005)]),
                reference_acceptance_eligible=False,unverified=['tail extrapolation assumptions','mixed-axis tail not bounded','continuum error not certified']))
        checks=dict(corrected_capacity_gate=gate['all_criteria_pass'],fourteen_fields=len(meshes)==14 and sum(r['origin']=='fine' for r in meshes)==8,
            equations_and_conservation=all(r['eligible'] for r in meshes),four_series=len(series)==4,
            two_stored_field_bounds=len(pairs)==12 and all(min(p['margins'].values())>=-1e-10 and p['contraction_excess']<=1e-10 for p in pairs),
            resource_budget=time.monotonic()-begin<120 and memory.memory_bytes()<1024**3)
    except Exception as exc:stop=repr(exc);raise
    finally:
        result=dict(protocol='docs/verification/avl-span1280-fine-experiment.md',protocol_sha256=digest(F/'avl-span1280-fine-experiment.md'),
            source_sha256=digest(__file__),dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [expansion,observation,equivalence,memory]},
            parent_sha256={n:digest(F/n) for n in ['avl-span1280-fine-ks-runs.json','avl-span1280-capacity-ks-runs.json','avl-span1280-capacity-ks-artifact-audit.json']},
            meshes=meshes,series=series,pairs=pairs,budgets=budgets,checks=checks,methods_criteria_pass=bool(checks) and all(checks.values()),
            stop_reason=stop,reference_convergence_accepted=False,continuum_error_bound_available=False,seconds=time.monotonic()-begin,peak_process_working_bytes=memory.memory_bytes())
        observation.write(output,result)
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],budgets=[{k:v for k,v in b.items() if k!='channels'} for b in budgets]),indent=2))
if __name__=='__main__':main()
