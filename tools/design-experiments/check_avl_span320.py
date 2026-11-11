"""Audit the added span level and retain conditional error-budget semantics."""
from pathlib import Path
import hashlib,json,math,time
import numpy as np
import check_avl_conservative_observation as observation
import check_q9_plate_convergence as memory

ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    begin=time.monotonic();path=FOLDER/'avl-span320-runs.json';raw=json.loads(path.read_text());assert raw['completed']
    oldruns=json.loads((FOLDER/'avl-band-aligned-grid-runs.json').read_text())
    parent_path=FOLDER/'avl-conservative-observation-result.json';parent=json.loads(parent_path.read_text())
    new=[];series=[];budgets=[];pairs=[]
    for case in ['FW-V01','FW-V02']:
        models=[];observed=[];metadata=[]
        for n in [40,80,160,320]:
            if n==320:
                r=next(r for r in raw['runs'] if r['kind']=='new' and r['case']==case)
                run=dict(r,directory=r['directory']+'/out')
            else:run=next(r for r in oldruns['runs'] if (r['case'],r['m'],r['n'],r['alpha'])==(case,24,n,1))
            model,obs,record=observation.native_model(run,'aligned');models.append(model);observed.append(obs);metadata.append(record)
            if n==320:
                folder=ROOT/run['directory'];rows=np.loadtxt(folder/'uv-rows.txt',skiprows=2);forces=np.loadtxt(folder/'uv-forces.txt',skiprows=2)
                assert rows.shape==(7680,21) and forces.shape==(23040,9) and np.all(np.isfinite(rows)) and np.all(np.isfinite(forces))
                coefs=observation.prior.old.coefficients(observation.prior.old.fields(folder/'total.txt'))
                aggregated=observation.prior.old.aggregate(forces[:,3:6],forces[:,6:9])
                force=np.array([math.fsum(aggregated['force'][:,j]) for j in range(3)])/observation.S
                moment=np.array([math.fsum(aggregated['moment'][:,j]) for j in range(3)])/(observation.S*np.array([8.,observation.C,8.]))
                record.update(residual_relative=float(np.max(abs(rows[:,20]))/np.max(np.linalg.norm(rows[:,14:17],axis=1))),
                    force_sum_difference=float(np.max(abs(force-coefs['force']))),moment_sum_difference=float(np.max(abs(moment-coefs['moment']))))
                new.append(record)
        vectors=[observation.channels(o) for o in observed]
        channels={k:observation.norm_diagnostic([v[k] for v in vectors]) for k in vectors[0]}
        series.append(dict(case=case,m=24,n=[40,80,160,320],channels=channels))
        for j in range(3):
            p=observation.pair_bounds(models[j],models[j+1],observed[j],observed[j+1]);p.update(case=case,pair=j);pairs.append(p)
        chord=next(s for s in parent['series'] if (s['family'],s['case'],s['direction'])==('aligned',case,'chord'))
        channel_budgets={}
        for k,v in channels.items():
            ch=chord['channels'][k];ce=ch['conditional_extrapolation'];se=v['conditional_extrapolation']
            available=ce is not None and se is not None and v['orders'] is not None
            transferred=se*2**v['orders'][-1] if available else None
            combined=ce+transferred if available else None
            channel_budgets[k]=dict(chord_status=ch['status'],span_status=v['status'],chord_base='M24N160',span_base='M24N320',
                chord_conditional=ce,span_conditional_at320=se,span_conditional_transferred_to160=transferred,
                combined_conditional_at160=combined,relative_at160=combined/ch['scale'] if combined is not None and ch['scale']>0 else None,
                endpoint_transfer='2^p; assumes same span leading error expansion',
                state='conditional_estimate_with_unverified_mixed_term' if available else 'unassessed_missing_axis')
        baseline=metadata[2]
        tau=channel_budgets['tau_p1']['combined_conditional_at160'];jp=channel_budgets['j_point']['combined_conditional_at160']
        distribution=channel_budgets['tau0']['relative_at160']
        lift_relative=2*math.sqrt(8)/observation.S*tau/abs(baseline['lift']) if tau is not None else None
        moment_relative=2*math.sqrt(8)/(observation.S*observation.C)*jp/max(abs(baseline['moment_point']),.1) if jp is not None else None
        budgets.append(dict(case=case,base='M24N160',channels=channel_budgets,
            conditional_lift_relative=lift_relative,conditional_moment_relative=moment_relative,conditional_distribution_relative=distribution,
            conditional_targets_met=all(v is not None and v<=tol for v,tol in [(lift_relative,.003),(moment_relative,.003),(distribution,.005)]),
            unverified=['mixed chord/span term and endpoint transfer assumptions','continuum error remains unbounded'],reference_acceptance_eligible=False))
    checks=dict(capacity_regression=len(raw['runs'])==6 and all(all(r['before_after_equal'].values()) for r in raw['runs']) and all(all(r['previous_outputs_equal'].values()) for r in raw['runs'] if r['kind']=='regression'),
        two_new_cases=len(new)==2 and all(r['residual_relative']<=1e-11 and max(r['force_sum_difference'],r['moment_sum_difference'])<=1e-10 for r in new),
        conservative_observations=all(max(r['load_identity_difference'],r['independent_band_difference'],r['p1_moment_preservation'])<=1e-10 for r in new),
        two_field_bounds=all(min(p['margins'].values())>=-1e-10 and p['contraction_excess']<=1e-10 for p in pairs),
        resource_budget=time.monotonic()-begin<120 and memory.memory_bytes()<512*1024**2)
    result=dict(protocol='docs/verification/avl-span320-reference-experiment.md',source_sha256=digest(Path(__file__)),
        dependency_sha256={Path(m.__file__).name:digest(Path(m.__file__)) for m in [observation,memory]},
        parent_sha256={name:digest(FOLDER/name) for name in ['avl-span320-runs.json','avl-conservative-observation-result.json','avl-band-aligned-grid-runs.json']},
        new_cases=new,series=series,pairs=pairs,budgets=budgets,checks=checks,methods_criteria_pass=all(checks.values()),
        reference_convergence_accepted=False,continuum_error_bound_available=False,seconds=time.monotonic()-begin,peak_process_working_bytes=memory.memory_bytes())
    observation.write(FOLDER/'avl-span320-result.json',result)
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],budgets=[{k:b[k] for k in ['case','conditional_lift_relative','conditional_moment_relative','conditional_distribution_relative','conditional_targets_met']} for b in budgets]),indent=2))

if __name__=='__main__':main()
