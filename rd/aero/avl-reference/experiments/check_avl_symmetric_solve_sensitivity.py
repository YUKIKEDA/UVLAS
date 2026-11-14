"""Map local solve variations to the same observations and allocation responses."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,math,time
import numpy as np
import avl_symmetric_observation_weights as weights
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;STATES=['direct','ir1','ir2']
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=F/'avl-symmetric-solve-sensitivity-result.json'
    if output.exists():raise FileExistsError(output)
    begin=time.monotonic();raw=json.loads((F/'avl-symmetric-solve-main-runs.json').read_text());assert raw['completed']
    original=json.loads((F/'avl-tip-redistribution-result.json').read_text());assert original['methods_criteria_pass']
    records=[];mapped={};original_mapped={};width_map={};comparisons=[];maximum_reproduction=0.;maximum_subtraction=0.;maximum_quad=0.;maximum_map_response=0.
    for run in raw['runs']:
        path=ROOT/run['weight_path'];assert digest(path)==run['weight_sha256'];cache=np.load(path)
        matrix=cache['matrix'];widths=cache['widths'];saved=cache['original_observations'];folder=ROOT/run['directory']/'uq'
        vectors=np.loadtxt(folder/'uq-vectors.txt',skiprows=1);n=run['native_unknowns'];assert vectors.shape==(n,10) and np.all(np.isfinite(vectors))
        q=vectors[:,2:6];r=vectors[:,6:10];rhs=vectors[:,1]
        baseline=weights.apply(matrix,q[:,0]);error=float(np.max(abs(baseline-saved)/np.maximum(1.,abs(saved))));maximum_reproduction=max(maximum_reproduction,error)
        assert error<=1e-12 and np.max(abs(q[:,0]-cache['q0']))<=1e-12*max(1.,float(np.max(abs(cache['q0']))))
        native=np.loadtxt(ROOT/run['directory']/'out/uv-rows.txt',skiprows=2);assert np.array_equal(rhs,-native[:,6])
        base_channels=weights.channels(baseline,widths);base_totals=weights.totals(baseline,widths);state_records=[]
        key=(run['case'],run['allocation']);mapped[key]={};original_mapped[key]=baseline;width_map[key]=widths
        for j,label in enumerate(STATES,1):
            dq=q[:,j]-q[:,0];change=weights.apply(matrix,dq);separate=weights.apply(matrix,q[:,j])-baseline
            subtraction=float(np.max(abs(change-separate)));maximum_subtraction=max(maximum_subtraction,subtraction)
            channels=weights.channels(change,widths);totals=weights.totals(change,widths);mapped[key][label]=channels
            state_records.append(dict(state=label,circulation_change_norm=float(np.linalg.norm(dq)),circulation_change_max=float(np.max(abs(dq))),
                circulation_relative=float(np.linalg.norm(dq)/np.linalg.norm(q[:,0])),observation_change=change.tolist(),separate_subtraction_difference=subtraction,
                channels={c:dict(norm=float(np.linalg.norm(v)),relative_to_own_original=float(np.linalg.norm(v)/np.linalg.norm(base_channels[c]))) for c,v in channels.items()},
                total_changes={c:dict(absolute=totals[c],relative=abs(totals[c])/max(abs(base_totals[c]),.1 if c.startswith('moment') else 0.)) for c in totals}))
        cross=np.loadtxt(folder/'uq-residual-cross.txt',skiprows=1);assert cross.shape==(32,5)
        assert set(cross[:,0].astype(int))=={1,2,n//4,n//2,n//2+1,3*n//4,n-1,n}
        quad_scaled=float(np.max(abs(cross[:,2]-cross[:,3])/np.maximum(1.,cross[:,4])));maximum_quad=max(maximum_quad,quad_scaled)
        for row in cross:assert row[3]==r[int(row[0])-1,int(row[1])-1]
        records.append(dict(case=run['case'],allocation=run['allocation'],native_unknowns=n,original_observation_reproduction=error,
            residual_infinity=[float(np.max(abs(r[:,j]))) for j in range(4)],sampled_quad_compensated_difference=quad_scaled,
            sampled_rows=sorted(set(cross[:,0].astype(int).tolist())),states=state_records,original_observations=baseline.tolist(),widths=widths.tolist()))
        if time.monotonic()-begin>300 or memory.memory_bytes()>512*1024**2:raise RuntimeError('postprocess_resource_limit')
    for pair in original['pairs']:
        if pair['from_allocation']!='uniform':continue
        case=pair['case'];allocation=pair['to_allocation'];baseline_key=(case,'uniform');target_key=(case,allocation)
        original_response=weights.channels(original_mapped[target_key]-original_mapped[baseline_key],width_map[baseline_key])
        reproduction={}
        for c,signal in pair['channels'].items():
            difference=float(np.linalg.norm(original_response[c]-signal['difference']));maximum_map_response=max(maximum_map_response,difference)
            reproduction[c]=dict(norm_difference=difference,relative_to_response=difference/signal['norm'] if signal['norm'] else None)
        rec=[]
        for label in STATES:
            channels={}
            for c,signal in pair['channels'].items():
                du=mapped[baseline_key][label][c];dt=mapped[target_key][label][c];change=dt-du;scale=signal['norm']
                channels[c]=dict(pair_response_change_norm=float(np.linalg.norm(change)),individual_change_norm_sum=float(np.linalg.norm(du)+np.linalg.norm(dt)),
                    original_response_norm=scale,response_change_ratio=float(np.linalg.norm(change)/scale) if scale else None,
                    individual_sum_ratio=float((np.linalg.norm(du)+np.linalg.norm(dt))/scale) if scale else None,
                    scope='sensitivity ratio; not an error bound')
            rec.append(dict(state=label,channels=channels))
        comparisons.append(dict(case=case,allocation=allocation,states=rec,original_response_map_difference=reproduction))
    checks=dict(six_systems=len(records)==6,eighteen_state_changes=sum(len(r['states']) for r in records)==18,
        original_observation_reproduction=maximum_reproduction<=1e-12,delta_map_subtraction_consistency=maximum_subtraction<=1e-12,
        sampled_quad_residuals=maximum_quad<=1e-12,four_allocation_comparisons=len(comparisons)==4,allocation_response_map_consistency=maximum_map_response<=1e-12,
        resource_budget=time.monotonic()-begin<300 and memory.memory_bytes()<512*1024**2)
    result=dict(protocol='rd/aero/avl-reference/records/avl-symmetric-solve-sensitivity-experiment.md',protocol_sha256=digest(F/'avl-symmetric-solve-sensitivity-experiment.md'),
        source_sha256=digest(__file__),dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [weights,memory]},
        parent_sha256={n:digest(F/n) for n in ['avl-symmetric-solve-main-runs.json','avl-symmetric-solve-small-check.json','avl-tip-redistribution-result.json','avl-symmetric-solve-weight-result.json']},
        systems=records,comparisons=comparisons,checks=checks,all_criteria_pass=all(checks.values()),seconds=time.monotonic()-begin,peak_bytes=memory.memory_bytes(),
        maximum_original_reproduction=maximum_reproduction,maximum_subtraction_difference=maximum_subtraction,maximum_sampled_quad_difference=maximum_quad,maximum_original_response_map_difference=maximum_map_response,
        forward_error_bound=None,coefficient_generation_error_bound=None,m24_error_budget=None,reference_acceptance=False)
    output.write_text(json.dumps(result,indent=2)+'\n');assert output.stat().st_size<16*1024**2 and result['all_criteria_pass']
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],allocation_ratios=[dict(case=c['case'],allocation=c['allocation'],states={s['state']:{k:v['response_change_ratio'] for k,v in s['channels'].items()} for s in c['states']}) for c in comparisons]),indent=2))
if __name__=='__main__':main()
