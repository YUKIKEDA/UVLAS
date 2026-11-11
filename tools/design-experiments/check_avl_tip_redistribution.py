"""Six stored fields and allocation sensitivities; no Richardson or tail estimate."""
from pathlib import Path
import json,math,time
import numpy as np
import avl_reference_resources as provenance
import avl_symmetry_reconstruction as expansion
import check_avl_conservative_observation as observation
import check_avl_symmetry_equivalence as equivalence
import check_avl_span_order_diagnostic as diagnostic
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification';digest=provenance.digest
def main():
    output=F/'avl-tip-redistribution-result.json'
    if output.exists():raise FileExistsError(output)
    start=time.monotonic();raw=json.loads((F/'avl-tip-redistribution-runs.json').read_text());assert raw['completed']
    baseline=json.loads((F/'avl-span1280-observation-result.json').read_text());assert baseline['methods_criteria_pass']
    models={};observed={};records=[];pairs=[];responses=[];checks={};stop=None;maximum_partition=0.
    def total_bytes():return sum(p.stat().st_size for p in (ROOT/raw['root']).rglob('*') if p.is_file())
    def budget():
        if time.monotonic()-start>120:raise TimeoutError('postprocess_time_limit')
        if memory.memory_bytes()>512*1024**2:raise MemoryError('postprocess_memory_limit')
        if total_bytes()>128*1024**2:raise RuntimeError('combined_output_limit')
    selections=[]
    for case in ['FW-V01','FW-V02']:
        rec=next(r for r in baseline['meshes'] if (r['case'],r['m'],r['n'])==(case,12,1280))
        selections.append(dict(case=case,allocation='uniform',origin='baseline',interval_counts=[64]*10,raw_directory=rec['raw_directory'],derived_directory=rec['derived_directory']))
    selections += [dict(r,origin='new',raw_directory=r['directory']+'/out',derived_directory=r['directory']+'/physical-fullwing') for r in raw['runs']]
    try:
        for r in selections:
            folder=ROOT/r['raw_directory'];derived=ROOT/r['derived_directory'];key=(r['case'],r['allocation'])
            if r['origin']=='new':expansion.expand(folder,12,1280,derived)
            native=np.loadtxt(folder/'uv-rows.txt',skiprows=2);forces=np.loadtxt(derived/'uv-forces.txt',skiprows=2);strips=np.loadtxt(folder/'uv-strips.txt',skiprows=2)
            assert native.shape==(7680,21) and forces.shape==(46080,9) and strips.shape==(640,5)
            assert np.all(np.isfinite(native)) and np.all(np.isfinite(forces)) and np.all(strips[:,3]==12)
            source=folder.parent/'wing.avl';lines=source.read_text().splitlines();sections=[list(map(float,lines[j+1].split()[:5])) for j,line in enumerate(lines) if line=='SECTION']
            assert len(sections)==11
            boundaries=np.array([s[1] for s in sections]);centers=(strips[:,1]+strips[:,2])/2;interval=np.searchsorted(boundaries[1:],centers,side='right')
            assert np.min(interval)==0 and np.max(interval)==9 and np.all(strips[:,2]>strips[:,1])
            counts=np.bincount(interval,minlength=10).tolist();assert counts==r['interval_counts']
            assert np.all(strips[:,1]>=boundaries[interval]-1e-11) and np.all(strips[:,2]<=boundaries[interval+1]+1e-11)
            assert max(abs(strips[0,1]),abs(strips[-1,2]-4),float(np.max(abs(strips[1:,1]-strips[:-1,2]))))<=1e-11
            meta={s.split()[0]:int(s.split()[1]) for s in (folder/'uh-meta.txt').read_text().splitlines()};assert meta['native_unknowns']==7680 and meta['physical_panels']==15360 and meta['real_bits']==64
            models[key],observed[key],record=observation.native_model(dict(case=r['case'],m=12,n=1280,directory=r['derived_directory']),'tip-allocation-'+r['allocation'])
            coeff=observation.prior.old.coefficients(observation.prior.old.fields(folder/'total.txt'));sums=observation.prior.old.aggregate(forces[:,3:6],forces[:,6:9])
            cf=np.array([math.fsum(sums['force'][:,j]) for j in range(3)])/observation.S
            cm=np.array([math.fsum(sums['moment'][:,j]) for j in range(3)])/(observation.S*np.array([8.,observation.C,8.]))
            record.update(allocation=r['allocation'],origin=r['origin'],actual_interval_counts=counts,native_unknowns=7680,physical_panels=15360,
                raw_directory=r['raw_directory'],derived_directory=r['derived_directory'],guard_pass=equivalence.guards(folder),
                residual_relative=float(np.max(abs(native[:,20]))/np.max(np.linalg.norm(native[:,14:17],axis=1))),normal_x_max=float(np.max(abs(native[:,4]))),
                force_sum_difference=float(np.max(abs(cf-coeff['force']))),moment_sum_difference=float(np.max(abs(cm-coeff['moment']))),
                files={p.name:dict(sha256=digest(p),bytes=p.stat().st_size) for p in derived.iterdir()})
            record['eligible']=record['guard_pass'] and record['residual_relative']<=1e-11 and record['normal_x_max']==0 and max(record[n] for n in ['force_sum_difference','moment_sum_difference','load_identity_difference','independent_band_difference','p1_moment_preservation'])<=1e-10
            records.append(record);budget()
        for case in ['FW-V01','FW-V02']:
            vector={a:observation.channels(observed[(case,a)]) for a in ['uniform','tip-rich','tip-poor']};masks=diagnostic.regions(observed[(case,'uniform')])
            for a,b in [('uniform','tip-rich'),('uniform','tip-poor'),('tip-rich','tip-poor')]:
                ka=(case,a);kb=(case,b);item=observation.pair_bounds(models[ka],models[kb],observed[ka],observed[kb]);channels={}
                for c in vector['uniform']:
                    delta=vector[b][c]-vector[a][c];norm=float(np.linalg.norm(delta));energy=delta[:20]**2+(delta[20:]**2 if c=='tau_p1' else 0)
                    regional={r:float(np.sum(energy[mask])) for r,mask in masks.items()};maximum_partition=max(maximum_partition,abs(sum(regional.values())-norm**2))
                    channels[c]=dict(difference=delta.tolist(),norm=norm,baseline_norm=float(np.linalg.norm(vector['uniform'][c])),relative_to_uniform=diagnostic.ratio(norm,np.linalg.norm(vector['uniform'][c])),
                        band_energy=energy.tolist(),regional_energy=regional,regional_fractions={r:diagnostic.ratio(v,norm**2) for r,v in regional.items()})
                scalars={}
                for c in ['lift','moment_line','moment_point']:
                    delta=[float(y[c]-x[c]) for x,y in zip(observed[ka],observed[kb])];s=diagnostic.signed(delta);base=math.fsum(o[c] for o in observed[(case,'uniform')]);den=max(abs(base),.1 if c.startswith('moment') else 0.)
                    scalars[c]=dict(bands=delta,signed=s,baseline=base,relative_to_uniform=abs(s['total'])/den)
                item.update(case=case,from_allocation=a,to_allocation=b,channels=channels,scalar_changes=scalars);pairs.append(item);budget()
            responses.append(dict(case=case,channels={c:dict(direction_cosine=diagnostic.cosine(vector['tip-rich'][c]-vector['uniform'][c],vector['tip-poor'][c]-vector['uniform'][c]),
                rich_response_norm=float(np.linalg.norm(vector['tip-rich'][c]-vector['uniform'][c])),poor_response_norm=float(np.linalg.norm(vector['tip-poor'][c]-vector['uniform'][c]))) for c in vector['uniform']},
                scope='allocation sensitivity only; simultaneous inner and tip changes'))
        checks=dict(six_fields=len(records)==6 and sum(r['origin']=='new' for r in records)==4,equations_and_conservation=all(r['eligible'] for r in records),
            actual_allocation_counts=all(sum(r['actual_interval_counts'])==640 for r in records),six_pairs=len(pairs)==6,
            stored_field_bounds=all(min(p['margins'].values())>=-1e-10 and p['contraction_excess']<=1e-10 for p in pairs),
            regional_partition=maximum_partition<=1e-12,resource_budget=time.monotonic()-start<120 and memory.memory_bytes()<512*1024**2 and total_bytes()<=128*1024**2)
    except Exception as exc:stop=repr(exc);raise
    finally:
        result=dict(protocol='docs/verification/avl-tip-redistribution-experiment.md',protocol_sha256=digest(F/'avl-tip-redistribution-experiment.md'),
            execution_protocol_sha256=digest(F/'avl-tip-redistribution-execution-addendum.md'),source_sha256=digest(__file__),
            dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [provenance,expansion,observation,equivalence,diagnostic,memory]},
            parent_sha256={n:digest(F/n) for n in ['avl-tip-redistribution-runs.json','avl-span1280-observation-result.json']},
            meshes=records,pairs=pairs,responses=responses,checks=checks,methods_criteria_pass=bool(checks) and all(checks.values()),stop_reason=stop,
            maximum_partition_difference=maximum_partition,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),combined_output_bytes=total_bytes(),
            richardson_order=None,continuum_tail_estimate=None,m24_error_budget=None,reference_acceptance=False)
        observation.write(output,result)
    assert all(checks.values()),checks
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],combined_output_bytes=result['combined_output_bytes'],responses=responses,
        relative_changes=[dict(case=p['case'],comparison=[p['from_allocation'],p['to_allocation']],channels={c:d['relative_to_uniform'] for c,d in p['channels'].items()},
            scalars={c:d['relative_to_uniform'] for c,d in p['scalar_changes'].items()}) for p in pairs]),indent=2))
if __name__=='__main__':main()
