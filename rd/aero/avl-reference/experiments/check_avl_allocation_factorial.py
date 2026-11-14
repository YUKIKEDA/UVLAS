"""Eight stored fields and independent allocation effects; no Richardson or tail estimate."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,math,time
import numpy as np
import avl_reference_resources as provenance
import avl_symmetry_reconstruction as expansion
import check_avl_conservative_observation as observation
import check_avl_symmetry_equivalence as equivalence
import check_avl_span_order_diagnostic as diagnostic
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;digest=provenance.digest
def main():
    output=F/'avl-allocation-factorial-result.json'
    if output.exists():raise FileExistsError(output)
    start=time.monotonic();raw=json.loads((F/'avl-allocation-factorial-runs.json').read_text());assert raw['completed']
    baseline=json.loads((F/'avl-span1280-observation-result.json').read_text());assert baseline['methods_criteria_pass']
    models={};observed={};records=[];pairs=[];responses=[];checks={};stop=None;max_reconstruction=0.;max_energy_partition=0.
    def total_bytes():return sum(p.stat().st_size for p in (ROOT/raw['root']).rglob('*') if p.is_file())
    def budget():
        if time.monotonic()-start>180:raise TimeoutError('postprocess_time_limit')
        if memory.memory_bytes()>512*1024**2:raise MemoryError('postprocess_memory_limit')
        if total_bytes()>128*1024**2:raise RuntimeError('combined_output_limit')
    selections=[]
    for case in ['FW-V01','FW-V02']:
        rec=next(r for r in baseline['meshes'] if (r['case'],r['m'],r['n'])==(case,12,1280))
        selections.append(dict(case=case,allocation='uniform',origin='baseline',interval_counts=[64]*10,raw_directory=rec['raw_directory'],derived_directory=rec['derived_directory']))
    tip=json.loads((F/'avl-tip-redistribution-result.json').read_text());assert tip['methods_criteria_pass']
    selections += [dict(r,origin='existing-tip') for r in tip['meshes'] if r['allocation']=='tip-rich']
    selections += [dict(r,origin='new',raw_directory=r['directory']+'/out',derived_directory=r['directory']+'/physical-fullwing') for r in raw['runs']]
    try:
        for r in selections:
            folder=ROOT/r['raw_directory'];derived=ROOT/r['derived_directory'];key=(r['case'],r['allocation'])
            n=2*sum(r.get('interval_counts',r.get('actual_interval_counts',[64]*10)));size=12*n//2
            if r['origin']=='new':expansion.expand(folder,12,n,derived)
            native=np.loadtxt(folder/'uv-rows.txt',skiprows=2);forces=np.loadtxt(derived/'uv-forces.txt',skiprows=2);strips=np.loadtxt(folder/'uv-strips.txt',skiprows=2)
            assert native.shape==(size,21) and forces.shape==(6*size,9) and strips.shape==(n//2,5)
            assert np.all(np.isfinite(native)) and np.all(np.isfinite(forces)) and np.all(strips[:,3]==12)
            source=folder.parent/'wing.avl';lines=source.read_text().splitlines();sections=[list(map(float,lines[j+1].split()[:5])) for j,line in enumerate(lines) if line=='SECTION']
            assert len(sections)==11
            boundaries=np.array([s[1] for s in sections]);centers=(strips[:,1]+strips[:,2])/2;interval=np.searchsorted(boundaries[1:],centers,side='right')
            assert np.min(interval)==0 and np.max(interval)==9 and np.all(strips[:,2]>strips[:,1])
            counts=np.bincount(interval,minlength=10).tolist();assert counts==r.get('interval_counts',r.get('actual_interval_counts'))
            assert np.all(strips[:,1]>=boundaries[interval]-1e-11) and np.all(strips[:,2]<=boundaries[interval+1]+1e-11)
            assert max(abs(strips[0,1]),abs(strips[-1,2]-4),float(np.max(abs(strips[1:,1]-strips[:-1,2]))))<=1e-11
            meta={s.split()[0]:int(s.split()[1]) for s in (folder/'uh-meta.txt').read_text().splitlines()};assert meta['native_unknowns']==size and meta['physical_panels']==2*size and meta['real_bits']==64
            models[key],observed[key],record=observation.native_model(dict(case=r['case'],m=12,n=n,directory=r['derived_directory']),'tip-allocation-'+r['allocation'])
            coeff=observation.prior.old.coefficients(observation.prior.old.fields(folder/'total.txt'));sums=observation.prior.old.aggregate(forces[:,3:6],forces[:,6:9])
            cf=np.array([math.fsum(sums['force'][:,j]) for j in range(3)])/observation.S
            cm=np.array([math.fsum(sums['moment'][:,j]) for j in range(3)])/(observation.S*np.array([8.,observation.C,8.]))
            record.update(allocation=r['allocation'],origin=r['origin'],actual_interval_counts=counts,native_unknowns=size,physical_panels=2*size,
                raw_directory=r['raw_directory'],derived_directory=r['derived_directory'],guard_pass=equivalence.guards(folder),
                residual_relative=float(np.max(abs(native[:,20]))/np.max(np.linalg.norm(native[:,14:17],axis=1))),normal_x_max=float(np.max(abs(native[:,4]))),
                force_sum_difference=float(np.max(abs(cf-coeff['force']))),moment_sum_difference=float(np.max(abs(cm-coeff['moment']))),
                files={p.name:dict(sha256=digest(p),bytes=p.stat().st_size) for p in derived.iterdir()})
            record['eligible']=record['guard_pass'] and record['residual_relative']<=1e-11 and record['normal_x_max']==0 and max(record[n] for n in ['force_sum_difference','moment_sum_difference','load_identity_difference','independent_band_difference','p1_moment_preservation'])<=1e-10
            records.append(record);budget()
        max_reconstruction=0.;max_energy_partition=0.
        for case in ['FW-V01','FW-V02']:
            labels=['uniform','inner-only','tip-only','tip-rich']
            vectors={a:observation.channels(observed[(case,a)]) for a in labels}
            masks=diagnostic.regions(observed[(case,'uniform')]);channels={}
            for channel in vectors['uniform']:
                A,B,C,D=[vectors[a][channel] for a in labels]
                effects={'inner':B-A,'tip':C-A,'interaction':D-B-C+A,'combined':D-A}
                error=float(np.max(abs(effects['inner']+effects['tip']+effects['interaction']-effects['combined'])))
                scale=max(1.,*(float(np.max(abs(v))) for v in [A,B,C,D]))
                max_reconstruction=max(max_reconstruction,error/scale)
                info={}
                for name,delta in effects.items():
                    energy=delta[:20]**2+(delta[20:]**2 if channel=='tau_p1' else 0)
                    regional={r:float(np.sum(energy[mask])) for r,mask in masks.items()}
                    norm=float(np.linalg.norm(delta));max_energy_partition=max(max_energy_partition,abs(sum(regional.values())-norm**2))
                    info[name]=dict(difference=delta.tolist(),norm=norm,relative_to_uniform=diagnostic.ratio(norm,np.linalg.norm(A)),
                        band_energy=energy.tolist(),regional_energy=regional,
                        cosine_to_combined=diagnostic.cosine(delta,effects['combined']))
                cross={x+'_'+y:2*float(np.dot(effects[x],effects[y])) for x,y in [('inner','tip'),('inner','interaction'),('tip','interaction')]}
                channels[channel]=dict(effects=info,cross_terms=cross,reconstruction_error=error)
            scalars={}
            for channel in ['lift','moment_line','moment_point']:
                values={a:np.array([float(o[channel]) for o in observed[(case,a)]]) for a in labels}
                A,B,C,D=[values[a] for a in labels]
                effects={'inner':B-A,'tip':C-A,'interaction':D-B-C+A,'combined':D-A}
                err=float(np.max(abs(effects['inner']+effects['tip']+effects['interaction']-effects['combined'])))
                max_reconstruction=max(max_reconstruction,err/max(1.,*(float(np.max(abs(v))) for v in values.values())))
                scalars[channel]=dict(values={a:math.fsum(v) for a,v in values.items()},effects={name:dict(bands=v.tolist(),signed=diagnostic.signed(v)) for name,v in effects.items()})
            responses.append(dict(case=case,channels=channels,scalars=scalars));budget()
        checks=dict(eight_fields=len(records)==8 and sum(r['origin']=='new' for r in records)==4,
            equations_and_conservation=all(r['eligible'] for r in records),
            actual_allocation_counts=all(r['actual_interval_counts']==({'uniform':[64]*10,'inner-only':[48]*8+[64]*2,'tip-only':[64]*8+[128]*2,'tip-rich':[48]*8+[128]*2}[r['allocation']]) for r in records),
            factorial_reconstruction=max_reconstruction<=1e-12,regional_partition=max_energy_partition<=1e-12,
            resource_budget=time.monotonic()-start<180 and memory.memory_bytes()<512*1024**2 and total_bytes()<=128*1024**2)
    except Exception as exc:stop=repr(exc);raise
    finally:
        result=dict(protocol='rd/aero/avl-reference/records/avl-allocation-factorial-experiment.md',protocol_sha256=digest(F/'avl-allocation-factorial-experiment.md'),
            source_sha256=digest(__file__),
            dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [provenance,expansion,observation,equivalence,diagnostic,memory]},
            parent_sha256={n:digest(F/n) for n in ['avl-allocation-factorial-runs.json','avl-span1280-observation-result.json','avl-tip-redistribution-result.json']},
            meshes=records,pairs=pairs,responses=responses,checks=checks,methods_criteria_pass=bool(checks) and all(checks.values()),stop_reason=stop,
            maximum_partition_difference=max_energy_partition,maximum_reconstruction_difference=max_reconstruction,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),combined_output_bytes=total_bytes(),
            richardson_order=None,continuum_tail_estimate=None,m24_error_budget=None,reference_acceptance=False)
        observation.write(output,result)
    assert all(checks.values()),checks
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],peak_bytes=result['peak_bytes'],combined_output_bytes=result['combined_output_bytes']),indent=2))
if __name__=='__main__':main()
