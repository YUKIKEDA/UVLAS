"""Fine even-static reference observations, with explicit conditional-budget scope."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,math,time
import numpy as np
import avl_symmetry_reconstruction as expansion
import check_avl_conservative_observation as observation
import check_avl_symmetry_equivalence as equivalence
import check_q9_plate_convergence as memory

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):observation.write(p,x)
def main():
    begin=time.monotonic();output=FOLDER/'avl-symmetry-fine-observation-result.json'
    if output.exists():raise FileExistsError(output)
    raw_path=FOLDER/'avl-symmetry-pi-fine-runs.json';raw=json.loads(raw_path.read_text());assert raw['completed']
    gate_path=FOLDER/'avl-symmetry-pi-equivalence-result.json';gate=json.loads(gate_path.read_text());assert gate['all_criteria_pass']
    spacing=json.loads((FOLDER/'avl-spacing-pi-result.json').read_text())
    models={};observed={};meshes=[];series=[];pairs=[];budgets=[];stop=None;checks={}
    try:
        for rec in gate['records']:
            if rec['m']!=24:continue
            run=dict(case=rec['case'],m=24,n=rec['physical_n'],directory=rec['derived_directory'])
            key=(rec['case'],24,rec['physical_n']);models[key],observed[key],record=observation.native_model(run,'spacing-pi64')
            record['origin']='validated-equivalence';meshes.append(record)
        for run in raw['runs']:
            folder=ROOT/run['directory']/'out';derived=ROOT/run['directory']/'physical-fullwing'
            rows,forces,info=expansion.expand(folder,run['m'],run['physical_n'],derived)
            n=run['m']*run['physical_n'];native=np.loadtxt(folder/'uv-rows.txt',skiprows=2)
            assert rows.shape==(n,21) and forces.shape==(3*n,9) and len(native)*2==n
            assert np.all(np.isfinite(rows)) and np.all(np.isfinite(forces))
            key=(run['case'],run['m'],run['physical_n'])
            models[key],observed[key],record=observation.native_model(dict(case=key[0],m=key[1],n=key[2],directory=derived.relative_to(ROOT).as_posix()),'spacing-pi64')
            coeff=observation.prior.old.coefficients(observation.prior.old.fields(folder/'total.txt'))
            sums=observation.prior.old.aggregate(forces[:,3:6],forces[:,6:9])
            cf=np.array([math.fsum(sums['force'][:,j]) for j in range(3)])/observation.S
            cm=np.array([math.fsum(sums['moment'][:,j]) for j in range(3)])/(observation.S*np.array([8.,observation.C,8.]))
            record.update(origin='fine-stage',native_unknowns=len(native),physical_panels=n,derived_directory=derived.relative_to(ROOT).as_posix(),
                raw_directory=folder.relative_to(ROOT).as_posix(),derivation=info,guard_pass=equivalence.guards(folder),
                residual_relative=float(np.max(abs(native[:,20]))/np.max(np.linalg.norm(native[:,14:17],axis=1))),
                normal_x_max=float(np.max(abs(native[:,4]))),force_sum_difference=float(np.max(abs(cf-coeff['force']))),
                moment_sum_difference=float(np.max(abs(cm-coeff['moment']))),
                files={p.name:dict(sha256=digest(p),bytes=p.stat().st_size) for p in derived.iterdir()})
            legacy=None
            if run['physical_n']<640:
                old=ROOT/run['comparison_directory'];a=np.loadtxt(old/'uv-rows.txt',skiprows=2);f=expansion.canonical_forces(forces)
                g=expansion.canonical_forces(np.loadtxt(old/'uv-forces.txt',skiprows=2))
                legacy=dict(geometry_difference=float(np.max(abs(rows[:,1:13]-a[:,1:13]))),point_difference=float(np.max(abs(f[:,3:6]-g[:,3:6]))),
                    gamma_difference=float(np.max(abs(rows[:,13]-a[:,13]))),velocity_difference=float(np.max(abs(rows[:,14:20]-a[:,14:20]))),
                    force_difference=float(np.max(abs(f[:,6:9]-g[:,6:9]))),reference_directory=run['comparison_directory'])
                legacy['eligible']=max(legacy['geometry_difference'],legacy['point_difference'])<=spacing['old_new_coordinate_envelope'] and legacy['gamma_difference']<=1e-11 and max(legacy['force_difference'],legacy['velocity_difference'])<=1e-10
            record['legacy_sensitivity']=legacy
            record['eligible']=record['guard_pass'] and record['residual_relative']<=1e-11 and record['normal_x_max']==0 and max(record[k] for k in ['force_sum_difference','moment_sum_difference','load_identity_difference','independent_band_difference','p1_moment_preservation'])<=1e-10 and (legacy is None or legacy['eligible'])
            meshes.append(record)
        for case in ['FW-V01','FW-V02']:
            definitions=[('span640',[(24,n) for n in [80,160,320,640]]),('chord160',[(m,160) for m in [3,6,12,24]]),('chord640',[(m,640) for m in [3,6,12,24]])]
            for label,indices in definitions:
                keys=[(case,m,n) for m,n in indices];vectors=[observation.channels(observed[k]) for k in keys]
                series.append(dict(case=case,label=label,meshes=indices,channels={c:observation.norm_diagnostic([v[c] for v in vectors]) for c in vectors[0]}))
                for j in range(3):
                    b=observation.pair_bounds(models[keys[j]],models[keys[j+1]],observed[keys[j]],observed[keys[j+1]])
                    b.update(case=case,label=label,pair=j);pairs.append(b)
            span=next(s for s in series if s['case']==case and s['label']=='span640')
            for base,label in [(160,'chord160'),(640,'chord640')]:
                chord=next(s for s in series if s['case']==case and s['label']==label);channel_budgets={}
                for c,d in span['channels'].items():
                    ce=chord['channels'][c]['conditional_extrapolation'];se=d['conditional_extrapolation']
                    power=None;triangle=None;transferred=None
                    if se is not None:
                        if base==640:transferred=se;power=se;triangle=se
                        elif d['orders'] is not None:
                            v0=observation.channels(observed[(case,24,160)])[c];v1=observation.channels(observed[(case,24,640)])[c]
                            power=se*4**d['orders'][-1];triangle=se+float(np.linalg.norm(v1-v0));transferred=max(power,triangle)
                    available=ce is not None and transferred is not None;combined=ce+transferred if available else None
                    scale=chord['channels'][c]['scale']
                    channel_budgets[c]=dict(chord_status=chord['channels'][c]['status'],span_status=d['status'],
                        chord_conditional=ce,span_tail_at640=se,span_power_transferred=power,span_triangle_transferred=triangle,
                        span_selected=transferred,combined_conditional=combined,relative=combined/scale if combined is not None and scale>0 else None,
                        state='conditional_estimate_only' if available else 'unassessed_missing_axis')
                baseline=next(r for r in meshes if (r['case'],r['m'],r['n'])==(case,24,base))
                tau=channel_budgets['tau_p1']['combined_conditional'];moment=channel_budgets['j_point']['combined_conditional']
                lift_rel=2*math.sqrt(8)/observation.S*tau/abs(baseline['lift']) if tau is not None else None
                moment_rel=2*math.sqrt(8)/(observation.S*observation.C)*moment/max(abs(baseline['moment_point']),.1) if moment is not None else None
                distribution=channel_budgets['tau0']['relative']
                # Same observables at both span resolutions quantify, but do not bound, mixed-axis effects.
                mixed={}
                for c in span['channels']:
                    values={n:[observation.channels(observed[(case,m,n)])[c] for m in [3,6,12,24]] for n in [160,640]}
                    mixed[c]=[dict(norm=float(np.linalg.norm((values[640][j+1]-values[640][j])-(values[160][j+1]-values[160][j]))),
                        relative_to_fine_chord=float(np.linalg.norm((values[640][j+1]-values[640][j])-(values[160][j+1]-values[160][j]))/max(np.linalg.norm(values[640][j+1]-values[640][j]),1e-12))) for j in range(3)]
                budgets.append(dict(case=case,base=f'M24N{base}',same_base_without_transfer=base==640,channels=channel_budgets,
                    conditional_lift_relative=lift_rel,conditional_moment_relative=moment_rel,conditional_distribution_relative=distribution,
                    conditional_targets_met=all(v is not None and v<=t for v,t in [(lift_rel,.003),(moment_rel,.003),(distribution,.005)]),
                    mixed_diagnostics=mixed,unverified=['tail extrapolation assumptions','continuum error not certified','mixed-axis tail not bounded'],reference_acceptance_eligible=False))
        checks=dict(gate_passed=gate['all_criteria_pass'],twenty_fields=len(meshes)==20 and len([r for r in meshes if r['origin']=='fine-stage'])==16,
            fine_equations_and_conservation=all(r.get('eligible',True) for r in meshes),six_series=len(series)==6,
            two_stored_field_bounds=len(pairs)==18 and all(min(p['margins'].values())>=-1e-10 and p['contraction_excess']<=1e-10 for p in pairs),
            resource_budget=time.monotonic()-begin<120 and memory.memory_bytes()<512*1024**2)
    except Exception as exc:stop=repr(exc);raise
    finally:
        result=dict(protocol='rd/aero/avl-reference/records/avl-symmetry-fine-observation-addendum.md',protocol_sha256=digest(FOLDER/'avl-symmetry-fine-observation-addendum.md'),
            source_sha256=digest(Path(__file__)),dependency_sha256={Path(m.__file__).name:digest(Path(m.__file__)) for m in [expansion,observation,equivalence,memory]},
            parent_sha256={name:digest(FOLDER/name) for name in ['avl-symmetry-pi-fine-runs.json','avl-symmetry-pi-equivalence-result.json','avl-spacing-pi-result.json']},
            meshes=meshes,series=series,pairs=pairs,budgets=budgets,checks=checks,methods_criteria_pass=bool(checks) and all(checks.values()),
            stop_reason=stop,reference_convergence_accepted=False,continuum_error_bound_available=False,
            seconds=time.monotonic()-begin,peak_process_working_bytes=memory.memory_bytes())
        write(output,result)
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],budgets=[{k:b[k] for k in ['case','base','conditional_lift_relative','conditional_moment_relative','conditional_distribution_relative','conditional_targets_met']} for b in budgets]),indent=2))

if __name__=='__main__':main()
