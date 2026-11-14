"""Hash/resource/field audit of corrected N1280 solves and conditional observations."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,math
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;GIB=1024**3
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(n):return json.loads((F/n).read_text())
def main():
    prep=read('avl-span1280-fine-ks-preparation.json');raw=read('avl-span1280-fine-ks-runs.json');obs=read('avl-span1280-observation-result.json')
    assert prep['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/finalize_avl_span1280_ksmax.py')
    assert raw['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/run_avl_span1280_fine_ks.py')
    assert obs['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/check_avl_span1280_observations.py')
    assert prep['protocol_sha256']==raw['protocol_sha256']==obs['protocol_sha256']==digest(ROOT/prep['protocol'])
    assert prep['correction_protocol_sha256']==digest(F/'avl-span1280-ksmax-addendum.md')
    assert prep['failed_runs_sha256']==digest(F/'avl-span1280-fine-runs.json')
    assert prep['capacity_audit_sha256']==digest(F/'avl-span1280-capacity-ks-artifact-audit.json')
    assert prep['parent_runs_sha256']==digest(F/'avl-span1280-capacity-ks-runs.json')
    assert raw['preparation_sha256']==digest(F/'avl-span1280-fine-ks-preparation.json')
    for n,h in raw['gate_sha256'].items():assert digest(F/n)==h and read(n)['all_criteria_pass']
    for d in [raw,obs]:
        for n,h in d['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/n)==h
    for n,h in obs['parent_sha256'].items():assert digest(F/n)==h
    exe=ROOT/'build/reference-avl/source-symmetry16000-ks1000-spacing-pi64/bin/avl'
    assert digest(exe)==prep['executable_sha256']==raw['executable_sha256']
    for deck in prep['decks']:
        donor=ROOT/deck['donor'];p=ROOT/deck['path'];assert digest(donor)==deck['donor_sha256'] and digest(p)==deck['sha256']
        expected=donor.read_bytes().replace(b' 32 1.0\r\n',b' 64 1.0\r\n').replace(b'24 1.0\r\n',f'{deck["m"]} 1.0\r\n'.encode())
        assert p.read_bytes()==expected and deck['native_strips']==640 and deck['native_unknowns']==640*deck['m']
    assert raw['completed'] and raw['stop_reason'] is None and len(raw['runs'])==8 and raw['seconds']<2400
    assert {(r['case'],r['m']) for r in raw['runs']}=={(c,m) for c in ['FW-V01','FW-V02'] for m in [3,6,12,24]}
    size=0;files=0
    for r in raw['runs']:
        assert r['physical_n']==1280 and r['native_unknowns']==r['m']*640
        assert r['exit_code']==0 and r['stopped'] is None and r['seconds']<1200
        assert max(r['peak_resident_bytes'],r['maximum_group_live_bytes'])<16*GIB
        assert r['memory_before']['cgroup_known'] and r['memory_before']['effective_available']>=20*GIB
        assert r['minimum_effective_available']>=2*GIB and r['memory_after']['cgroup_known']
        assert r['maximum_swap_used']-min(r['memory_before']['swap_used'],raw['stage_memory_before']['swap_used'])<=256*1024**2
        folder=ROOT/r['directory'];assert digest(folder/'wing.avl')==r['sha256']
        assert len(r['before_after_equal'])==8 and all(r['before_after_equal'].values())
        for n in r['before_after_equal']:assert digest(folder/'out'/n)==digest(folder/'after'/n)
        for n,h in r['files'].items():
            p=folder/n;assert digest(p)==h['sha256'] and p.stat().st_size==h['bytes'];size+=h['bytes'];files+=1
    assert size==raw['total_bytes'] and size<512*1024**2
    assert obs['methods_criteria_pass'] and all(obs['checks'].values()) and obs['stop_reason'] is None
    assert len(obs['meshes'])==14 and len(obs['series'])==4 and len(obs['pairs'])==12 and len(obs['budgets'])==2
    assert obs['seconds']<120 and obs['peak_process_working_bytes']<GIB
    for r in obs['meshes']:
        assert r['eligible'] and r['residual_relative']<=1e-11 and r['native_unknowns']*2==r['physical_panels']
        folder=ROOT/r['derived_directory']
        for n,h in r['files'].items():assert digest(folder/n)==h['sha256'] and (folder/n).stat().st_size==h['bytes']
    for s in obs['series']:
        assert len(s['channels'])==5
        for c in s['channels'].values():
            value=c['conditional_extrapolation']
            if c['status']=='asymptotic_candidate':
                assert all(c['differences'][j+1]<c['differences'][j] for j in range(2))
                assert .5<=c['orders'][-1]<=6 and abs(c['orders'][-1]-c['orders'][-2])<=.5 and min(c['difference_cosines'])>=.9
                assert value is not None and math.isfinite(value) and value>=0
            elif c['status']!='roundoff_candidate':assert value is None
    for b in obs['budgets']:
        assert b['base']=='M24N1280' and b['same_base_without_transfer'] and not b['reference_acceptance_eligible']
        for c in b['channels'].values():
            if c['chord_conditional'] is None or c['span_conditional'] is None:assert c['combined_conditional'] is None
            else:assert math.isclose(c['combined_conditional'],c['chord_conditional']+c['span_conditional'],rel_tol=1e-14,abs_tol=1e-18)
    assert not raw['reference_acceptance'] and not obs['reference_convergence_accepted'] and not obs['continuum_error_bound_available']
    preview=read('avl-span1280-first-wing-result.json');checkpoint=read('avl-span1280-first-wing-checkpoint.json')
    assert not checkpoint['completed'] and sum(r['case']=='FW-V01' for r in checkpoint['runs'])==4
    assert preview['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/check_avl_span1280_first_wing.py')
    assert preview['methods_criteria_pass'] and not preview['reference_convergence_accepted']
    for n,h in preview['parent_sha256'].items():assert digest(F/n)==h
    assert preview['series']==[s for s in obs['series'] if s['case']=='FW-V01']
    assert preview['budgets']==[b for b in obs['budgets'] if b['case']=='FW-V01']
    result=dict(source_sha256=digest(__file__),runs_sha256=digest(F/'avl-span1280-fine-ks-runs.json'),observations_sha256=digest(F/'avl-span1280-observation-result.json'),
        completed_fine_cases=8,physical_fields=14,series=4,stored_field_pairs=12,raw_files=files,
        seconds=raw['seconds'],peak_resident_bytes=max(r['peak_resident_bytes'] for r in raw['runs']),raw_bytes=size,
        initial_failure_preserved=True,separate_first_wing_checkpoint_preserved=True,reference_acceptance=False,all_criteria_pass=True)
    (F/'avl-span1280-fine-artifact-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
