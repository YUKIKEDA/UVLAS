"""Provenance and numerical-identity audit; no AVL process is launched."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,math

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS;BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(name):return json.loads((FOLDER/name).read_text())
def finite(x):
    if isinstance(x,float):assert math.isfinite(x)
    elif isinstance(x,dict):
        for v in x.values():finite(v)
    elif isinstance(x,list):
        for v in x:finite(v)
def source(data,name):assert data['source_sha256']==digest(rd_paths.sibling(name))
def dependencies(data):
    for name,value in data['dependency_sha256'].items():assert digest(rd_paths.sibling(name))==value
    for name,value in data['parent_sha256'].items():assert digest(FOLDER/name)==value

def main():
    r=read('avl-conservative-observation-result.json');source(r,'check_avl_conservative_observation.py');dependencies(r)
    assert r['protocol_sha256']==digest(FOLDER/'avl-conservative-observation-experiment.md')
    assert r['methods_criteria_pass'] and len(r['checks'])==9 and all(r['checks'].values()) and r['stop_reason'] is None
    assert not r['reference_convergence_accepted'] and not r['continuum_error_bound_available']
    assert len(r['meshes'])==40 and len(r['pairs'])==36 and len(r['series'])==12 and len(r['budgets'])==4
    unit=r['manufactured'];assert len(unit['cases'])==4 and len(unit['checks'])==6 and all(unit['checks'].values())
    assert len(unit['failures'])==7 and all(f['refused'] and f['input_unchanged'] for f in unit['failures'])
    assert max(c['omit_T1_difference'] for c in unit['cases'])>1e-5 and max(c['omit_couple_difference'] for c in unit['cases'])>1e-5
    for mesh in r['meshes']:
        assert mesh['observation_id']=='tau-k-band-moments-v1' and len(mesh['bands'])==20
        assert max(mesh[k] for k in ['load_identity_difference','previous_ledger_difference','independent_band_difference','p1_moment_preservation'])<=1e-10
        for info in mesh['inputs'].values():assert digest(ROOT/info['path'])==info['sha256']
        for b in mesh['bands']:
            assert b['width']>0 and abs(b['tau0']*b['width']-b['T0'])<=1e-10 and abs(b['tau1']*b['width']**2/6-b['T1'])<=1e-10
    for pair in r['pairs']:
        assert min(pair['margins'].values())>=-1e-10 and pair['contraction_excess']<=1e-10
        assert set(pair['bounds'])==set(pair['margins']) and len(pair['bounds'])==8
    for b in r['budgets']:
        assert not b['reference_acceptance_eligible']
        for c in b['channels'].values():
            if c['missing_axes']:assert c['combined_conditional_relative'] is None and c['state']=='unassessed_missing_axis'
    initial=read('avl-conservative-observation-initial-result.json')
    assert initial['source_sha256']==digest(ROOT/'build/analysis/check_avl_conservative_observation_initial.py')
    assert not initial['methods_criteria_pass'] and 'outside_wing' in initial['stop_reason'] and not initial['meshes']
    assert r['seconds']<300 and r['peak_process_working_bytes']<512*1024**2
    finite(r)
    report=dict(source_sha256=digest(Path(__file__)),result_sha256=digest(FOLDER/'avl-conservative-observation-result.json'),
        method_checks=9,manufactured_cases=4,rejected_cases=7,meshes=40,pair_bounds=36,bounds_per_pair=8,
        initial_failure_preserved=True,all_criteria_pass=True)
    (FOLDER/'avl-conservative-observation-artifact-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if not (FOLDER/'avl-span320-result.json').exists():return
    prep=read('avl-span320-preparation.json');source(prep,'prepare_avl_span320.py')
    assert sorted(prep['changed_source_files'])==['ADIMEN.INC','avl_heap_inc.f90'] and prep['capacity']==8000
    for name,value in prep['source_hashes'].items():
        assert digest(BASE/'source-core0-grid/src'/name)==value['original']
        assert digest(BASE/'source-core0-grid8000/src'/name)==value['capacity8000']
    assert prep['parent_executable_sha256']==digest(BASE/'source-core0-grid/bin/avl')
    assert prep['parent_input_sha256']==digest(FOLDER/'avl-band-aligned-inputs.json')
    assert prep['build_script_sha256']==digest(BASE/'build-grid8000.sh')
    assert prep['patch_sha256']==digest(FOLDER/'fixtures/avl-measurement/capacity8000.patch')
    parent=read('avl-band-aligned-inputs.json')
    for deck in prep['decks']:
        assert digest(ROOT/deck['path'])==deck['sha256'] and len(deck['sections'])==21 and abs(deck['area_xy']-8.4)<1e-13
        case=next(c for c in parent['cases'] if c['id']==deck['case']);assert deck['sections']==case['sections']
        old=next(d for d in case['decks'] if (d['m'],d['n'])==(24,160));text=(FOLDER/old['path']).read_bytes()
        assert text.count(b' 8 1.0\r\n')==21
        assert text.replace(b' 8 1.0\r\n',b' 16 1.0\r\n')==(ROOT/deck['path']).read_bytes()
    runs=read('avl-span320-runs.json');source(runs,'run_avl_span320.py')
    assert runs['completed'] and runs['stop_reason'] is None and len(runs['runs'])==6
    assert runs['protocol_sha256']==digest(FOLDER/'avl-span320-reference-experiment.md')
    assert runs['preparation_sha256']==digest(FOLDER/'avl-span320-preparation.json')
    assert runs['parent_runs_sha256']==digest(FOLDER/'avl-band-aligned-grid-runs.json')
    assert runs['executable_sha256']==digest(BASE/'source-core0-grid8000/bin/avl')
    assert runs['executable_sha256']!=prep['parent_executable_sha256']
    files=0
    for run in runs['runs']:
        assert run['exit_code']==0 and run['stopped'] is None and run['seconds']<300 and run['peak_resident_bytes']<4*1024**3
        assert run['available_memory_before']>=6*1024**3 and all(run['before_after_equal'].values())
        if run['kind']=='regression':
            assert all(run['previous_outputs_equal'].values())
            for name in run['previous_outputs_equal']:assert digest(ROOT/run['previous_directory']/name)==digest(ROOT/run['directory']/'out'/name)
        for name,value in run['files'].items():
            p=ROOT/run['directory']/name;assert digest(p)==value['sha256'] and p.stat().st_size==value['bytes'];files+=1
    assert runs['seconds']<900 and runs['total_bytes']<128*1024**2
    s=read('avl-span320-result.json');source(s,'check_avl_span320.py');dependencies(s)
    assert len(s['checks'])==5 and all(s['checks'].values()) and s['methods_criteria_pass']
    assert len(s['new_cases'])==2 and len(s['series'])==2 and len(s['pairs'])==6
    assert not s['reference_convergence_accepted'] and not s['continuum_error_bound_available']
    for b in s['budgets']:assert not b['reference_acceptance_eligible'] and b['unverified']
    assert s['seconds']<120 and s['peak_process_working_bytes']<512*1024**2
    finite(runs);finite(s)
    audit=dict(source_sha256=digest(Path(__file__)),result_sha256=digest(FOLDER/'avl-span320-result.json'),
        runs_sha256=digest(FOLDER/'avl-span320-runs.json'),regression_cases=4,new_cases=2,raw_files=files,
        changed_source_parameters=2,method_checks=5,all_criteria_pass=True)
    (FOLDER/'avl-span320-artifact-audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))

if __name__=='__main__':main()
