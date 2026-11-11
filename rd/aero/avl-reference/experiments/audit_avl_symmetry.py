"""Audit original failed gate, corrected spacing gate and native/derived fine data."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from decimal import Decimal
import hashlib,json,math

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS;BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(name):return json.loads((FOLDER/name).read_text())
def source(d,name):assert d['source_sha256']==digest(rd_paths.sibling(name))
def dependencies(d):
    for name,value in d['dependency_sha256'].items():assert digest(rd_paths.sibling(name))==value
def finite(x):
    if isinstance(x,float):assert math.isfinite(x)
    elif isinstance(x,dict):
        for v in x.values():finite(v)
    elif isinstance(x,list):
        for v in x:finite(v)
def files(folder,records):
    for name,v in records.items():
        p=folder/name;assert digest(p)==v['sha256'] and p.stat().st_size==v['bytes']
def runs(name,script,count,limit):
    d=read(name);source(d,script);assert d['completed'] and len(d['runs'])==count and d['stop_reason'] is None
    assert d['protocol_sha256']==digest(ROOT/d['protocol']) and d['seconds']<limit
    assert d['resource_runner_sha256']==digest(rd_paths.sibling('run_avl_span320.py'))
    for r in d['runs']:
        assert r['exit_code']==0 and r['stopped'] is None and r['seconds']<300 and r['peak_resident_bytes']<4*1024**3
        assert r['available_memory_before']>=6*1024**3 and all(r['before_after_equal'].values())
        if r['kind']=='half_measured':assert all(r['parent_half_outputs_equal'].values())
        folder=ROOT/r['directory'];files(folder,r['files'])
        matrix=folder/'out/uh-matrix.bin'
        if matrix.exists():assert digest(matrix)==digest(folder/'after/uh-matrix.bin') and matrix.stat().st_size==r['native_unknowns']**2*8
    finite(d);return d

def main():
    original=read('avl-symmetry-preparation.json');source(original,'prepare_avl_symmetry.py')
    for name,sha in original['adapter_sha256'].items():assert digest(rd_paths.sibling(name))==sha
    for name,h in original['source_hashes'].items():
        assert digest(BASE/'source-core0-grid8000/src'/name)==h['parent']
        assert digest(BASE/'source-symmetry8000/src'/name)==h['measured']
    assert original['build_script_sha256']==digest(BASE/'build-symmetry.sh')
    assert original['patch_sha256']==digest(FOLDER/'fixtures/avl-measurement/symmetry-export.patch')
    failed_runs=runs('avl-symmetry-regression-runs.json','run_avl_symmetry.py',14,600)
    assert failed_runs['preparation_sha256']==digest(FOLDER/'avl-symmetry-preparation.json')
    failed=read('avl-symmetry-equivalence-initial-result.json');source(failed,'check_avl_symmetry_equivalence.py');dependencies(failed)
    assert not failed['all_criteria_pass'] and failed['completed'] and not failed['checks']['six_full_equivalences']
    assert len(failed['records'])==6 and all(r['geometry_difference']>1e-11 for r in failed['records'])
    for r in failed['records']:files(ROOT/r['derived_directory'],r['files'])
    assert digest(FOLDER/'avl-symmetry-equivalence-result.json')==digest(FOLDER/'avl-symmetry-equivalence-initial-result.json')
    spacing=read('avl-spacing-pi-result.json');source(spacing,'check_avl_spacing_pi.py')
    assert len(spacing['checks'])==4 and all(spacing['checks'].values()) and spacing['all_criteria_pass']
    assert Decimal(spacing['pi_difference'])<=Decimal('1e-75')
    for name,sha in spacing['files'].items():assert digest(BASE/name)==sha
    prep=read('avl-symmetry-pi-preparation.json');source(prep,'prepare_avl_symmetry_pi.py')
    assert prep['parent_preparation_sha256']==digest(FOLDER/'avl-symmetry-preparation.json')
    assert prep['build_script_sha256']==digest(BASE/'build-symmetry-pi64.sh')
    assert prep['driver_sha256']==digest(rd_paths.sibling('reference_avl_spacing_driver.f'))
    for v in prep['variants']:
        changed=[]
        for name,h in v['sources'].items():
            assert digest(BASE/v['parent']/'src'/name)==h['parent']
            assert digest(BASE/v['target']/'src'/name)==h['pi64']
            if h['parent']!=h['pi64']:changed.append(name)
        assert changed==['sgutil.f'] and digest(ROOT/v['patch_path'])==v['patch_sha256']
    corrected=runs('avl-symmetry-pi-regression-runs.json','run_avl_symmetry_pi.py',18,900)
    assert corrected['preparation_sha256']==digest(FOLDER/'avl-symmetry-pi-preparation.json')
    assert corrected['spacing_result_sha256']==digest(FOLDER/'avl-spacing-pi-result.json')
    for r in corrected['runs']:
        expected='source-core0-grid8000-spacing-pi64' if r['kind']=='half_baseline' else 'source-symmetry8000-spacing-pi64'
        assert r['executable_sha256']==digest(BASE/expected/'bin/avl')
    gate=read('avl-symmetry-pi-equivalence-result.json');source(gate,'check_avl_symmetry_pi_equivalence.py');dependencies(gate)
    assert gate['regression_runs_sha256']==digest(FOLDER/'avl-symmetry-pi-regression-runs.json')
    assert gate['preparation_sha256']==digest(FOLDER/'avl-symmetry-pi-preparation.json')
    assert gate['spacing_result_sha256']==digest(FOLDER/'avl-spacing-pi-result.json')
    assert gate['completed'] and gate['all_criteria_pass'] and len(gate['checks'])==8 and all(gate['checks'].values())
    assert len(gate['records'])==6 and len(gate['matrix'])==2 and all(m['eligible'] for m in gate['matrix'])
    for r in gate['records']:
        assert r['eligible'] and r['legacy_sensitivity']['eligible'] and r['guard_cases_pass']
        assert max(r['geometry_difference'],r['point_difference'])<=1e-11 and r['gamma_difference']<=1e-11
        files(ROOT/r['derived_directory'],r['files'])
    fineprep=read('avl-symmetry-pi-fine-preparation.json');source(fineprep,'prepare_avl_symmetry_pi_fine.py')
    assert fineprep['parent_preparation_sha256']==digest(FOLDER/'avl-symmetry-pi-preparation.json')
    for d in fineprep['decks']:assert digest(ROOT/d['path'])==d['sha256']
    fine=runs('avl-symmetry-pi-fine-runs.json','run_avl_symmetry_pi_fine.py',16,600)
    assert fine['preparation_sha256']==digest(FOLDER/'avl-symmetry-pi-fine-preparation.json')
    assert fine['gate_sha256']==digest(FOLDER/'avl-symmetry-pi-equivalence-result.json')
    assert fine['spacing_result_sha256']==digest(FOLDER/'avl-spacing-pi-result.json')
    result=read('avl-symmetry-fine-observation-result.json');source(result,'check_avl_symmetry_fine_observations.py');dependencies(result)
    for name,sha in result['parent_sha256'].items():assert digest(FOLDER/name)==sha
    assert result['protocol_sha256']==digest(FOLDER/'avl-symmetry-fine-observation-addendum.md')
    assert result['methods_criteria_pass'] and len(result['checks'])==6 and all(result['checks'].values()) and result['stop_reason'] is None
    assert len(result['meshes'])==20 and len(result['series'])==6 and len(result['pairs'])==18
    for m in result['meshes']:
        for info in m['inputs'].values():assert digest(ROOT/info['path'])==info['sha256']
        if m['origin']=='fine-stage':
            assert m['eligible'] and m['native_unknowns']*2==m['physical_panels'] and m['residual_relative']<=1e-11
            files(ROOT/m['derived_directory'],m['files'])
    assert not result['reference_convergence_accepted'] and not result['continuum_error_bound_available']
    assert len(result['budgets'])==4 and all(not b['reference_acceptance_eligible'] and b['unverified'] for b in result['budgets'])
    assert result['seconds']<120 and result['peak_process_working_bytes']<512*1024**2
    for d in [gate,result,spacing]:finite(d)
    report=dict(source_sha256=digest(Path(__file__)),result_sha256=digest(FOLDER/'avl-symmetry-fine-observation-result.json'),
        initial_failed_comparisons=6,spacing_checks=4,corrected_gate_checks=8,fine_method_checks=6,
        original_processes=14,corrected_regression_processes=18,fine_processes=16,even_static_scope_only=True,
        reference_acceptance=False,all_criteria_pass=True)
    (FOLDER/'avl-symmetry-artifact-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
