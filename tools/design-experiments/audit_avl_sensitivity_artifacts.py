"""Recheck provenance and preserved failed attempts without rerunning AVL."""
from pathlib import Path
import hashlib
import json
import math

ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification';BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(name):return json.loads((FOLDER/name).read_text())
def finite(x):
    if isinstance(x,float):assert math.isfinite(x)
    elif isinstance(x,dict):
        for v in x.values():finite(v)
    elif isinstance(x,list):
        for v in x:finite(v)

def main():
    runs=read('avl-solve-sensitivity-runs.json');result=read('avl-solve-sensitivity-result.json')
    prep=read('avl-sensitivity-preparation-corrected.json')
    assert runs['source_sha256']==digest(Path(__file__).with_name('run_avl_sensitivity.py'))
    assert prep['source_sha256']==digest(Path(__file__).with_name('prepare_avl_sensitivity.py'))
    assert prep['adapter_sha256']==digest(Path(__file__).with_name('reference_avl_sensitivity.f'))
    assert result['source_sha256']==digest(Path(__file__).with_name('check_avl_sensitivity.py'))
    assert result['runs_sha256']==digest(FOLDER/'avl-solve-sensitivity-runs.json')
    assert result['parent_linear_result_sha256']==digest(FOLDER/'avl-linear-load-wrench-result.json')
    assert runs['preparation_sha256']==digest(FOLDER/'avl-sensitivity-preparation-corrected.json')
    assert prep['build_script_sha256']==digest(BASE/'build-sensitivity-corrected.sh')
    assert prep['patch_sha256']==digest(FOLDER/'fixtures/avl-measurement/solve-sensitivity-corrected.patch')
    assert prep['modified_source_sha256']==digest(BASE/'source-core0-sensitivity-corrected/src/aoper.f')
    assert prep['unchanged_kernel_sha256']==digest(BASE/'source-core0-grid/src/aic.f')==digest(BASE/'source-core0-sensitivity-corrected/src/aic.f')
    assert (BASE/'source-core0-sensitivity-corrected/src/uvl_sensitivity.f').read_bytes()==Path(__file__).with_name('reference_avl_sensitivity.f').read_bytes()
    assert runs['executable_sha256']==digest(BASE/'source-core0-sensitivity-corrected/bin/avl')
    for parent in runs['parents'].values():assert digest(ROOT/parent['path'])==parent['sha256']
    for name,value in result['dependency_sha256'].items():assert digest(Path(__file__).with_name(name))==value
    assert len(runs['runs'])==8 and runs['completed'] and runs['stop_reason'] is None
    files=0;solve_status=[]
    adapter=Path(__file__).with_name('reference_avl_sensitivity.f').read_text()
    assert 'IF(UINFO.NE.0) STOP 83' in adapter
    for run in runs['runs']:
        assert run['exit_code']==0 and run['stopped'] is None and run['seconds']<180
        assert run['peak_resident_bytes']<2*1024**3 and all(run['unchanged'].values())
        summary={k:float(v) for k,v in (line.split() for line in (ROOT/run['directory']/'out/us-summary.txt').read_text().splitlines())}
        assert 'UVL SENSITIVITY: read-only diagnostics written' in (ROOT/run['directory']/'console.txt').read_text()
        solve_status.append(dict(directory=run['directory'],dgetrs_calls=47+int(summary['NORM_SOLVES_1']+summary['NORM_SOLVES_INF']),
            info=0,evidence='Inferred from the unconditional nonzero-INFO abort guard, normal US completion, and successful post-US output; not a separately logged per-call INFO array.',refinement_steps=2))
        for name,info in run['files'].items():
            p=ROOT/run['directory']/name;assert digest(p)==info['sha256'] and p.stat().st_size==info['bytes'];files+=1
        for name in run['unchanged']:
            assert digest(ROOT/run['original_directory']/name)==digest(ROOT/run['directory']/'out'/name)==digest(ROOT/run['directory']/'after'/name)
    assert runs['seconds']<900 and runs['total_bytes']<256*1024**2
    assert len(result['checks'])==11 and all(result['checks'].values()) and result['all_criteria_pass']
    assert not result['reference_convergence_accepted'] and len(result['circulation_ledger'])==40 and len(result['grid_vs_solve'])==6
    assert result['seconds']<180 and result['peak_process_working_bytes']<512*1024**2
    old=read('avl-solve-sensitivity-uncompensated-runs.json')
    assert old['source_sha256']==digest(BASE/'run_avl_sensitivity.v1.py') and old['completed']
    initialprep=read('avl-sensitivity-preparation-initial.json')
    assert initialprep['adapter_sha256']==digest(BASE/'reference_avl_sensitivity.initial.f')
    firstprep=read('avl-sensitivity-preparation.json')
    assert firstprep['adapter_sha256']==digest(BASE/'reference_avl_sensitivity.v1.f')
    assert firstprep['source_sha256']==digest(BASE/'prepare_avl_sensitivity.v1.py')
    assert 'Unterminated character constant' in (BASE/'build-sensitivity.initial.log').read_text()
    assert not read('avl-solve-sensitivity-run-initial-status.json')['completed']
    assert read('avl-sensitivity-analysis-initial-status.json')['status']=='postprocessing_serialization_failed'
    assert 'UC=UC+((US-UT)+UX)' in Path(__file__).with_name('reference_avl_sensitivity.f').read_text()
    lapack=json.loads((BASE/'norm-estimator/source-manifest.json').read_text())
    assert prep['lapack_manifest_sha256']==digest(BASE/'norm-estimator/source-manifest.json')
    # The manifest includes each official source's URL, byte count and digest.
    for item in lapack['files']:
        assert digest(BASE/'norm-estimator'/item['name'])==item['sha256']
    finite(result);finite(runs)
    report=dict(source_sha256=digest(Path(__file__)),result_sha256=digest(FOLDER/'avl-solve-sensitivity-result.json'),
        runs_sha256=digest(FOLDER/'avl-solve-sensitivity-runs.json'),raw_files_verified=files,
        systems=8,criteria=11,solve_status=solve_status,previous_attempts_preserved=True,all_criteria_pass=True)
    (FOLDER/'avl-solve-sensitivity-artifact-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
