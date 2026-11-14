"""Generate a separately versioned staged runner from the validated archive runner."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;HERE=Path(__file__).parent
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    target=HERE/'run_avl_symmetric_solve.py'
    if target.exists():raise FileExistsError(target)
    s=(HERE/'run_avl_tip_redistribution.py').read_text()
    s=s.replace('import datetime,hashlib,json,os,re,shutil,subprocess,time,zipfile','import datetime,hashlib,json,math,os,re,shutil,subprocess,sys,time,zipfile')
    s=s.replace('import avl_tip_watchdog as guard','import avl_symmetric_solve_watchdog as guard')
    start=s.index('    output=');end=s.index('    root=ROOT/',start)
    header='''    stage='main' if '--main' in sys.argv[1:] else 'small'
    output=F/f'avl-symmetric-solve-{stage}-runs.json'
    if output.exists():raise FileExistsError(output)
    prep_path=F/'avl-symmetric-solve-preparation.json';prep=json.loads(prep_path.read_text())
    assert digest(ROOT/prep['protocol'])==prep['protocol_sha256']
    gates=['avl-tip-redistribution-artifact-audit.json','avl-symmetric-solve-watchdog-result.json']
    for n in gates:assert json.loads((F/n).read_text())['all_criteria_pass']
    watchdog=json.loads((F/gates[1]).read_text());assert watchdog['watchdog_sha256']==digest(guard.__file__)
    assert watchdog['source_sha256']==digest(rd_paths.sibling('check_avl_symmetric_solve_watchdog.py')) and watchdog['termination_sha256']==digest(termination.__file__)
    assert prep['weights_sha256']==digest(F/'avl-symmetric-solve-weight-result.json')
    for n,h in prep['source_files'].items():assert digest(ROOT/prep['target']/'src'/n)==h['probe']
    assert digest(ROOT/prep['target']/'src'/prep['new_adapter'])==prep['adapter_sha256']
    exe=ROOT/prep['target']/'bin/avl'
    if stage=='main':
        gate=F/'avl-symmetric-solve-small-check.json';g=json.loads(gate.read_text());assert g['all_criteria_pass']
        assert g['preparation_sha256']==digest(prep_path) and g['small_runs_sha256']==digest(F/'avl-symmetric-solve-small-runs.json')
        gates.append(gate.name)
    selections=[r for r in prep['cases'] if (r['allocation']=='small')==(stage=='small')]
    assert len(selections)==(2 if stage=='small' else 6)
'''
    s=s[:start]+header+s[end:]
    s=s.replace("'avl-tip-redistribution-'+datetime", "'avl-symmetric-solve-'+stage+'-'+datetime")
    s=s.replace("execution_protocol='rd/aero/avl-reference/records/avl-tip-redistribution-execution-addendum.md'", "execution_protocol=prep['protocol']")
    s=s.replace("digest(F/'avl-tip-redistribution-execution-addendum.md')", "prep['protocol_sha256']")
    s=s.replace("capacity_preparation_sha256=digest(F/'avl-span1280-capacity-ks-preparation.json')", "stage=stage,build_script_sha256=prep['build_script_sha256'],build_log_sha256=digest(ROOT/'build/reference-avl/build-symmetric-solve-probe.log'),build_resource_sha256=digest(ROOT/'build/reference-avl/build-symmetric-solve-probe-resource.txt')")
    old="        for deck in sorted(prep['decks'],key=lambda d:(d['case'],0 if d['allocation']=='tip-rich' else 1)):"
    assert old in s;s=s.replace(old,"        for deck in sorted(selections,key=lambda d:(d['case'],{'small':0,'uniform':0,'tip-rich':1,'tip-poor':2}[d['allocation']])):\n            names=NAMES+(['uh-matrix.bin'] if deck['native_unknowns']<=256 else [])")
    s=s.replace("time.monotonic()-start>1200", "time.monotonic()-start>guard.LIMITS[stage][1]")
    s=s.replace("source=ROOT/deck['path'];assert digest(source)==deck['sha256']", "source=ROOT/deck['input'];assert digest(source)==deck['input_sha256']")
    s=s.replace("['out','after']", "['out','after','uq']")
    old="            for name in ['out/','after/']:commands += [f'ft {name}total.txt',f'fs {name}strips.txt',f'fe {name}elements.txt',f'uh {name}']"
    assert old in s;s=s.replace(old,"            for name in ['out/','after/']:\n                if name=='after/':commands+=['uq uq/']\n                commands += [f'ft {name}total.txt',f'fs {name}strips.txt',f'fe {name}elements.txt',f'uh {name}']")
    s=s.replace('guard.stop_reason(', 'guard.stop_reason(stage,')
    s=s.replace('in NAMES','in names')
    marker="            record['archive_members']=";pos=s.index(marker)
    check='''            record['parent_outputs_equal']={n:digest(directory/'out'/n)==digest(ROOT/deck['raw_directory']/n) for n in names}
            assert all(record['parent_outputs_equal'].values())
            meta={v.split()[0]:int(v.split()[1]) for v in (directory/'uq/uq-meta.txt').read_text().splitlines()}
            assert meta['N']==deck['native_unknowns'] and meta['REAL_BITS']==64 and meta['QUAD_BITS']==128 and meta['QUAD_DIGITS']>=30 and meta['STATE_UNCHANGED']==1
            tests=[list(map(int,line.split())) for line in (directory/'uq/uq-guard-tests.txt').read_text().splitlines()[1:]]
            assert len(tests)==8 and all((status!=0)==bool(expected) for _,status,expected in tests)
            vectors=(directory/'uq/uq-vectors.txt').read_text().splitlines()[1:]
            assert len(vectors)==deck['native_unknowns']
            assert all(len(v.split())==10 and all(math.isfinite(float(x)) for x in v.split()) for v in vectors)
            record['probe_meta']=meta
'''
    s=s[:pos]+check+s[pos:]
    s=s.replace('len(runs)==4','len(runs)==len(selections)').replace("len(r.get('before_after_equal',{}))==8", "len(r.get('before_after_equal',{}))==(9 if r['native_unknowns']<=256 else 8)")
    target.write_text(s)
    (F/'avl-symmetric-solve-runner-preparation.json').write_text(json.dumps(dict(source_sha256=digest(__file__),template_sha256=digest(HERE/'run_avl_tip_redistribution.py'),runner_sha256=digest(target)),indent=2)+'\n')
    print('Generated separately staged UQ runner')
if __name__=='__main__':main()
