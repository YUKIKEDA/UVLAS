"""Audit the retained failure, one local-capacity change, eight loads and eight regressions."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,re,hashlib,shlex
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;GIB=1024**3
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(n):return json.loads((F/n).read_text())
def files(folder,manifest):
    for n,h in manifest.items():assert digest(folder/n)==h['sha256'] and (folder/n).stat().st_size==h['bytes']
def main():
    p=read('avl-span1280-capacity-ks-preparation.json');r=read('avl-span1280-capacity-ks-runs.json');load=read('avl-span1280-load-result.json');failed=read('avl-span1280-fine-runs.json')
    assert not failed['completed'] and len(failed['runs'])==1 and failed['stop_reason']
    assert p['failed_runs_sha256']==digest(F/'avl-span1280-fine-runs.json')
    assert failed['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/run_avl_span1280_fine.py')
    for run in failed['runs']:
        folder=ROOT/run['directory'];files(folder,run['files']);assert 'Increase KSMAX to' in (folder/'console.txt').read_text()
        assert not (folder/'out/uv-rows.txt').exists()
    assert p['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/prepare_avl_span1280_ksmax.py')
    assert r['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/run_avl_span1280_capacity_ks.py')
    assert load['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/check_avl_span1280_load.py')
    assert p['protocol_sha256']==r['protocol_sha256']==load['protocol_sha256']==digest(ROOT/p['protocol'])
    assert load['display_protocol_sha256']==digest(F/'avl-span1280-load-diagnostics-addendum.md')
    parent=ROOT/p['parent'];target=ROOT/p['target'];changes=[]
    for n,h in p['source_files'].items():
        assert digest(parent/'src'/n)==h['parent'] and digest(target/'src'/n)==h['capacity']
        if h['parent']!=h['capacity']:changes.append(n)
    assert changes==['amake.f'] and p['changed_constants']==1
    assert (parent/'src/amake.f').read_bytes().replace(b'KSMAX=500',b'KSMAX=1000')==(target/'src/amake.f').read_bytes()
    assert len(p['reused_object_sha256'])==40
    for n,h in p['reused_object_sha256'].items():assert digest(parent/'bin'/n)==h==digest(target/'bin'/n)
    assert digest(parent/'bin/avl')==p['parent_executable_sha256']
    assert digest(target/'bin/avl')==r['executable_sha256']==load['executable_sha256']
    base=ROOT/'build/reference-avl';script=base/'build-symmetry16000-ks1000.sh'
    assert digest(script)==p['build_script_sha256']==r['build_script_sha256']
    assert digest(F/'fixtures/avl-measurement/capacity-ks1000.patch')==p['patch_sha256']
    commands=[shlex.split(s) for s in script.read_text().splitlines() if ' -c ' in s]
    assert len(commands)==1 and commands[0][-1]=='amake.o' and '-fdefault-real-8' in commands[0]
    linux='<repo>'
    for n,h in p['link_inputs'].items():assert digest(ROOT/n[len(linux)+1:])==h
    assert digest(base/'build-symmetry16000-ks1000.log')==r['build_log_sha256']
    # The generated runner retained the parent resource filename; preserve that record
    # and independently bind/validate the actual correction-build resource below.
    assert digest(base/'build-symmetry16000-resource.txt')==r['build_resource_sha256']
    build=(base/'build-symmetry16000-ks1000-resource.txt').read_text();assert re.search(r'Exit status:\s*0',build)
    assert int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',build)[1])*1024<4*GIB
    assert r['preparation_sha256']==digest(F/'avl-span1280-capacity-ks-preparation.json')
    for n,h in r['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/n)==h
    assert r['completed'] and len(r['runs'])==8 and r['seconds']<900 and r['total_bytes']<256*1024**2
    expected={(c,m,n) for c in ['FW-V01','FW-V02'] for m,n in [(6,40),(24,160),(24,320),(24,640)]}
    assert {(v['case'],v['m'],v['physical_n']) for v in r['runs']}==expected
    count=0;size=0
    for v in r['runs']:
        folder=ROOT/v['directory'];files(folder,v['files']);size+=sum(h['bytes'] for h in v['files'].values())
        assert v['exit_code']==0 and v['stopped'] is None and v['seconds']<300 and max(v['peak_resident_bytes'],v['maximum_group_live_bytes'])<5*GIB
        assert v['memory_before']['effective_available']>=8*GIB
        assert digest(folder/'wing.avl')==v['input_sha256']==digest(ROOT/v['input'])
        assert all(v['before_after_equal'].values()) and all(v['parent_outputs_equal'].values())
        assert len(v['parent_file_sha256'])==(9 if v['native_unknowns']<=256 else 8)
        for n,h in v['parent_file_sha256'].items():assert digest(ROOT/v['parent_directory']/n)==h==digest(folder/'out'/n)==digest(folder/'after'/n);count+=1
    assert size==r['total_bytes'] and count==66
    assert load['all_criteria_pass'] and load['no_solve'] and len(load['cases'])==8 and load['seconds']<60 and load['total_bytes']<8*1024**2
    assert sum(v['i4_display_overflow'] for v in load['cases'])==2
    for v in load['cases']:
        files(ROOT/v['directory'],v['files']);assert v['strips']==640 and v['peak_bytes']<256*1024**2
    result=dict(source_sha256=digest(__file__),preparation_sha256=digest(F/'avl-span1280-capacity-ks-preparation.json'),
        runs_sha256=digest(F/'avl-span1280-capacity-ks-runs.json'),load_sha256=digest(F/'avl-span1280-load-result.json'),
        failed_runs_sha256=digest(F/'avl-span1280-fine-runs.json'),recompiled_objects=1,reused_objects=40,changed_constants=1,
        load_cases=8,regression_cases=8,compared_outputs=count,regression_seconds=r['seconds'],
        peak_resident_bytes=max(v['peak_resident_bytes'] for v in r['runs']),
        runner_recorded_parent_build_resource=True,corrected_build_resource_sha256=digest(base/'build-symmetry16000-ks1000-resource.txt'),
        all_criteria_pass=True,reference_acceptance=False)
    (F/'avl-span1280-capacity-ks-artifact-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
