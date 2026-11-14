"""Independent provenance/equality audit for capacity-only regressions; no solves."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,math,re,shlex

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS;GIB=1024**3
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(name):return json.loads((FOLDER/name).read_text())

def main():
    prep=read('avl-span1280-capacity-preparation.json');runs=read('avl-span1280-capacity-runs.json');resource=read('avl-span1280-resource-audit.json')
    assert prep['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/prepare_avl_span1280_capacity.py')
    assert runs['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/run_avl_span1280_capacity.py')
    assert resource['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/avl_reference_resources.py')
    guards=read('avl-span1280-resource-guard-result.json')
    assert guards['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/check_avl_reference_resource_guards.py')
    assert guards['resource_source_sha256']==resource['source_sha256'] and guards['all_criteria_pass'] and guards['criteria_passed']==6
    assert all(r['passed'] for r in guards['cases']) and sum(r['fine_start_eligible'] for r in guards['cases'])==1
    for n,h in runs['dependency_sha256'].items():assert h==digest(rd_paths.EXPERIMENTS/n)
    assert runs['preparation_sha256']==digest(FOLDER/'avl-span1280-capacity-preparation.json')
    assert runs['resource_audit_sha256']==digest(FOLDER/'avl-span1280-resource-audit.json')
    assert prep['protocol_sha256']==runs['protocol_sha256']==resource['protocol_sha256']==digest(ROOT/prep['protocol'])
    assert prep['capacity']==16000 and prep['strips']==1000 and prep['changed_constants']==3
    assert prep['changed_source_files']==['ADIMEN.INC','avl_heap_inc.f90']
    parent=ROOT/prep['parent'];target=ROOT/prep['target'];changed=[]
    for n,h in prep['source_files'].items():
        a=parent/'src'/n;b=target/'src'/n;assert digest(a)==h['parent'] and digest(b)==h['capacity']
        if h['parent']!=h['capacity']:changed.append(n)
    assert sorted(changed)==prep['changed_source_files']
    for name,replace in {'ADIMEN.INC':[(b'NVMAX=8000',b'NVMAX=16000'),(b'NSMAX=500',b'NSMAX=1000')],
        'avl_heap_inc.f90':[(b'NVX=8000',b'NVX=16000')]}.items():
        expected=(parent/'src'/name).read_bytes()
        for old,new in replace:assert expected.count(old)==1;expected=expected.replace(old,new)
        assert expected==(target/'src'/name).read_bytes()
    assert digest(parent/'bin/avl')==prep['parent_executable_sha256']
    assert digest(target/'bin/avl')==runs['executable_sha256']
    base=ROOT/'build/reference-avl';script=base/'build-symmetry16000.sh'
    assert digest(script)==prep['build_script_sha256']==runs['build_script_sha256']
    assert digest(FOLDER/'fixtures/avl-measurement/capacity16000.patch')==prep['patch_sha256']
    linux='<repo>'
    def local(p):assert p.startswith(linux+'/');return ROOT/p[len(linux)+1:]
    commands=[shlex.split(line) for line in script.read_text().splitlines() if ' -c ' in line]
    assert len(commands)==len(prep['objects'])==41
    assert [c[c.index('-o')+1] for c in commands]==prep['objects']
    objects={}
    for c in commands:
        obj=c[c.index('-o')+1];source=next(a for a in c if a.endswith(('.f','.f90')))
        assert local(source).parent==target/'src' and '-O' in c and '-fallow-argument-mismatch' in c
        assert ('-fdefault-real-8' in c)==(obj!='matrix-lapacksubs.o')
        assert digest(local(c[0]))==prep['compiler_sha256'];objects[obj]=digest(target/'bin'/obj)
    for n,h in prep['link_inputs'].items():assert digest(local(n))==h
    assert digest(base/'build-symmetry16000.log')==runs['build_log_sha256']
    timing=base/'build-symmetry16000-resource.txt';assert digest(timing)==runs['build_resource_sha256']
    text=timing.read_text();assert re.search(r'Exit status:\s*0',text)
    build_peak=int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',text)[1])*1024
    clock=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)',text)[1]
    parts=[float(v) for v in clock.split(':')];build_seconds=sum(v*60**j for j,v in enumerate(reversed(parts)))
    assert build_peak<4*GIB and build_seconds<600
    assert digest(FOLDER/'avl-symmetry-pi-equivalence-result.json')==prep['gate_sha256'] and read('avl-symmetry-pi-equivalence-result.json')['all_criteria_pass']
    assert digest(FOLDER/'avl-symmetry-pi-fine-runs.json')==prep['parent_runs_sha256']==resource['parent_runs_sha256']
    for n,h in resource['source_files'].items():assert digest(parent/'src'/n)==h
    assert resource['allocated_dense_bytes']==16000**2*8*7 and resource['active_dense_bytes']==15360**2*8*7
    assert math.isclose(resource['cubic_time_scenario'],resource['parent_max_seconds']*8,rel_tol=1e-15)
    assert runs['completed'] and runs['stop_reason'] is None and len(runs['runs'])==8 and runs['seconds']<900
    total=0;compared=0
    assert {(r['case'],r['m'],r['physical_n']) for r in runs['runs']}=={(c,m,n) for c in ['FW-V01','FW-V02'] for m,n in [(6,40),(24,160),(24,320),(24,640)]}
    for r in runs['runs']:
        assert r['physical_n']<=640 and r['native_unknowns']==r['m']*r['physical_n']//2
        assert r['exit_code']==0 and r['stopped'] is None and r['seconds']<300
        assert max(r['peak_resident_bytes'],r['maximum_group_live_bytes'])<5*GIB and r['memory_before']['effective_available']>=8*GIB
        assert r['major_faults'] is not None and r['minor_faults'] is not None
        folder=ROOT/r['directory'];old=ROOT/r['parent_directory']
        assert digest(ROOT/r['input'])==r['input_sha256']==digest(folder/'wing.avl')
        assert set(r['before_after_equal'])==set(r['parent_outputs_equal'])==set(r['parent_file_sha256'])
        assert len(r['parent_outputs_equal'])==(9 if r['native_unknowns']<=256 else 8)
        assert all(r['before_after_equal'].values()) and all(r['parent_outputs_equal'].values())
        for n,h in r['parent_file_sha256'].items():
            assert digest(old/n)==h==digest(folder/'out'/n)==digest(folder/'after'/n);compared+=1
        for n,meta in r['files'].items():
            p=folder/n;assert digest(p)==meta['sha256'] and p.stat().st_size==meta['bytes'];total+=meta['bytes']
    assert total==runs['total_bytes'] and total<256*1024**2
    assert runs['no_n1280_run_performed'] and resource['no_fine_run_performed'] and not runs['reference_acceptance']
    result=dict(source_sha256=digest(__file__),preparation_sha256=digest(FOLDER/'avl-span1280-capacity-preparation.json'),
        runs_sha256=digest(FOLDER/'avl-span1280-capacity-runs.json'),resource_sha256=digest(FOLDER/'avl-span1280-resource-audit.json'),
        changed_source_files=2,changed_constants=3,recompiled_objects=41,object_sha256=objects,resource_guard_cases=6,
        regression_cases=8,compared_outputs=compared,build_seconds=build_seconds,build_peak_bytes=build_peak,
        regression_seconds=runs['seconds'],regression_peak_bytes=max(r['peak_resident_bytes'] for r in runs['runs']),
        regression_bytes=total,no_n1280_run_performed=True,reference_acceptance=False,all_criteria_pass=True)
    (FOLDER/'avl-span1280-capacity-artifact-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='object_sha256'},indent=2))

if __name__=='__main__':main()
