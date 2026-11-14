"""Four fixed-count allocation solves with verified lossless duplicate archiving."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,subprocess,time,zipfile
import avl_reference_resources as policy
import avl_factorial_watchdog as guard
import avl_span1280_watchdog as termination
import run_avl_span320 as monitor
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;digest=policy.digest
NAMES=['total.txt','strips.txt','elements.txt','uv-rows.txt','uv-forces.txt','uv-strips.txt','uh-meta.txt','uh-guard-tests.txt']
def main():
    output=F/'avl-allocation-factorial-runs.json'
    if output.exists():raise FileExistsError(output)
    prep_path=F/'avl-allocation-factorial-preparation.json';prep=json.loads(prep_path.read_text())
    assert digest(ROOT/prep['protocol'])==prep['protocol_sha256']
    gates=['avl-allocation-factorial-preparation-audit.json','avl-symmetric-solve-artifact-audit.json','avl-factorial-watchdog-result.json']
    for n in gates:assert json.loads((F/n).read_text())['all_criteria_pass'],n
    assert json.loads((F/gates[0]).read_text())['preparation_sha256']==digest(prep_path)
    assert prep['solve_audit_sha256']==digest(F/gates[1])
    for name,key in [('avl-span1280-fine-ks-runs.json','parent_sha256'),('avl-tip-redistribution-runs.json','tip_parent_sha256')]:
        assert digest(F/name)==prep[key]
        parent=json.loads((F/name).read_text());assert parent['completed']
        for r in parent['runs']:
            for n,h in r.get('files',{}).items():assert digest(ROOT/r['directory']/n)==h['sha256']
    watchdog=json.loads((F/gates[2]).read_text());assert watchdog['watchdog_sha256']==digest(guard.__file__)
    assert watchdog['source_sha256']==digest(rd_paths.sibling('check_avl_factorial_watchdog.py')) and watchdog['termination_source_sha256']==digest(termination.__file__)
    capacity=json.loads((F/'avl-span1280-capacity-ks-preparation.json').read_text())
    for n,h in capacity['source_files'].items():assert digest(ROOT/capacity['target']/'src'/n)==h['capacity']
    exe=ROOT/capacity['target']/'bin/avl';assert digest(exe)==prep['executable_sha256']
    root=ROOT/'build/analysis'/('avl-allocation-factorial-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S'));root.mkdir()
    environment=os.environ.copy();environment['LD_LIBRARY_PATH']=str(ROOT/'build/reference-avl/toolchain/lib')
    start=time.monotonic();runs=[];stop=None;stage_memory=policy.memory_snapshot()
    def total_bytes():return sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    def save(completed=False):
        result=dict(protocol=prep['protocol'],protocol_sha256=prep['protocol_sha256'],execution_protocol=prep['protocol'],
            execution_protocol_sha256=prep['protocol_sha256'],source_sha256=digest(__file__),
            dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [policy,guard,termination,monitor]},
            preparation_sha256=digest(prep_path),gate_sha256={n:digest(F/n) for n in gates},executable_sha256=digest(exe),
            capacity_preparation_sha256=digest(F/'avl-span1280-capacity-ks-preparation.json'),root=root.relative_to(ROOT).as_posix(),runs=runs,
            completed=completed,stop_reason=stop,stage_memory_before=stage_memory,seconds=time.monotonic()-start,total_bytes=total_bytes(),reference_acceptance=False)
        temporary=output.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(output)
    try:
        save()
        for deck in sorted(prep['decks'],key=lambda d:(0 if d['allocation']=='inner-only' else 1,d['case'])):
            before=policy.memory_snapshot();reason=guard.start_reason(before)
            if reason:raise MemoryError(reason)
            if time.monotonic()-start>1600:raise TimeoutError('stage_time_limit')
            source=ROOT/deck['path'];assert digest(source)==deck['sha256']
            directory=root/f'{deck["case"]}-{deck["allocation"]}';directory.mkdir();shutil.copyfile(source,directory/'wing.avl')
            for name in ['out','after']:(directory/name).mkdir()
            commands=['load wing.avl','oper','o','a','r','t','b','','m','v 20','d 1.2','mn 0','','mrf','a a 1','x']
            for name in ['out/','after/']:commands += [f'ft {name}total.txt',f'fs {name}strips.txt',f'fe {name}elements.txt',f'uh {name}']
            commands+=['','quit',''];(directory/'session.txt').write_text('\n'.join(commands))
            begin=time.monotonic();reason=None;peak_group=0;minimum=before['effective_available'];max_swap=before['swap_used'];last_heartbeat=begin
            with (directory/'session.txt').open('rb') as inp,(directory/'console.txt').open('wb') as console:
                proc=subprocess.Popen(['/usr/bin/time','-v','-o','resource.txt',str(exe)],cwd=directory,stdin=inp,stdout=console,stderr=subprocess.STDOUT,env=environment,start_new_session=True)
                try:
                    while proc.poll() is None:
                        snapshot=policy.memory_snapshot();peak_group=max(peak_group,monitor.live_rss(proc.pid));minimum=min(minimum,snapshot['effective_available']);max_swap=max(max_swap,snapshot['swap_used'])
                        reason=guard.stop_reason(time.monotonic()-begin,time.monotonic()-start,peak_group,(directory/'console.txt').stat().st_size,total_bytes(),snapshot,before['swap_used'],stage_memory['swap_used'])
                        if reason:guard.terminate(proc);break
                        if time.monotonic()-last_heartbeat>=45:
                            print(json.dumps(dict(event='running',case=deck['case'],allocation=deck['allocation'],seconds=time.monotonic()-begin,group_rss=peak_group,effective_available=snapshot['effective_available'])),flush=True);last_heartbeat=time.monotonic()
                        time.sleep(.2)
                finally:
                    if proc.poll() is None:guard.terminate(proc)
            after=policy.memory_snapshot();minimum=min(minimum,after['effective_available']);max_swap=max(max_swap,after['swap_used'])
            timing=(directory/'resource.txt').read_text() if (directory/'resource.txt').exists() else ''
            def counter(label):
                match=re.search(re.escape(label)+r':\s*(\d+)',timing);return int(match[1]) if match else None
            peak=counter('Maximum resident set size (kbytes)');peak=peak*1024 if peak is not None else None
            reason=reason or guard.stop_reason(time.monotonic()-begin,time.monotonic()-start,max(peak_group,peak or 0),(directory/'console.txt').stat().st_size,total_bytes(),after,before['swap_used'],stage_memory['swap_used'])
            record=dict(deck,directory=directory.relative_to(ROOT).as_posix(),seconds=time.monotonic()-begin,memory_before=before,memory_after=after,
                minimum_effective_available=minimum,maximum_swap_used=max_swap,peak_resident_bytes=peak,maximum_group_live_bytes=peak_group,
                minor_faults=counter('Minor (reclaiming a frame) page faults'),major_faults=counter('Major (requiring I/O) page faults'),exit_code=proc.returncode,stopped=reason,files={})
            runs.append(record);save();assert proc.returncode==0 and reason is None and peak is not None
            assert not re.search(r'unsupported|unrecognized|not recognized|NaN|Infinity|REJECTED|too many|overflow|increase.*(MAX|NVX)',(directory/'console.txt').read_text(),re.I)
            record['before_after_equal']={n:digest(directory/'out'/n)==digest(directory/'after'/n) for n in NAMES};assert all(record['before_after_equal'].values())
            record['archive_members']={n:dict(sha256=digest(directory/'after'/n),bytes=(directory/'after'/n).stat().st_size) for n in NAMES}
            with zipfile.ZipFile(directory/'after.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
                for n in NAMES:archive.write(directory/'after'/n,arcname=n)
            with zipfile.ZipFile(directory/'after.zip') as archive:
                for n,h in record['archive_members'].items():
                    content=archive.read(n);assert hashlib.sha256(content).hexdigest()==h['sha256'] and len(content)==h['bytes']
            for n in NAMES:
                p=directory/'after'/n;assert p.resolve().is_relative_to(root.resolve());p.unlink()
            for p in directory.rglob('*'):
                if p.is_file():record['files'][p.relative_to(directory).as_posix()]=dict(sha256=digest(p),bytes=p.stat().st_size)
            assert total_bytes()<=128*1024**2;save();print(json.dumps(dict(event='completed',case=deck['case'],allocation=deck['allocation'],seconds=record['seconds'],peak_resident_bytes=peak)),flush=True)
    except Exception as exc:stop=repr(exc);raise
    finally:save(completed=len(runs)==4 and stop is None and all(len(r.get('files',{}))>0 and len(r.get('before_after_equal',{}))==8 and all(r['before_after_equal'].values()) for r in runs))
if __name__=='__main__':main()
