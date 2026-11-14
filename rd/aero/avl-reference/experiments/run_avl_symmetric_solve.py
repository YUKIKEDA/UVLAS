"""Four fixed-count allocation solves with verified lossless duplicate archiving."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import datetime,hashlib,json,math,os,re,shutil,subprocess,sys,time,zipfile
import avl_reference_resources as policy
import avl_symmetric_solve_watchdog as guard
import avl_span1280_watchdog as termination
import run_avl_span320 as monitor
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;digest=policy.digest
NAMES=['total.txt','strips.txt','elements.txt','uv-rows.txt','uv-forces.txt','uv-strips.txt','uh-meta.txt','uh-guard-tests.txt']
def main():
    stage='main' if '--main' in sys.argv[1:] else 'small'
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
    root=ROOT/'build/analysis'/('avl-symmetric-solve-'+stage+'-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S'));root.mkdir()
    environment=os.environ.copy();environment['LD_LIBRARY_PATH']=str(ROOT/'build/reference-avl/toolchain/lib')
    start=time.monotonic();runs=[];stop=None;stage_memory=policy.memory_snapshot()
    def total_bytes():return sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    def save(completed=False):
        result=dict(protocol=prep['protocol'],protocol_sha256=prep['protocol_sha256'],execution_protocol=prep['protocol'],
            execution_protocol_sha256=prep['protocol_sha256'],source_sha256=digest(__file__),
            dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [policy,guard,termination,monitor]},
            preparation_sha256=digest(prep_path),gate_sha256={n:digest(F/n) for n in gates},executable_sha256=digest(exe),
            stage=stage,build_script_sha256=prep['build_script_sha256'],build_log_sha256=digest(ROOT/'build/reference-avl/build-symmetric-solve-probe.log'),build_resource_sha256=digest(ROOT/'build/reference-avl/build-symmetric-solve-probe-resource.txt'),root=root.relative_to(ROOT).as_posix(),runs=runs,
            completed=completed,stop_reason=stop,stage_memory_before=stage_memory,seconds=time.monotonic()-start,total_bytes=total_bytes(),reference_acceptance=False)
        temporary=output.with_suffix('.tmp');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(output)
    try:
        save()
        for deck in sorted(selections,key=lambda d:(d['case'],{'small':0,'uniform':0,'tip-rich':1,'tip-poor':2}[d['allocation']])):
            names=NAMES+(['uh-matrix.bin'] if deck['native_unknowns']<=256 else [])
            before=policy.memory_snapshot();reason=guard.start_reason(before)
            if reason:raise MemoryError(reason)
            if time.monotonic()-start>guard.LIMITS[stage][1]:raise TimeoutError('stage_time_limit')
            source=ROOT/deck['input'];assert digest(source)==deck['input_sha256']
            directory=root/f'{deck["case"]}-{deck["allocation"]}';directory.mkdir();shutil.copyfile(source,directory/'wing.avl')
            for name in ['out','after','uq']:(directory/name).mkdir()
            commands=['load wing.avl','oper','o','a','r','t','b','','m','v 20','d 1.2','mn 0','','mrf','a a 1','x']
            for name in ['out/','after/']:
                if name=='after/':commands+=['uq uq/']
                commands += [f'ft {name}total.txt',f'fs {name}strips.txt',f'fe {name}elements.txt',f'uh {name}']
            commands+=['','quit',''];(directory/'session.txt').write_text('\n'.join(commands))
            begin=time.monotonic();reason=None;peak_group=0;minimum=before['effective_available'];max_swap=before['swap_used'];last_heartbeat=begin
            with (directory/'session.txt').open('rb') as inp,(directory/'console.txt').open('wb') as console:
                proc=subprocess.Popen(['/usr/bin/time','-v','-o','resource.txt',str(exe)],cwd=directory,stdin=inp,stdout=console,stderr=subprocess.STDOUT,env=environment,start_new_session=True)
                try:
                    while proc.poll() is None:
                        snapshot=policy.memory_snapshot();peak_group=max(peak_group,monitor.live_rss(proc.pid));minimum=min(minimum,snapshot['effective_available']);max_swap=max(max_swap,snapshot['swap_used'])
                        reason=guard.stop_reason(stage,time.monotonic()-begin,time.monotonic()-start,peak_group,(directory/'console.txt').stat().st_size,total_bytes(),snapshot,before['swap_used'],stage_memory['swap_used'])
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
            reason=reason or guard.stop_reason(stage,time.monotonic()-begin,time.monotonic()-start,max(peak_group,peak or 0),(directory/'console.txt').stat().st_size,total_bytes(),after,before['swap_used'],stage_memory['swap_used'])
            record=dict(deck,directory=directory.relative_to(ROOT).as_posix(),seconds=time.monotonic()-begin,memory_before=before,memory_after=after,
                minimum_effective_available=minimum,maximum_swap_used=max_swap,peak_resident_bytes=peak,maximum_group_live_bytes=peak_group,
                minor_faults=counter('Minor (reclaiming a frame) page faults'),major_faults=counter('Major (requiring I/O) page faults'),exit_code=proc.returncode,stopped=reason,files={})
            runs.append(record);save();assert proc.returncode==0 and reason is None and peak is not None
            assert not re.search(r'unsupported|unrecognized|not recognized|NaN|Infinity|REJECTED|too many|overflow|increase.*(MAX|NVX)',(directory/'console.txt').read_text(),re.I)
            record['before_after_equal']={n:digest(directory/'out'/n)==digest(directory/'after'/n) for n in names};assert all(record['before_after_equal'].values())
            record['parent_outputs_equal']={n:digest(directory/'out'/n)==digest(ROOT/deck['raw_directory']/n) for n in names}
            assert all(record['parent_outputs_equal'].values())
            meta={v.split()[0]:int(v.split()[1]) for v in (directory/'uq/uq-meta.txt').read_text().splitlines()}
            assert meta['N']==deck['native_unknowns'] and meta['REAL_BITS']==64 and meta['QUAD_BITS']==128 and meta['QUAD_DIGITS']>=30 and meta['STATE_UNCHANGED']==1
            tests=[list(map(int,line.split())) for line in (directory/'uq/uq-guard-tests.txt').read_text().splitlines()[1:]]
            assert len(tests)==8 and all((status!=0)==bool(expected) for _,status,expected in tests)
            vectors=(directory/'uq/uq-vectors.txt').read_text().splitlines()[1:]
            assert len(vectors)==deck['native_unknowns']
            assert all(len(v.split())==10 and all(math.isfinite(float(x)) for x in v.split()) for v in vectors)
            record['probe_meta']=meta
            record['archive_members']={n:dict(sha256=digest(directory/'after'/n),bytes=(directory/'after'/n).stat().st_size) for n in names}
            with zipfile.ZipFile(directory/'after.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
                for n in names:archive.write(directory/'after'/n,arcname=n)
            with zipfile.ZipFile(directory/'after.zip') as archive:
                for n,h in record['archive_members'].items():
                    content=archive.read(n);assert hashlib.sha256(content).hexdigest()==h['sha256'] and len(content)==h['bytes']
            for n in names:
                p=directory/'after'/n;assert p.resolve().is_relative_to(root.resolve());p.unlink()
            for p in directory.rglob('*'):
                if p.is_file():record['files'][p.relative_to(directory).as_posix()]=dict(sha256=digest(p),bytes=p.stat().st_size)
            assert total_bytes()<=128*1024**2;save();print(json.dumps(dict(event='completed',case=deck['case'],allocation=deck['allocation'],seconds=record['seconds'],peak_resident_bytes=peak)),flush=True)
    except Exception as exc:stop=repr(exc);raise
    finally:save(completed=len(runs)==len(selections) and stop is None and all(len(r.get('files',{}))>0 and len(r.get('before_after_equal',{}))==(9 if r['native_unknowns']<=256 else 8) and all(r['before_after_equal'].values()) for r in runs))
if __name__=='__main__':main()
