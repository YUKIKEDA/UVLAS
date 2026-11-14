"""Eight capacity-only regressions. Deliberately contains no N1280 solve."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import datetime,json,os,re,shutil,signal,subprocess,time
import avl_reference_resources as policy
import run_avl_span320 as monitor

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS;GIB=1024**3
digest=policy.digest

def main():
    output=FOLDER/'avl-span1280-capacity-ks-runs.json'
    if output.exists():raise FileExistsError(output)
    prep_path=FOLDER/'avl-span1280-capacity-ks-preparation.json';prep=json.loads(prep_path.read_text())
    protocol=ROOT/prep['protocol'];assert digest(protocol)==prep['protocol_sha256']
    for name,hashes in prep['source_files'].items():assert digest(ROOT/prep['target']/'src'/name)==hashes['capacity']
    exe=ROOT/prep['target']/'bin/avl';assert exe.exists()
    gate_path=FOLDER/'avl-symmetry-pi-equivalence-result.json';gate=json.loads(gate_path.read_text());assert gate['all_criteria_pass'] and digest(gate_path)==prep['gate_sha256']
    root=ROOT/'build/analysis'/('avl-capacity16000-ks1000-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S'));root.mkdir()
    selections=sorted(prep['inputs'],key=lambda r:(r['native_unknowns'],r['case']))
    environment=os.environ.copy();environment['LD_LIBRARY_PATH']=str(ROOT/'build/reference-avl/toolchain/lib')
    start=time.monotonic();runs=[];stop=None
    def total_bytes():return sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    try:
        for item in selections:
            before=policy.memory_snapshot()
            if before['effective_available']<8*GIB:raise MemoryError('less_than_8GiB_effective_available')
            parent=ROOT/item['parent_directory'];inp=ROOT/item['input'];assert digest(inp)==item['input_sha256']
            directory=root/f'{item["case"]}-M{item["m"]}-N{item["physical_n"]}';directory.mkdir();shutil.copyfile(inp,directory/'wing.avl')
            for folder in ['out','after']:(directory/folder).mkdir()
            commands=['load wing.avl','oper','o','a','r','t','b','','m','v 20','d 1.2','mn 0','','mrf','a a 1','x']
            for folder in ['out/','after/']:commands += [f'ft {folder}total.txt',f'fs {folder}strips.txt',f'fe {folder}elements.txt',f'uh {folder}']
            commands+=['','quit',''];(directory/'session.txt').write_text('\n'.join(commands))
            begin=time.monotonic();reason=None;peak_group=0
            with (directory/'session.txt').open('rb') as source,(directory/'console.txt').open('wb') as console:
                proc=subprocess.Popen(['/usr/bin/time','-v','-o','resource.txt',str(exe)],cwd=directory,stdin=source,stdout=console,stderr=subprocess.STDOUT,env=environment,start_new_session=True)
                try:
                    while proc.poll() is None:
                        peak_group=max(peak_group,monitor.live_rss(proc.pid))
                        if time.monotonic()-begin>300:reason='per_process_time_limit'
                        if time.monotonic()-start>900:reason='stage_time_limit'
                        if peak_group>5*GIB:reason='process_group_rss_limit'
                        if (directory/'console.txt').stat().st_size>8*1024**2:reason='console_limit'
                        if total_bytes()>256*1024**2:reason='disk_limit'
                        if reason:os.killpg(proc.pid,signal.SIGKILL);proc.wait();break
                        time.sleep(.2)
                finally:
                    if proc.poll() is None:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
            timing=(directory/'resource.txt').read_text() if (directory/'resource.txt').exists() else ''
            def counter(label):
                match=re.search(re.escape(label)+r':\s*(\d+)',timing);return int(match[1]) if match else None
            peak=counter('Maximum resident set size (kbytes)');peak=peak*1024 if peak is not None else None
            record=dict(item,directory=directory.relative_to(ROOT).as_posix(),seconds=time.monotonic()-begin,
                memory_before=before,memory_after=policy.memory_snapshot(),peak_resident_bytes=peak,maximum_group_live_bytes=peak_group,
                minor_faults=counter('Minor (reclaiming a frame) page faults'),major_faults=counter('Major (requiring I/O) page faults'),
                exit_code=proc.returncode,stopped=reason,files={})
            runs.append(record)
            for p in directory.rglob('*'):
                if p.is_file():record['files'][p.relative_to(directory).as_posix()]=dict(sha256=digest(p),bytes=p.stat().st_size)
            assert proc.returncode==0 and reason is None and peak is not None and peak<5*GIB
            assert peak_group<5*GIB and record['seconds']<300 and time.monotonic()-start<900 and total_bytes()<256*1024**2
            assert not re.search(r'unsupported|unrecognized|not recognized|NaN|Infinity|REJECTED|too many|increase.*(MAX|NVX)',(directory/'console.txt').read_text(),re.I)
            names=['total.txt','strips.txt','elements.txt','uv-rows.txt','uv-forces.txt','uv-strips.txt','uh-meta.txt','uh-guard-tests.txt']
            if item['native_unknowns']<=256:names+=['uh-matrix.bin']
            record['before_after_equal']={n:digest(directory/'out'/n)==digest(directory/'after'/n) for n in names}
            record['parent_outputs_equal']={n:digest(directory/'out'/n)==digest(parent/n) for n in names}
            record['parent_file_sha256']={n:digest(parent/n) for n in names}
            assert all(record['before_after_equal'].values()) and all(record['parent_outputs_equal'].values())
            print(json.dumps({k:record[k] for k in ['case','m','physical_n','seconds','peak_resident_bytes','major_faults']}),flush=True)
    except Exception as exc:stop=repr(exc);raise
    finally:
        result=dict(protocol=prep['protocol'],protocol_sha256=digest(protocol),source_sha256=digest(Path(__file__)),
            dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [policy,monitor]},preparation_sha256=digest(prep_path),
            resource_audit_sha256=digest(FOLDER/'avl-span1280-resource-audit.json'),executable_sha256=digest(exe),
            build_script_sha256=digest(ROOT/'build/reference-avl/build-symmetry16000-ks1000.sh'),
            build_log_sha256=digest(ROOT/'build/reference-avl/build-symmetry16000-ks1000.log'),build_resource_sha256=digest(ROOT/'build/reference-avl/build-symmetry16000-resource.txt'),
            runs=runs,completed=len(runs)==8 and stop is None,stop_reason=stop,seconds=time.monotonic()-start,total_bytes=total_bytes(),
            no_n1280_run_performed=True,reference_acceptance=False)
        output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
