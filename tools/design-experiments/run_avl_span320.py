"""WSL runner for four capacity regressions then two predeclared fine grids."""
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,signal,subprocess,time

ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def available():
    return int(re.search(r'MemAvailable:\s*(\d+)',Path('/proc/meminfo').read_text())[1])*1024
def live_rss(pid):
    todo=[pid];seen=set();total=0
    while todo:
        p=todo.pop()
        if p in seen:continue
        seen.add(p)
        try:
            text=Path(f'/proc/{p}/status').read_text();match=re.search(r'VmRSS:\s*(\d+)',text)
            if match:total+=int(match[1])*1024
            todo.extend(int(x) for x in Path(f'/proc/{p}/task/{p}/children').read_text().split())
        except FileNotFoundError:pass
    return total

def main():
    result_path=FOLDER/'avl-span320-runs.json'
    if result_path.exists():raise FileExistsError(result_path)
    prep_path=FOLDER/'avl-span320-preparation.json';prep=json.loads(prep_path.read_text())
    old_path=FOLDER/'avl-band-aligned-grid-runs.json';old=json.loads(old_path.read_text())
    exe=ROOT/'build/reference-avl/source-core0-grid8000/bin/avl'
    root=ROOT/'build/analysis'/('avl-span320-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S'));root.mkdir()
    environment=os.environ.copy();environment['LD_LIBRARY_PATH']=str(ROOT/'build/reference-avl/toolchain/lib')
    start=time.monotonic();runs=[];stop=None
    selections=[]
    for case in ['FW-V01','FW-V02']:
        for m,n in [(6,40),(24,160)]:
            r=next(r for r in old['runs'] if (r['case'],r['m'],r['n'],r['alpha'])==(case,m,n,1))
            selections.append(dict(case=case,m=m,n=n,kind='regression',previous=r['directory'],input=str(ROOT/r['directory']).rsplit('/a1',1)[0]+'/wing.avl'))
    for deck in prep['decks']:selections.append(dict(case=deck['case'],m=24,n=320,kind='new',input=str(ROOT/deck['path']),input_sha256=deck['sha256']))
    try:
        for item in selections:
            if item['kind']=='new':assert sum(r['kind']=='regression' and all(r['previous_outputs_equal'].values()) for r in runs)==4
            free=available()
            if free<6*1024**3:raise MemoryError('less_than_6GiB_available')
            directory=root/f'{item["kind"]}-{item["case"]}-M{item["m"]}-N{item["n"]}';directory.mkdir()
            shutil.copyfile(item['input'],directory/'wing.avl')
            if 'input_sha256' in item:assert digest(directory/'wing.avl')==item['input_sha256']
            for name in ['out','after']:(directory/name).mkdir()
            commands=['load wing.avl','oper','o','a','r','t','b','','m','v 20','d 1.2','mn 0','','mrf','a a 1','x']
            for prefix in ['out/','after/']:commands += [f'ft {prefix}total.txt',f'fs {prefix}strips.txt',f'fe {prefix}elements.txt',f'ug {prefix}']
            commands+=['','quit',''];(directory/'session.txt').write_text('\n'.join(commands))
            begin=time.monotonic();reason=None;max_live=0
            with (directory/'session.txt').open('rb') as inp,(directory/'console.txt').open('wb') as console:
                proc=subprocess.Popen(['/usr/bin/time','-v','-o','resource.txt',str(exe)],cwd=directory,stdin=inp,
                    stdout=console,stderr=subprocess.STDOUT,env=environment,start_new_session=True)
                while proc.poll() is None:
                    max_live=max(max_live,live_rss(proc.pid))
                    if time.monotonic()-begin>300:reason='per_run_time_limit'
                    if time.monotonic()-start>900:reason='total_time_limit'
                    if max_live>4*1024**3:reason='rss_limit'
                    if (directory/'console.txt').stat().st_size>8*1024**2:reason='console_limit'
                    if reason:os.killpg(proc.pid,signal.SIGKILL);proc.wait();break
                    time.sleep(.2)
            resource=(directory/'resource.txt').read_text() if (directory/'resource.txt').exists() else ''
            peak=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',resource)
            peak_bytes=int(peak[1])*1024 if peak else None
            record={k:item[k] for k in ['case','m','n','kind']};record.update(directory=directory.relative_to(ROOT).as_posix(),
                alpha=1,seconds=time.monotonic()-begin,exit_code=proc.returncode,stopped=reason,
                available_memory_before=free,peak_resident_bytes=peak_bytes,maximum_group_live_bytes=max_live,files={})
            runs.append(record)
            for p in directory.rglob('*'):
                if p.is_file():record['files'][p.relative_to(directory).as_posix()]=dict(sha256=digest(p),bytes=p.stat().st_size)
            assert proc.returncode==0 and reason is None,record
            assert peak_bytes is not None and peak_bytes<4*1024**3
            assert not re.search(r'unsupported|unrecognized|not recognized|NaN|Infinity',(directory/'console.txt').read_text(),re.I)
            names=['total.txt','strips.txt','elements.txt','uv-rows.txt','uv-forces.txt','uv-strips.txt']
            record['before_after_equal']={name:digest(directory/'out'/name)==digest(directory/'after'/name) for name in names}
            assert all(record['before_after_equal'].values())
            if item['kind']=='regression':
                record['previous_directory']=item['previous'];previous=ROOT/item['previous']
                record['previous_outputs_equal']={name:digest(previous/name)==digest(directory/'out'/name) for name in names}
                assert all(record['previous_outputs_equal'].values()),record['previous_outputs_equal']
            if sum(p.stat().st_size for p in root.rglob('*') if p.is_file())>128*1024**2:raise RuntimeError('disk_limit')
            print(json.dumps({k:record[k] for k in ['kind','case','m','n','seconds','peak_resident_bytes']}),flush=True)
    except Exception as exc:stop=repr(exc);raise
    finally:
        result=dict(protocol='docs/verification/avl-span320-reference-experiment.md',protocol_sha256=digest(FOLDER/'avl-span320-reference-experiment.md'),
            source_sha256=digest(Path(__file__)),preparation_sha256=digest(prep_path),parent_runs_sha256=digest(old_path),
            executable_sha256=digest(exe),runs=runs,completed=len(runs)==6 and stop is None,stop_reason=stop,
            seconds=time.monotonic()-start,total_bytes=sum(p.stat().st_size for p in root.rglob('*') if p.is_file()))
        result_path.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
