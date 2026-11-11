"""WSL staged native symmetry runs. Fine runs require a completed equivalence gate."""
from pathlib import Path
import datetime,hashlib,json,os,re,shutil,signal,subprocess,sys,time
import run_avl_span320 as resource

ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    fine=True;stage='fine' if fine else 'regression'
    output=FOLDER/f'avl-symmetry-pi-{stage}-runs.json'
    if output.exists():raise FileExistsError(output)
    prep_path=FOLDER/'avl-symmetry-pi-fine-preparation.json';prep=json.loads(prep_path.read_text())
    gate_path=FOLDER/'avl-symmetry-pi-equivalence-result.json'
    if fine:
        gate=json.loads(gate_path.read_text());assert gate['all_criteria_pass'] and gate['completed']
        assert gate['preparation_sha256']==prep['parent_preparation_sha256']
    root=ROOT/'build/analysis'/('avl-symmetry-pi-'+stage+'-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S'));root.mkdir()
    measured=ROOT/'build/reference-avl/source-symmetry8000-spacing-pi64/bin/avl'
    parent=ROOT/'build/reference-avl/source-core0-grid8000-spacing-pi64/bin/avl'
    selections=[]
    for deck in prep['decks']:
        if fine:
            if deck['physical_n'] in [80,640] or deck.get('fine_role')=='chord_anchor':selections.append((deck,'half_fine',measured,True,ROOT/deck['path']))
        elif deck['physical_n']<640 and deck['physical_n']!=80:
            selections.append((deck,'full_control',measured,True,ROOT/deck['full_input']))
            selections.extend([(deck,'half_baseline',parent,False,ROOT/deck['path']),(deck,'half_measured',measured,True,ROOT/deck['path'])])
    spacing=json.loads((FOLDER/'avl-spacing-pi-result.json').read_text());assert spacing['all_criteria_pass']
    if fine:selections.sort(key=lambda v:v[0]['physical_n'])
    start=time.monotonic();runs=[];stop=None;environment=os.environ.copy();environment['LD_LIBRARY_PATH']=str(ROOT/'build/reference-avl/toolchain/lib')
    try:
        for deck,kind,exe,export,input_path in selections:
            free=resource.available()
            if free<6*1024**3:raise MemoryError('less_than_6GiB_available')
            directory=root/f'{kind}-{deck["case"]}-M{deck["m"]}-N{deck["physical_n"]}';directory.mkdir()
            shutil.copyfile(input_path,directory/'wing.avl')
            for name in ['out','after']:(directory/name).mkdir()
            commands=['load wing.avl','oper','o','a','r','t','b','','m','v 20','d 1.2','mn 0','','mrf','a a 1','x']
            for prefix in ['out/','after/']:
                commands += [f'ft {prefix}total.txt',f'fs {prefix}strips.txt',f'fe {prefix}elements.txt']
                if export:commands += [f'uh {prefix}']
            commands+=['','quit',''];(directory/'session.txt').write_text('\n'.join(commands))
            begin=time.monotonic();reason=None;peak_group=0
            with (directory/'session.txt').open('rb') as inp,(directory/'console.txt').open('wb') as console:
                proc=subprocess.Popen(['/usr/bin/time','-v','-o','resource.txt',str(exe)],cwd=directory,stdin=inp,
                    stdout=console,stderr=subprocess.STDOUT,env=environment,start_new_session=True)
                while proc.poll() is None:
                    peak_group=max(peak_group,resource.live_rss(proc.pid))
                    if time.monotonic()-begin>300:reason='per_process_time_limit'
                    if time.monotonic()-start>(600 if fine else 900):reason='stage_time_limit'
                    if peak_group>4*1024**3:reason='rss_limit'
                    if (directory/'console.txt').stat().st_size>8*1024**2:reason='console_limit'
                    if reason:os.killpg(proc.pid,signal.SIGKILL);proc.wait();break
                    time.sleep(.2)
            text=(directory/'resource.txt').read_text() if (directory/'resource.txt').exists() else ''
            match=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',text)
            rss=int(match[1])*1024 if match else None
            record=dict(case=deck['case'],m=deck['m'],physical_n=deck['physical_n'],kind=kind,
                comparison_directory=deck.get('comparison_directory'),directory=directory.relative_to(ROOT).as_posix(),
                native_unknowns=deck['m']*deck['physical_n'] if kind=='full_control' else deck['m']*deck['native_n'],
                executable_sha256=digest(exe),seconds=time.monotonic()-begin,exit_code=proc.returncode,stopped=reason,
                peak_resident_bytes=rss,maximum_group_live_bytes=peak_group,available_memory_before=free,files={})
            runs.append(record)
            for p in directory.rglob('*'):
                if p.is_file():record['files'][p.relative_to(directory).as_posix()]=dict(sha256=digest(p),bytes=p.stat().st_size)
            assert proc.returncode==0 and reason is None,record
            assert rss is not None and rss<4*1024**3
            console=(directory/'console.txt').read_text()
            assert not re.search(r'unsupported|unrecognized|not recognized|NaN|Infinity|REJECTED',console,re.I)
            names=['total.txt','strips.txt','elements.txt']
            if export:names+=['uv-rows.txt','uv-forces.txt','uv-strips.txt','uh-meta.txt','uh-guard-tests.txt']
            record['before_after_equal']={name:digest(directory/'out'/name)==digest(directory/'after'/name) for name in names}
            assert all(record['before_after_equal'].values())
            if kind=='full_control':
                previous=ROOT/deck['comparison_directory']
                record['legacy_outputs_equal']={name:digest(directory/'out'/name)==digest(previous/name) for name in names if name.startswith('uv-') or name in ['total.txt','strips.txt','elements.txt']}
            if kind=='half_measured':
                base=next(r for r in runs if (r['kind'],r['case'],r['m'],r['physical_n'])==('half_baseline',deck['case'],deck['m'],deck['physical_n']))
                record['parent_half_outputs_equal']={name:digest(directory/'out'/name)==digest(ROOT/base['directory']/'out'/name) for name in ['total.txt','strips.txt','elements.txt']}
                assert all(record['parent_half_outputs_equal'].values())
            assert sum(p.stat().st_size for p in root.rglob('*') if p.is_file())<384*1024**2
            print(json.dumps({k:record[k] for k in ['kind','case','m','physical_n','seconds','peak_resident_bytes']}),flush=True)
    except Exception as exc:stop=repr(exc);raise
    finally:
        result=dict(protocol='docs/verification/avl-symmetry-fine-observation-addendum.md',protocol_sha256=digest(FOLDER/'avl-symmetry-fine-observation-addendum.md'),
            source_sha256=digest(Path(__file__)),resource_runner_sha256=digest(Path(resource.__file__)),preparation_sha256=digest(prep_path),
            spacing_result_sha256=digest(FOLDER/'avl-spacing-pi-result.json'),
            gate_sha256=digest(gate_path) if fine else None,stage=stage,runs=runs,
            completed=len(runs)==(16 if fine else 18) and stop is None,stop_reason=stop,
            seconds=time.monotonic()-start,total_bytes=sum(p.stat().st_size for p in root.rglob('*') if p.is_file()))
        output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
