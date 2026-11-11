"""Frozen-grid reference runs; common LU reused for five prescribed angles."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    output=ROOT/'docs/verification/avl-fixed-band-grid-runs.json'
    if output.exists():raise FileExistsError('preserve earlier grid experiment')
    start=time.monotonic();manifest=ROOT/'docs/verification/finite-wing-reference-inputs.json'
    fixtures=json.loads(manifest.read_text())
    parent_path=ROOT/'docs/verification/avl-stable-reference-runs.json'
    previous=json.loads(parent_path.read_text())
    existing=[r for r in previous['runs'] if r['variant']=='core0']
    assert len(existing)==20
    root=ROOT/'build/analysis'/('avl-fixed-bands-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S'));root.mkdir()
    exe=ROOT/'build/reference-avl/source-core0-grid/bin/avl'
    env=os.environ.copy();env['LD_LIBRARY_PATH']=str(ROOT/'build/reference-avl/toolchain/lib')
    runs=[dict(r,reused=True) for r in existing];processes=[];probe=[];stop_reason=None
    def execute(case,deck,label):
        directory=root/label;directory.mkdir()
        src=ROOT/'docs/verification'/deck['path'];assert digest(src)==deck['sha256']
        shutil.copyfile(src,directory/'wing.avl')
        commands=['load wing.avl','oper','o','a','r','t','b','','m','v 20','d 1.2','mn 0','','mrf']
        for alpha in [-2,-1,0,1,2]:
            prefix=f'a{alpha}/';(directory/prefix).mkdir()
            commands += [f'a a {alpha}','x',f'ft {prefix}total.txt',f'fs {prefix}strips.txt',f'fe {prefix}elements.txt',
                f'ug {prefix}',f'ft {prefix}after-total.txt']
        commands+=['','quit',''];(directory/'session.txt').write_text('\n'.join(commands))
        begin=time.monotonic();reason=None
        with (directory/'session.txt').open('rb') as inp,(directory/'console.txt').open('wb') as out:
            proc=subprocess.Popen(['/usr/bin/time','-v','-o','resource.txt',str(exe)],cwd=directory,
                stdin=inp,stdout=out,stderr=subprocess.STDOUT,env=env,start_new_session=True)
            while proc.poll() is None:
                if time.monotonic()-begin>300:reason='grid_time_limit'
                if (directory/'console.txt').stat().st_size>8*1024**2:reason='console_size_limit'
                if reason:
                    os.killpg(proc.pid,signal.SIGKILL);proc.wait();break
                time.sleep(.1)
        resource=(directory/'resource.txt').read_text() if (directory/'resource.txt').exists() else ''
        peak=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',resource)
        rss=int(peak[1])*1024 if peak else None
        process=dict(case=case['id'],m=deck['m'],n=deck['n'],directory=directory.relative_to(ROOT).as_posix(),
            seconds=time.monotonic()-begin,exit_code=proc.returncode,stopped=reason,peak_resident_bytes=rss,
            input_sha256=digest(directory/'wing.avl'),session_sha256=digest(directory/'session.txt'),
            console_sha256=digest(directory/'console.txt'),resource_sha256=digest(directory/'resource.txt') if resource else None)
        processes.append(process)
        if reason or proc.returncode!=0:raise RuntimeError('AVL grid failure:'+str(reason))
        if rss is None or rss>2*1024**3:raise RuntimeError('memory_limit')
        log=(directory/'console.txt').read_text()
        if re.search(r'unsupported|unrecognized|not recognized|forrtl:|NaN|Infinity',log,re.I):raise RuntimeError('invalid_console')
        records=[]
        for alpha in [-2,-1,0,1,2]:
            d=directory/f'a{alpha}'
            names=['total.txt','strips.txt','elements.txt','after-total.txt','uv-rows.txt','uv-forces.txt','uv-strips.txt']
            files={name:dict(sha256=digest(d/name),bytes=(d/name).stat().st_size) for name in names}
            records.append(dict(variant='core0-grid',case=case['id'],m=deck['m'],n=deck['n'],alpha=alpha,
                directory=d.relative_to(ROOT).as_posix(),reused=False,executable_sha256=digest(exe),
                after_export_equal=files['total.txt']['sha256']==files['after-total.txt']['sha256'],files=files))
        return records
    try:
        case=fixtures['cases'][0];deck=next(d for d in case['decks'] if d['m']==4 and d['n']==12)
        probe=execute(case,deck,'probe-'+case['id'])
        for p in probe:
            old=next(r for r in existing if r['case']==p['case'] and r['m']==p['m'] and r['n']==p['n'] and r['alpha']==p['alpha'])
            p['previous_outputs_equal']=all(p['files'][name]['sha256']==old['files'][name]['sha256'] for name in
                ['total.txt','strips.txt','elements.txt','uv-rows.txt','uv-forces.txt'])
        if not all(p['previous_outputs_equal'] and p['after_export_equal'] for p in probe):raise RuntimeError('coarse_export_regression')
        print(json.dumps(dict(probe='passed',angles=5)),flush=True)
        for case in fixtures['cases']:
            wing_begin=time.monotonic()
            for deck in case['decks']:
                if (deck['m'],deck['n']) in [(4,12),(8,24)]:continue
                if time.monotonic()-wing_begin>1800:raise TimeoutError('wing_runtime_limit')
                records=execute(case,deck,f'{case["id"]}-M{deck["m"]}-N{deck["n"]}')
                runs.extend(records)
                (root/'progress.json').write_text(json.dumps(dict(completed=len(runs),last=records[-1]['directory']))+'\n')
                size=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
                if size>2*1024**3:raise RuntimeError('output_size_limit')
                print(json.dumps(dict(case=case['id'],m=deck['m'],n=deck['n'],completed=len(runs),seconds=processes[-1]['seconds'])),flush=True)
    except Exception as exc:
        stop_reason=repr(exc);raise
    finally:
        result=dict(protocol='docs/verification/avl-fixed-band-convergence-experiment.md',source_sha256=digest(Path(__file__)),
            input_manifest_sha256=digest(manifest),parent_result_sha256=digest(parent_path),
            preparation_sha256=digest(ROOT/'docs/verification/avl-grid-export-preparation.json'),
            executable_sha256=digest(exe),reference_id='AVL3.52/DP/zero-core-exterior-stable-v1',
            probe=probe,processes=processes,runs=runs,completed=len(runs)==100 and stop_reason is None,
            stop_reason=stop_reason,seconds=time.monotonic()-start)
        output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
