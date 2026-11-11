"""Execute only the eight predeclared reference sensitivity systems in WSL."""
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
    output=ROOT/'docs/verification/avl-solve-sensitivity-runs.json'
    if output.exists():raise FileExistsError(output)
    manifests={f:ROOT/'docs/verification'/name for f,name in [
        ('original','avl-fixed-band-grid-runs.json'),('aligned','avl-band-aligned-grid-runs.json')]}
    parents={f:json.loads(p.read_text()) for f,p in manifests.items()}
    selected=([('original','FW-V01',4,12),('original','FW-V02',8,24),
        ('original','FW-V01',16,96),('original','FW-V01',32,96)]+
        [('aligned',case,m,160) for case in ['FW-V01','FW-V02'] for m in [12,24]])
    root=ROOT/'build/analysis'/('avl-sensitivity-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S'));root.mkdir()
    exe=ROOT/'build/reference-avl/source-core0-sensitivity-corrected/bin/avl'
    env=os.environ.copy();env['LD_LIBRARY_PATH']=str(ROOT/'build/reference-avl/toolchain/lib')
    begin=time.monotonic();runs=[];stop=None
    try:
        for family,case,m,n in selected:
            old=next(r for r in parents[family]['runs'] if (r['case'],r['m'],r['n'],r['alpha'])==(case,m,n,1))
            previous=ROOT/old['directory'];directory=root/f'{family}-{case}-M{m}-N{n}';directory.mkdir()
            source=previous.parent/'wing.avl'
            if not source.exists():source=previous/'wing.avl'
            shutil.copyfile(source,directory/'wing.avl')
            (directory/'out').mkdir();(directory/'after').mkdir()
            cmds=['load wing.avl','oper','o','a','r','t','b','','m','v 20','d 1.2','mn 0','','mrf','a a 1','x']
            for prefix in ['out/','after/']:
                cmds += [f'ft {prefix}total.txt',f'fs {prefix}strips.txt',f'fe {prefix}elements.txt',f'ug {prefix}']
                if prefix=='out/':cmds += ['us out/']
            cmds += ['','quit',''];(directory/'session.txt').write_text('\n'.join(cmds))
            start=time.monotonic();reason=None
            with (directory/'session.txt').open('rb') as inp,(directory/'console.txt').open('wb') as console:
                proc=subprocess.Popen(['/usr/bin/time','-v','-o','resource.txt',str(exe)],cwd=directory,
                    stdin=inp,stdout=console,stderr=subprocess.STDOUT,env=env,start_new_session=True)
                while proc.poll() is None:
                    if time.monotonic()-start>180:reason='per_system_time_limit'
                    if time.monotonic()-begin>900:reason='total_time_limit'
                    if (directory/'console.txt').stat().st_size>8*1024**2:reason='console_limit'
                    if reason:os.killpg(proc.pid,signal.SIGKILL);proc.wait();break
                    time.sleep(.1)
            resource=(directory/'resource.txt').read_text() if (directory/'resource.txt').exists() else ''
            peak=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',resource)
            rss=int(peak[1])*1024 if peak else None
            record=dict(family=family,case=case,m=m,n=n,original_directory=old['directory'],
                directory=directory.relative_to(ROOT).as_posix(),seconds=time.monotonic()-start,
                exit_code=proc.returncode,stopped=reason,peak_resident_bytes=rss,files={})
            runs.append(record)
            for p in directory.rglob('*'):
                if p.is_file():record['files'][p.relative_to(directory).as_posix()]=dict(sha256=digest(p),bytes=p.stat().st_size)
            assert proc.returncode==0 and reason is None,record
            assert rss is not None and rss<=2*1024**3,record
            assert not re.search(r'unsupported|unrecognized|not recognized|NaN|Infinity', (directory/'console.txt').read_text(),re.I)
            names=['total.txt','strips.txt','elements.txt','uv-rows.txt','uv-forces.txt']
            if (previous/'uv-strips.txt').exists():names+=['uv-strips.txt']
            record['unchanged']={name:digest(previous/name)==digest(directory/'out'/name)==digest(directory/'after'/name) for name in names}
            assert all(record['unchanged'].values()),record['unchanged']
            assert sum(p.stat().st_size for p in root.rglob('*') if p.is_file())<=256*1024**2
            print(json.dumps({k:record[k] for k in ['family','case','m','n','seconds','peak_resident_bytes']}),flush=True)
    except Exception as exc:
        stop=repr(exc);raise
    finally:
        result=dict(protocol='docs/verification/avl-solve-sensitivity-experiment.md',source_sha256=digest(Path(__file__)),
            preparation_sha256=digest(ROOT/'docs/verification/avl-sensitivity-preparation-corrected.json'),
            parents={f:dict(path=p.relative_to(ROOT).as_posix(),sha256=digest(p)) for f,p in manifests.items()},
            executable_sha256=digest(exe),runs=runs,completed=len(runs)==8 and stop is None,stop_reason=stop,
            seconds=time.monotonic()-begin,total_bytes=sum(p.stat().st_size for p in root.rglob('*') if p.is_file()))
        output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
