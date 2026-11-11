"""Run external AVL variants under WSL; numerical analysis is separate."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    manifest=ROOT/'docs/verification/finite-wing-reference-inputs.json'
    fixtures=json.loads(manifest.read_text())
    root=ROOT/'build/analysis'/('avl-stable-reference-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S'))
    root.mkdir()
    env=os.environ.copy();env['LD_LIBRARY_PATH']=str(ROOT/'build/reference-avl/toolchain/lib')
    initial_path=ROOT/'docs/verification/avl-source-measurement-runs.json'
    initial=json.loads(initial_path.read_text())
    runs=[r for r in initial['runs'] if r['variant']=='baseline']
    try:
        for variant in ['measured','core0']:
            exe=ROOT/f'build/reference-avl/source-{variant}-stable/bin/avl'
            for case in fixtures['cases']:
                for m,n in [(4,12),(8,24)]:
                    deck=next(d for d in case['decks'] if d['m']==m and d['n']==n)
                    for alpha in [-2,-1,0,1,2]:
                        if time.monotonic()-START>300:raise TimeoutError('global_runtime_limit')
                        label=f'{variant}-{case["id"]}-M{m}-N{n}-a{alpha}'
                        run=root/label;run.mkdir()
                        source=ROOT/'docs/verification'/deck['path']
                        assert digest(source)==deck['sha256']
                        shutil.copyfile(source,run/'wing.avl')
                        commands=['load wing.avl','oper','o','a','r','t','b','','m','v 20','d 1.2','mn 0','',
                            f'a a {alpha}','x','mrf','ft before-total.txt','fs before-strips.txt','fe before-elements.txt']
                        if variant!='baseline':commands.append('ud')
                        commands+=['ft total.txt','fs strips.txt','fe elements.txt','','quit','']
                        (run/'session.txt').write_text('\n'.join(commands))
                        begin=time.monotonic()
                        with (run/'session.txt').open('rb') as inp,(run/'console.txt').open('wb') as out:
                            proc=subprocess.Popen(['/usr/bin/time','-v','-o','resource.txt',str(exe)],
                                cwd=run,stdin=inp,stdout=out,stderr=subprocess.STDOUT,env=env)
                            try:proc.wait(timeout=30)
                            except subprocess.TimeoutExpired:
                                proc.kill();proc.wait();raise
                        resource=(run/'resource.txt').read_text()
                        rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',resource)
                        files={p.name:dict(sha256=digest(p),bytes=p.stat().st_size) for p in run.iterdir() if p.is_file()}
                        item=dict(variant=variant,case=case['id'],m=m,n=n,alpha=alpha,
                            directory=run.relative_to(ROOT).as_posix(),executable_sha256=digest(exe),
                            seconds=time.monotonic()-begin,exit_code=proc.returncode,
                            peak_resident_bytes=int(rss[1])*1024 if rss else None,files=files)
                        runs.append(item)
                        if proc.returncode!=0:raise RuntimeError('avl_process_failure:'+label)
                        if item['peak_resident_bytes'] is None or item['peak_resident_bytes']>2*1024**3:raise RuntimeError('memory_limit')
                        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file())>256*1024**2:raise RuntimeError('output_limit')
            print(json.dumps(dict(variant=variant,completed=len(runs))),flush=True)
    finally:
        result=dict(protocol='docs/verification/avl-source-measurement-experiment.md',
            source_sha256=digest(Path(__file__)),input_manifest_sha256=digest(manifest),
            preparation_sha256=digest(ROOT/'docs/verification/avl-stable-reference-preparation.json'),
            build_script_sha256=digest(ROOT/'build/reference-avl/build-stable-reference.sh'),
            compiler=(ROOT/'build/reference-avl/compiler-version.txt').read_text(),
            toolchain_lock_sha256=digest(ROOT/'build/reference-avl/toolchain-explicit.txt'),
            initial_result_sha256=digest(initial_path),baseline_runs_reused=20,
            runs=runs,seconds=time.monotonic()-START,completed=len(runs)==60)
        output=ROOT/'docs/verification/avl-stable-reference-runs.json'
        if output.exists():raise FileExistsError('preserve earlier source measurement runs')
        output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
