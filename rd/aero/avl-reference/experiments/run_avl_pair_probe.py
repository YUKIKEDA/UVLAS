import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess
import time

ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    start=time.monotonic();parent=ROOT/'rd/aero/avl-reference/records/avl-source-measurement-corrected-runs.json'
    previous=json.loads(parent.read_text());selected=[r for r in previous['runs'] if r['variant']=='core0' and r['alpha']==1]
    assert len(selected)==4
    root=ROOT/'build/analysis'/('avl-pair-probe-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S'));root.mkdir()
    exe=ROOT/'build/reference-avl/source-core0-pairprobe/bin/avl'
    env=os.environ.copy();env['LD_LIBRARY_PATH']=str(ROOT/'build/reference-avl/toolchain/lib')
    runs=[]
    for old in selected:
        run=root/f'{old["case"]}-M{old["m"]}-N{old["n"]}';run.mkdir()
        source=ROOT/old['directory']
        for name in ['wing.avl','session.txt']:shutil.copyfile(source/name,run/name)
        with (run/'session.txt').open('rb') as inp,(run/'console.txt').open('wb') as out:
            result=subprocess.run(['/usr/bin/time','-v','-o','resource.txt',str(exe)],cwd=run,
                stdin=inp,stdout=out,stderr=subprocess.STDOUT,env=env,timeout=20,check=True)
        resource=(run/'resource.txt').read_text();rss=int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',resource)[1])*1024
        names=['total.txt','strips.txt','elements.txt','uv-rows.txt','uv-matrix.txt','uv-forces.txt']
        unchanged=all(digest(source/n)==digest(run/n) for n in names)
        files={p.name:dict(sha256=digest(p),bytes=p.stat().st_size) for p in run.iterdir() if p.is_file()}
        runs.append(dict(case=old['case'],m=old['m'],n=old['n'],directory=run.relative_to(ROOT).as_posix(),
            previous_directory=old['directory'],original_outputs_unchanged=unchanged,peak_resident_bytes=rss,files=files))
    result=dict(protocol='rd/aero/avl-reference/records/avl-collinear-kernel-experiment.md',source_sha256=digest(Path(__file__)),
        parent_result_sha256=digest(parent),executable_sha256=digest(exe),runs=runs,seconds=time.monotonic()-start)
    (ROOT/'rd/aero/avl-reference/records/avl-pair-probe-runs.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(runs=4,unchanged=all(r['original_outputs_unchanged'] for r in runs),seconds=result['seconds'])))

if __name__=='__main__':main()
