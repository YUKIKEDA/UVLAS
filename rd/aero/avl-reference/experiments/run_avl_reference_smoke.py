"""Isolated AVL 3.52 output smoke runs, not finite-wing acceptance."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import datetime
import hashlib
import json
import re
import shutil
import subprocess
import time

ROOT=rd_paths.ROOT
EXE=ROOT/'build/reference-avl/avl352.exe'
EXPECTED='443520d255408491222a8df9060bd000f78da95a68845a2f6efdbc67e203f07a'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    assert digest(EXE)==EXPECTED
    manifest=ROOT/'rd/aero/finite-wing/records/finite-wing-reference-inputs.json'
    inputs=json.loads(manifest.read_text())
    directory=ROOT/'build/analysis'/('avl-smoke-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f'))
    directory.mkdir()
    runs=[]
    for label,case in [('planar','FW-V01'),('dihedral','FW-V02'),('repeat','FW-V01')]:
        run=directory/label;run.mkdir()
        deck=next(d for c in inputs['cases'] if c['id']==case for d in c['decks'] if d['m']==4 and d['n']==12)
        source=rd_paths.RECORDS/deck['path'];assert digest(source)==deck['sha256']
        shutil.copyfile(source,run/'wing.avl')
        session='\n'.join(['load wing.avl','oper','o','a','r','t','v','b','','m','v 20','d 1.2','mn 0','',
            'a a 1','x','mrf','ft total.txt','fs strips.txt','fsb strips-body.txt','fe elements.txt','','quit',''])
        (run/'session.txt').write_text(session,encoding='ascii')
        begin=time.monotonic();stopped=None
        with (run/'session.txt').open('rb') as stdin,(run/'console.txt').open('wb') as stdout:
            proc=subprocess.Popen([str(EXE)],cwd=run,stdin=stdin,stdout=stdout,stderr=subprocess.STDOUT,
                creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
            while proc.poll() is None:
                if time.monotonic()-begin>20:stopped='time_limit'
                if (run/'console.txt').stat().st_size>8*1024**2:stopped='log_limit'
                if sum(p.stat().st_size for p in run.iterdir() if p.is_file())>32*1024**2:stopped='output_limit'
                if stopped:
                    proc.kill();proc.wait(timeout=5);break
                time.sleep(.025)
        log=(run/'console.txt').read_text(encoding='latin1')
        outputs={}
        for name in ['total.txt','strips.txt','strips-body.txt','elements.txt']:
            p=run/name
            outputs[name]=dict(exists=p.exists(),bytes=p.stat().st_size if p.exists() else 0,
                sha256=digest(p) if p.exists() else None)
        checks=dict(normal_exit=proc.returncode==0 and stopped is None,
            version='Version  3.52' in log or bool(re.search(r'Version\s+3\.52',log)),
            all_outputs=all(v['exists'] and v['bytes']>0 for v in outputs.values()),
            mrf='machine readable format' in log.lower(),
            no_command_failure=not bool(re.search(r'not recognized|unrecognized|forrtl:|input file.*not found|NaN|Infinity',log,re.I)))
        runs.append(dict(label=label,case=case,directory=run.relative_to(ROOT).as_posix(),
            seconds=time.monotonic()-begin,exit_code=proc.returncode,stopped=stopped,
            input_sha256=digest(run/'wing.avl'),session_sha256=digest(run/'session.txt'),console_sha256=digest(run/'console.txt'),
            outputs=outputs,checks=checks))
        print(json.dumps(dict(label=label,checks=checks,seconds=runs[-1]['seconds'])),flush=True)
    repeated=all(runs[0]['outputs'][name]['sha256']==runs[2]['outputs'][name]['sha256'] for name in runs[0]['outputs'])
    result=dict(protocol='rd/aero/avl-reference/records/avl-reference-execution-protocol.md',
        source_sha256=digest(Path(__file__)),input_manifest_sha256=digest(manifest),
        executable=dict(url='https://web.mit.edu/drela/Public/web/avl/avl352.exe',sha256=EXPECTED,bytes=EXE.stat().st_size),
        runs=runs,repeat_outputs_equal=repeated,
        transport_checks_pass=all(all(r['checks'].values()) for r in runs) and repeated,
        limitations=['transport checks only; semantic parsing and option audit follow',
            'no core-zero source build','no quantitative residual or conservative element force reconstruction yet'])
    output=ROOT/'rd/aero/avl-reference/records/avl-reference-execution-result.json'
    if output.exists(): raise FileExistsError('Preserve earlier reference results; choose a new result name')
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')

if __name__=='__main__':main()
