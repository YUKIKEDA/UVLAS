"""Isolated AVL 3.52 output smoke runs, not finite-wing acceptance."""
from pathlib import Path
import datetime
import ctypes
import hashlib
import json
import re
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
EXE=ROOT/'build/reference-avl/avl352.exe'
EXPECTED='443520d255408491222a8df9060bd000f78da95a68845a2f6efdbc67e203f07a'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def process_peak(proc):
    class Counters(ctypes.Structure):
        _fields_=[('cb',ctypes.c_ulong),('faults',ctypes.c_ulong)]+[(name,ctypes.c_size_t) for name in
            ['PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage',
             'QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage']]
    counters=Counters();counters.cb=ctypes.sizeof(counters)
    query=ctypes.windll.psapi.GetProcessMemoryInfo
    query.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_ulong]
    query.restype=ctypes.c_int
    if not query(int(proc._handle),ctypes.byref(counters),counters.cb):raise ctypes.WinError()
    return int(counters.PeakWorkingSetSize)

def semantic(run,log):
    def final_flag(label):
        values=re.findall(re.escape(label)+r'\s*:\s*([TF])',log)
        return values[-1] if values else None
    total=(run/'total.txt').read_text()
    values={}
    for line in total.splitlines():
        if '|' in line:
            left,right=line.split('|',1)
            try:values[right.strip()]=[float(v) for v in left.split()]
            except ValueError:pass
    elements=(run/'elements.txt').read_text()
    strips=(run/'strips.txt').read_text()
    counts=[int(values[k][0]) for k in ['# surfaces','# strips','# vortices']]
    options=dict(trailing=final_flag('T rail.leg forces'),viscous=final_flag('V iscous forces'),body=final_flag('B ody forces'))
    velocities=re.findall(r'V  velocity\s*=\s*([0-9.Ee+\-]+)',log)
    densities=re.findall(r'D  air dens\.\s*=\s*([0-9.Ee+\-]+)',log)
    geometric='Geometric axis orientation,  X aft, Z up' in total
    body_axis='Rates,moments about Body Axes, X along geometric X' in log
    element_rows=[line for line in elements.splitlines() if '| I, X, Y, Z, DX, Slope, dCp' in line]
    strip_rows=[line for line in strips.splitlines() if re.match(r'^\s*\d+\s+[-+0-9.]',line) and '|' not in line]
    checks=dict(options=options==dict(trailing='T',viscous='F',body='F'),axes=geometric and body_axis,
        counts=counts==[1,12,48] and len(element_rows)==48 and len(strip_rows)==12,
        reference_quantities=values.get('Sref, Cref, Bref')==[8.4,1.114285714285714,8.],
        origin=values.get('Xref, Yref, Zref')==[.375,0.,0.],
        angle=values.get("Alpha, pb/2V, p'b/2V")==[1.,0.,0.],
        mach=values.get("Mach, rb/2V, r'b/2V")==[0.,0.,0.],
        flow_parameters=bool(velocities and densities and float(velocities[-1])==20. and float(densities[-1])==1.2))
    return dict(options=options,values=values,element_rows=len(element_rows),strip_rows=len(strip_rows),checks=checks)

def main():
    assert digest(EXE)==EXPECTED
    manifest=ROOT/'docs/verification/finite-wing-reference-inputs.json'
    inputs=json.loads(manifest.read_text())
    directory=ROOT/'build/analysis'/('avl-smoke-corrected-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f'))
    directory.mkdir()
    runs=[]
    for label,case in [('planar','FW-V01'),('dihedral','FW-V02'),('repeat','FW-V01')]:
        run=directory/label;run.mkdir()
        deck=next(d for c in inputs['cases'] if c['id']==case for d in c['decks'] if d['m']==4 and d['n']==12)
        source=ROOT/'docs/verification'/deck['path'];assert digest(source)==deck['sha256']
        shutil.copyfile(source,run/'wing.avl')
        session='\n'.join(['load wing.avl','oper','o','a','r','t','b','','m','v 20','d 1.2','mn 0','',
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
        peak=process_peak(proc)
        log=(run/'console.txt').read_text(encoding='latin1')
        meaning=semantic(run,log)
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
            outputs=outputs,checks=checks,semantic=meaning,peak_process_working_bytes=peak))
        print(json.dumps(dict(label=label,checks=checks,seconds=runs[-1]['seconds'])),flush=True)
    repeated=all(runs[0]['outputs'][name]['sha256']==runs[2]['outputs'][name]['sha256'] for name in runs[0]['outputs'])
    result=dict(protocol='docs/verification/avl-reference-execution-protocol.md',
        source_sha256=digest(Path(__file__)),input_manifest_sha256=digest(manifest),
        executable=dict(url='https://web.mit.edu/drela/Public/web/avl/avl352.exe',sha256=EXPECTED,bytes=EXE.stat().st_size),
        runs=runs,repeat_outputs_equal=repeated,
        initial_result_sha256=digest(ROOT/'docs/verification/avl-reference-execution-result.json'),
        semantic_checks_pass=all(all(r['semantic']['checks'].values()) and r['peak_process_working_bytes']<512*1024**2 for r in runs),
        transport_checks_pass=all(all(r['checks'].values()) for r in runs) and repeated,
        limitations=['output acquisition and option checks only; no finite-wing physical acceptance',
            'no core-zero source build','no quantitative residual or conservative element force reconstruction yet'])
    output=ROOT/'docs/verification/avl-reference-execution-corrected-result.json'
    if output.exists(): raise FileExistsError('Preserve earlier reference results; choose a new result name')
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')

if __name__=='__main__':main()
