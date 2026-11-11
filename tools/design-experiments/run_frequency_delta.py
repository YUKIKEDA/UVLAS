from pathlib import Path
import sys,json,time,ctypes,numpy as np
from uvlm_frequency_delta import Solver
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/frequency-delta'
def peak_bytes():
    class Counters(ctypes.Structure):
        _fields_=[('cb',ctypes.c_ulong),('faults',ctypes.c_ulong)]+[(name,ctypes.c_size_t) for name in ['peak_working','working','peak_paged','paged','peak_nonpaged','nonpaged','pagefile','peak_pagefile']]
    v=Counters();v.cb=ctypes.sizeof(v);h=ctypes.windll.kernel32.GetCurrentProcess;h.restype=ctypes.c_void_p;query=ctypes.windll.psapi.GetProcessMemoryInfo;query.argtypes=[ctypes.c_void_p,ctypes.POINTER(Counters),ctypes.c_ulong]
    if not query(h(),ctypes.byref(v),v.cb):raise RuntimeError('memory_query_failed')
    return int(v.peak_working)
def main():
    side=sys.argv[1];m,n=map(int,sys.argv[2:4]);OUT.mkdir(exist_ok=True);start=time.monotonic()
    if side=='candidate':
        from uvlm_steady_reference import geometry
        g=geometry(m,n)
    else:
        from ptera_frequency_reference import geometry
        g=geometry(m,n)
    solver=Solver(g,side)
    for k in [.05,.2,.8]:
        p=OUT/f'{side}-m{m}-n{n}-k{k}.json'
        if p.exists():continue
        r=solver.run(k);peak=peak_bytes();assert peak<24*2**30;np.savez_compressed(p.with_suffix('.npz'),gamma=r['gamma'],transfer=r['transfer'],W=r['W']);meta=dict(side=side,m=m,n=n,k=k,residuals=r['residuals'],refinement_steps=r['refinement_steps'],seconds=time.monotonic()-start,peak_bytes=peak,transfer=[[[float(z.real),float(z.imag)] for z in v] for v in r['transfer']]);p.write_text(json.dumps(meta,indent=2),encoding='utf-8');print(json.dumps(meta),flush=True)
if __name__=='__main__':main()
