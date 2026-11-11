from pathlib import Path
import json,hashlib,time,numpy as np
from ptera_frequency_reference import geometry,project
from ptera_scalar_reference import baseline,delta,line
import uvlm_frequency_delta as old
ROOT=Path(__file__).resolve().parents[2]
def main():
    start=time.monotonic();records=[]
    for d in [0,5]:
        g=geometry(4,24)
        if d:
            for key in ['X','P','left','right']:g[key][...,2]=abs(g[key][...,1])*np.tan(np.deg2rad(d))
            g['normal']=np.column_stack([np.zeros(len(g['P'])),-np.sign(g['P'][:,1])*np.sin(np.deg2rad(d)),np.full(len(g['P']),np.cos(np.deg2rad(d)))])
        for eps in [0.,1e-5,.03]:
            expected=old.baseline(g,'reference',eps,128.);actual=baseline(g,'reference',eps,128.);matrix_error=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected));assert matrix_error<1e-9
            for k in [.05,.2,.8]:
                w=40*k/1.114285714285714;expected=old.delta(g,'reference',w,eps,128.,16);actual=delta(g,'reference',w,eps,128.,16);error=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected));assert error<1e-9;records.append(dict(dihedral=d,eps=eps,k=k,matrix_error=matrix_error,wake_error=error))
    native=ROOT/'build/reference-ptera/env/Lib/site-packages/pterasoftware/_aerodynamics_functions.py';paths=[Path(__file__),ROOT/'tools/design-experiments/ptera_scalar_reference.py',native,ROOT/'docs/verification/ptera-scalar-reference-protocol.md'];result=dict(passed=True,records=records,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'docs/verification/ptera-scalar-reference-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
