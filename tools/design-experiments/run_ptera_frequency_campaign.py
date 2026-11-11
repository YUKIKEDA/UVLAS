from pathlib import Path
import json,time,gc,numpy as np
from ptera_frequency_reference import run
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/unsteady-2-3'
def main():
    tasks=[(f'ptera-frequency-m{m}-n{n}',m,n,256.,12,1e-5) for m,n in [(2,8),(4,24),(8,48),(16,48),(16,96)]]
    tasks += [('ptera-frequency-quadrature8',8,48,256.,8,1e-5),('ptera-frequency-length128',8,48,128.,12,1e-5),('ptera-frequency-length512',8,48,512.,12,1e-5),('ptera-frequency-core1e-4',8,48,256.,12,1e-4)]
    for label,m,n,L,order,eps in tasks:
        p=OUT/(label+'.json')
        if p.exists():continue
        print('START',label,flush=True);start=time.monotonic();r=run(m,n,length=L,order=order,eps=eps);coef=[r['F'][2],r['M'][1]];result=dict(label=label,m=m,n=n,length=L,order=order,eps=eps,transfer=[[float(z.real),float(z.imag)] for z in coef],residual=r['residual'],seconds=time.monotonic()-start);np.savez_compressed(p.with_suffix('.npz'),gamma=r['gamma'],F=r['F'],M=r['M'],W=r['W']);p.write_text(json.dumps(result,indent=2),encoding='utf-8');print('DONE',label,json.dumps(result),flush=True);del r;gc.collect()
if __name__=='__main__':main()
