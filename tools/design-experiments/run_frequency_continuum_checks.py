from pathlib import Path
import argparse,json,time,gc,numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/unsteady-2-3'
def main():
    p=argparse.ArgumentParser();p.add_argument('side');args=p.parse_args()
    if args.side=='candidate':from uvlm_frequency_fast import run
    else:from ptera_frequency_reference import run
    tasks=[(8,96,12,False),(16,192,12,False)]
    if args.side=='reference':tasks.append((64,192,16,True))
    for m,n,order,refine in tasks:
        label=f'{args.side}-continuum-m{m}-n{n}-q{order}';p=OUT/(label+'.json')
        if p.exists():continue
        print('START',label,flush=True);start=time.monotonic();kw=dict(length=128.,order=order)
        if args.side=='reference':kw['refine_near']=refine
        r=run(m,n,**kw);coef=[r['F'][2],r['M'][1]];result=dict(label=label,m=m,n=n,length=128.,order=order,refine_near=refine,eps=1e-5,transfer=[[float(z.real),float(z.imag)] for z in coef],residual=r['residual'],seconds=time.monotonic()-start);np.savez_compressed(p.with_suffix('.npz'),gamma=r['gamma'],F=r['F'],M=r['M'],W=r['W']);p.write_text(json.dumps(result,indent=2),encoding='utf-8');print('DONE',label,json.dumps(result),flush=True);del r;gc.collect()
if __name__=='__main__':main()
