import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np,json
from trefftz_observations import midpoint
ROOT=rd_paths.ROOT;OLD=ROOT/'build/analysis/finite-wing-2-3';NEW=ROOT/'build/analysis/drag-common-limit'
def trace(n,d):
    t=(1-np.cos(np.linspace(0,np.pi,n//2+1)))/2;y=np.r_[4*(t-1),4*t[1:]];return np.column_stack([np.zeros(n+1),y,abs(y)*np.tan(np.deg2rad(d))])
def circulation(raw,m,n):
    h=raw['horse'];mu=h.reshape(n,m).sum(axis=1)
    if 'left' in raw.files:
        y=raw['left'][:,1].reshape(n,m);assert max(np.ptp(y,axis=1))<1e-12;assert np.all(np.diff(y[:,0])>0)
    return mu
def main():
    records=[]
    for d in [0,5]:
        for m,n,key in [(16,96,'m16-n96'),(32,96,'m32-n96'),(16,192,'m16-n192'),(32,192,'m32-n192'),(16,384,'span384'),(8,192,'band-n192'),(8,384,'band-n384'),(8,768,'band-n768')]:
            p=OLD/f'reference-d{d}-{key}.npz';raw=np.load(p);mu=circulation(raw,m,n);D=midpoint(trace(n,d),mu);records.append(dict(dihedral=d,m=m,n=n,core=1e-8,far_drag=float(D),source=p.relative_to(ROOT).as_posix()));print(records[-1],flush=True)
        for core in [1e-9,1e-10]:
            p=NEW/f'reference-d{d}-m8-n768-core{core}.npz';raw=np.load(p);mu=circulation(raw,8,768);D=midpoint(trace(768,d),mu);records.append(dict(dihedral=d,m=8,n=768,core=core,far_drag=float(D),source=p.relative_to(ROOT).as_posix()));print(records[-1],flush=True)
    (NEW/'far-reference-series.json').write_text(json.dumps(dict(records=records),indent=2),encoding='utf-8')
if __name__=='__main__':main()
