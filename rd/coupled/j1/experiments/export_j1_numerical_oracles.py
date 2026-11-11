import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,numpy as np
from scipy.linalg import solve,eigvals
from pg_uvlm_free_step import FreeProblem
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/numeric-oracles'
def main():
    OUT.mkdir(exist_ok=True);p=FreeProblem();s=p.sstar;x=p.guess(s);_,_,J=p.evaluate(x,s,True);d=p.dimension;Rp=np.empty((16,d));Cp=np.empty((d,d));Cx=np.empty((d,16));h=2e-6
    for j in range(d):
        v=np.eye(d)[j];a,ca,_=p.evaluate(x,s+h*v);b,cb,_=p.evaluate(x,s-h*v);Rp[:,j]=(a-b)/(2*h);Cp[:,j]=(ca-cb)/(2*h)
    for j in range(16):
        v=np.eye(16)[j];Cx[:,j]=(p.evaluate(x+h*v,s)[1]-p.evaluate(x-h*v,s)[1])/(2*h)
    P=Cp-Cx@solve(J,Rp);path=OUT/'pg-free-step.txt';assert not path.exists()
    with path.open('w') as f:f.write(str(d)+'\n');np.savetxt(f,P,fmt='%.17g')
    values=eigvals(P);np.savez_compressed(OUT/'pg-free-oracle.npz',P=P,eigenvalues=values,M=p.mass,K=p.internal(p.star,True)[1]);meta=dict(source='actual deformed PG/free UVLM full history step',dimension=d,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),spectral_radius=float(max(abs(values))));(OUT/'oracle.json').write_text(json.dumps(meta,indent=2),encoding='utf-8');print(json.dumps(meta))
if __name__=='__main__':main()
