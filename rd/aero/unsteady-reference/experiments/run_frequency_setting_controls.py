import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import sys,time,json,numpy as np
from uvlm_frequency_delta import Solver
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/frequency-settings'
side=sys.argv[1];OUT.mkdir(exist_ok=True);start=time.monotonic()
if side=='candidate':
    from uvlm_steady_reference import geometry
else:
    from ptera_frequency_reference import geometry
g=geometry(4,24)
for eps,length,order,label in [(1e-5,128.,16,'base'),(1e-5,64.,16,'tail64'),(1e-5,256.,16,'tail256'),(1e-4,128.,16,'core1e-4'),(1e-6,128.,16,'core1e-6'),(1e-5,128.,24,'q24')]:
    solver=Solver(g,side,eps,length)
    for k in [.05,.2,.8]:
        dest=OUT/f'{side}-{label}-k{k}.npz';assert not dest.exists();r=solver.run(k,order);np.savez_compressed(dest,transfer=r['transfer'],gamma=r['gamma']);print(json.dumps(dict(side=side,setting=label,k=k,residuals=r['residuals'],seconds=time.monotonic()-start)),flush=True);assert time.monotonic()-start<180
