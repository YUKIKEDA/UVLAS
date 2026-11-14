import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import subprocess,sys
ROOT=rd_paths.ROOT
side=sys.argv[1];grids=[(128,192)] if side=='candidate' else [(8,96),(16,48),(16,192),(16,384),(32,192),(64,192)]
for m,n in grids:
    print('START',side,m,n,flush=True)
    subprocess.run([sys.executable,str(ROOT/'rd/aero/unsteady-reference/experiments/run_frequency_delta.py'),side,str(m),str(n)],check=True,timeout=1200,cwd=ROOT)
