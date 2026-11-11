from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parents[2]
side=sys.argv[1];grids=[(128,192)] if side=='candidate' else [(8,96),(16,48),(16,192),(16,384),(32,192),(64,192)]
for m,n in grids:
    print('START',side,m,n,flush=True)
    subprocess.run([sys.executable,str(ROOT/'tools/design-experiments/run_frequency_delta.py'),side,str(m),str(n)],check=True,timeout=1200,cwd=ROOT)
