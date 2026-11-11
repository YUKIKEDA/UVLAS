from pathlib import Path
import sys,subprocess
ROOT=Path(__file__).resolve().parents[2]
for m,n in [(8,192),(16,48),(16,96),(16,192),(16,384),(32,96),(32,192),(64,192)]:
    print('START',m,n,flush=True);subprocess.run([sys.executable,str(ROOT/'tools/design-experiments/run_ptera_ideal_frequency.py'),str(m),str(n)],check=True,timeout=1200,cwd=ROOT)
