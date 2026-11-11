import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
p=rd_paths.sibling('run_ptera_ideal_frequency.py');source=p.read_text(encoding='utf-8').replace("build/analysis/frequency-ideal-reference'","build/analysis/frequency-ideal-reference-q24'").replace('solver.run(k)','solver.run(k,order=24)');exec(compile(source,str(p),'exec'))
