import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import sys
order=int(sys.argv.pop(2));p=rd_paths.sibling('run_multiple_surface_harmonic.py');source=p.read_text(encoding='utf-8');source=source.replace("build/analysis/multiple-surface-harmonic'",f"build/analysis/multiple-surface-harmonic-ideal-q{order}'");source=source.replace('1e-5,128.','0.,128.').replace('128.,16)','128.,'+str(order)+')');exec(compile(source,str(p),'exec'))
