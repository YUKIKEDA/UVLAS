from pathlib import Path
import sys
order=int(sys.argv.pop(2));p=Path(__file__).with_name('run_multiple_surface_harmonic.py');source=p.read_text(encoding='utf-8');source=source.replace("build/analysis/multiple-surface-harmonic'",f"build/analysis/multiple-surface-harmonic-ideal-q{order}'");source=source.replace('1e-5,128.','0.,128.').replace('128.,16)','128.,'+str(order)+')');exec(compile(source,str(p),'exec'))
