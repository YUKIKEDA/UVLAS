"""Run the frozen native-core comparison with epsilon=0, leaving it untouched."""
from pathlib import Path
import sys
p=Path(__file__).with_name('run_multiple_surface_harmonic.py');source=p.read_text(encoding='utf-8');source=source.replace("build/analysis/multiple-surface-harmonic'","build/analysis/multiple-surface-harmonic-ideal'");source=source.replace('1e-5,128.','0.,128.');exec(compile(source,str(p),'exec'))
