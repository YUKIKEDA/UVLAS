from pathlib import Path
p=Path(__file__).with_name('run_ptera_ideal_frequency.py');source=p.read_text(encoding='utf-8').replace("build/analysis/frequency-ideal-reference'","build/analysis/frequency-ideal-reference-q24'").replace('solver.run(k)','solver.run(k,order=24)');exec(compile(source,str(p),'exec'))
