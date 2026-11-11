"""Repeat unchanged numerical checks to complete the missing peak-memory measure."""
from pathlib import Path
from run_frequency_delta import peak_bytes
p=Path(__file__).with_name('check_coupled_bdf2_prescribed.py');source=p.read_text(encoding='utf-8').replace("'docs/verification/coupled-bdf2-prescribed-result.json'","'docs/verification/coupled-bdf2-resource-result.json'").replace('paths=[Path(__file__),','paths=[Path(__file__),p,ROOT/\'docs/verification/coupled-bdf2-prescribed-result.json\',').replace('result=dict(passed=passed,','measured_peak=peak_bytes();result=dict(passed=passed and measured_peak<2*2**30,peak_bytes=measured_peak,');exec(compile(source,str(p),'exec'))
