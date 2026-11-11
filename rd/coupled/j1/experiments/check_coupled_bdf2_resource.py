"""Repeat unchanged numerical checks to complete the missing peak-memory measure."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from run_frequency_delta import peak_bytes
p=rd_paths.sibling('check_coupled_bdf2_prescribed.py');source=p.read_text(encoding='utf-8').replace("'rd/coupled/j1/records/coupled-bdf2-prescribed-result.json'","'rd/coupled/j1/records/coupled-bdf2-resource-result.json'").replace('paths=[Path(__file__),','paths=[Path(__file__),p,ROOT/\'rd/coupled/j1/records/coupled-bdf2-prescribed-result.json\',').replace('result=dict(passed=passed,','measured_peak=peak_bytes();result=dict(passed=passed and measured_peak<2*2**30,peak_bytes=measured_peak,');exec(compile(source,str(p),'exec'))
