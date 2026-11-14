import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,numpy as np
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/finite-wing-2-3'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    oldpath=ROOT/'rd/aero/avl-reference/records/fixed-band-final-extension-result.json';old=json.loads(oldpath.read_text());assert not old['passed'];records=[];artifacts={}
    def get(d,key):
        p=OUT/f'reference-d{d}-{key}.json';artifacts[p.relative_to(ROOT).as_posix()]=sha(p);return json.loads(p.read_text())
    for d in [0,5]:
        c8=get(d,'band-n768');c9=get(d,'band768-core1e-09');c10=get(d,'band768-core1e-10');last=abs(np.array(c10['bands'])/c9['bands']-1);prev=abs(np.array(c9['bands'])/c8['bands']-1);assert max(last)<max(prev) and max(last)<.0005
        p=OUT/f'candidate-d{d}-band-n768.json';artifacts[p.relative_to(ROOT).as_posix()]=sha(p);own=json.loads(p.read_text());diff=abs(np.array(own['bands'])/c10['bands']-1);rows=[r for r in old['records'] if r['case']=='fixed_band_space' and r['dihedral']==d];space=sum(np.array(r['sum_indicator']) for r in rows);correction=abs(np.array(c10['bands'])/c8['bands']-1);combined=space+correction+diff+last;assert max(combined)<.02
        records.append(dict(dihedral=d,core_1e8_to_1e9=float(max(prev)),core_1e9_to_1e10=float(max(last)),selected_reference_core=1e-10,max_independent_band_difference=float(max(diff)),max_combined_band_indicator=float(max(combined)),combined_band_indicators=combined.tolist(),passed=True,rigorous_error_bound=False))
    sources={p.relative_to(ROOT).as_posix():sha(p) for p in [Path(__file__),oldpath,ROOT/'rd/aero/avl-reference/experiments/probe_fixed_band_core_plateau.py',ROOT/'rd/aero/avl-reference/records/fixed-band-core-plateau-protocol.md']};r=dict(passed=True,records=records,sources=sources,artifacts=artifacts,scope='FW-V01 and FW-V02 first-order conservative 20-band loads',uvlm_2_3_complete=False);dest=ROOT/'rd/aero/avl-reference/records/fixed-band-core-plateau-result.json';assert not dest.exists();dest.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,records=records),indent=2))
if __name__=='__main__':main()
