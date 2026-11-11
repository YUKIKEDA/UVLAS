from pathlib import Path
import numpy as np,json
from probe_ptera_reference import run
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/unsteady-2-3'
def main():
    s,_,base=run();core=float(s._currentStackBoundRc0s[0]);audit,_,r=run(core=core);old=np.load(OUT/'original.npz');assert np.array_equal(old['F'],base['F']) and np.array_equal(old['M'],base['M']);same=np.array_equal(base['F'],r['F']) and np.array_equal(base['M'],r['M']);assert same
    label='original-core-exact-parity';np.savez_compressed(OUT/(label+'.npz'),**r);(OUT/(label+'.json')).write_text(json.dumps(dict(label=label,actual_original_core=core,rounded_initial_core=.0315,bitwise_loads=True,seconds=r['seconds']),indent=2),encoding='utf-8');print('BITWISE_PARITY',repr(core),same)
if __name__=='__main__':main()
