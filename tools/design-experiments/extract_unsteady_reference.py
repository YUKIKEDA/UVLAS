"""Correct pole shift from retained raw moments; do not overwrite exploratory records."""
from pathlib import Path
import json,numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/unsteady-2-3'
def extract(label):
    meta=json.loads((OUT/(label+'.json')).read_text());raw=np.load(OUT/(label+'.npz'));t=raw['t'];omega=float(raw['omega']);amp=float(raw['amplitude']);last=t>=t[-1]-2*np.pi/omega-1e-9;L=-raw['F'][:,2];My=raw['M'][:,1]+.375*L
    D=np.column_stack([np.ones(last.sum()),np.cos(omega*t[last]),np.sin(omega*t[last])]);fit=np.linalg.lstsq(D,np.column_stack([L[last],My[last]]),rcond=None)[0];coef=(fit[1]-1j*fit[2])/amp;res=np.linalg.norm(D@fit-np.column_stack([L[last],My[last]]),axis=0)/np.linalg.norm(np.column_stack([L[last],My[last]]),axis=0)
    return meta,coef,res
def main():
    for label in ['original','original-core-parity','time64','time128','time256','cycles6','space4x24','space8x48']:
        if (OUT/(label+'.json')).exists():
            m,c,r=extract(label);print(label,c.tolist(),'fit residual',r.tolist())
if __name__=='__main__':main()
