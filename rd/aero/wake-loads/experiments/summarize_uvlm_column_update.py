"""Derived F/M and fixed-window TE coefficient diagnostics; no new solve."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import hashlib
import json
from pathlib import Path
import numpy as np
import check_uvlm_local_pressure_pitch as reference

ROOT=rd_paths.ROOT


def main():
    path=ROOT/'rd/aero/wake-loads/records/uvlm-column-update-result.json'
    source=Path(__file__);result=json.loads(path.read_text());rows=[]
    for row in result['continuous_cases']:
        mode=next(v for v in reference.MODES if v['name']==row['mode']);k=row['k']
        b,_=reference.coefficients(k,mode);limit=(2/3)*1.2*8**2*(b[0]+4*b[1]+8*b[2])
        ps,fs,ms=reference.scales(mode);loads=np.array([complex(*v) for v in row['loads']])
        exact=reference.reference_load(k,mode)
        differences=(np.abs(loads-exact)/np.maximum(np.abs(exact),.01*np.array([fs,ms]))).tolist()
        windows=[]
        # The last four entries are the prescribed cumulative TE windows;
        # the ordinary bins separately contain the same [.95, 1] interval.
        for p in row['pressure_profile'][-4:]:
            width=p['upper']-p['lower']
            if p['upper']==1. and any(abs(width-w)<1e-12 for w in [.1,.05,.025,.0125]):
                value=complex(*p['pressure'])/np.sqrt(width)
                finite=complex(*p['reference'])/np.sqrt(width)
                windows.append(dict(width=width,divided_by_sqrt=reference.P.pair(value),
                    analytic_finite_window=reference.P.pair(finite),analytic_limit=reference.P.pair(limit),
                    normalized_limit_difference=float(abs(value-limit)/max(ps,abs(limit))),
                    normalized_finite_window_difference=float(abs(value-finite)/max(ps,abs(finite)))))
        assert len(windows)==4
        rows.append(dict(m=row['m'],k=k,mode=row['mode'],loads=row['loads'],
            reference_loads=[reference.P.pair(v) for v in exact],load_complex_errors=differences,trailing_windows=windows))
    data=dict(parent_result_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={'check_uvlm_local_pressure_pitch.py':hashlib.sha256(rd_paths.sibling('check_uvlm_local_pressure_pitch.py').read_bytes()).hexdigest()},
        cases=rows,limits='derived diagnostics from fixed windows; not a new solve or P0 pointwise pressure acceptance')
    (ROOT/'rd/aero/wake-loads/records/uvlm-column-update-observation-audit.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(cases=len(rows),fixed_TE_windows=sum(len(r['trailing_windows']) for r in rows),new_solves=0)))


if __name__=='__main__':main()
