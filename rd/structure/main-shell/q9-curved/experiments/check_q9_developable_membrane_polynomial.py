"""Independent Q2 cardinal derivatives versus the analytic P2 membrane defect."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json
import numpy as np
import q9_developable_bending_approximation as method
core=method.core;F=core.F;ROOT=core.ROOT
def main():
    target=F/'q9-developable-membrane-polynomial-result.json'
    if target.exists():raise FileExistsError(target)
    records=[];k=.4
    for n in [1,2,4,8]:
        maximum=0.;sample_zero=0.
        for cell in range(n):
            center=(cell+.5)/n;half=1/(2*n);nodes=center+half*np.array([-1.,0.,1.]);U=-k*nodes**3/3;W=nodes**2/2
            for s in [-.91,-1/np.sqrt(3),0.,1/np.sqrt(3),.83]:
                derivative=np.array([s-.5,-2*s,s+.5])/half;x=center+half*s;actual=derivative@U+k*x*(derivative@W);expected=k*half**2*(s*s-1/3);maximum=max(maximum,float(abs(actual-expected)))
                if abs(abs(s)-1/np.sqrt(3))<1e-15:sample_zero=max(sample_zero,float(abs(actual)))
        records.append(dict(n=n,polynomial_difference=maximum,two_gauss_point_value=sample_zero,checks_pass=bool(maximum<=1e-12 and sample_zero<=1e-12)))
    result=dict(source_sha256=core.digest(__file__),contract_sha256=core.digest(ROOT/'rd/structure/main-shell/q9-curved/q9-curved-membrane-compatibility-contract.md'),cases=records,all_criteria_pass=bool(all(r['checks_pass'] for r in records)),new_reduction_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
