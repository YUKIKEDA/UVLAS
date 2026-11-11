import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib
import json
import numpy as np
import check_avl_collinear_kernel as kernels

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    path=FOLDER/'avl-stable-reference-runs.json';runs=json.loads(path.read_text())
    previous_path=FOLDER/'avl-source-measurement-corrected-runs.json';previous=json.loads(previous_path.read_text())
    analysis_path=FOLDER/'avl-stable-reference-result.json';analysis=json.loads(analysis_path.read_text())
    comparisons=[];pair_comparisons=[];unchanged=[]
    for run in runs['runs']:
        if run['variant']=='baseline':continue
        old=next(r for r in previous['runs'] if all(r[k]==run[k] for k in ['variant','case','m','n','alpha']))
        folder=ROOT/run['directory'];old_folder=ROOT/old['directory']
        if run['variant']=='measured':
            unchanged.append(dict(case=run['case'],m=run['m'],n=run['n'],alpha=run['alpha'],
                outputs_equal=all(digest(folder/n)==digest(old_folder/n) for n in ['total.txt','strips.txt','elements.txt'])))
        else:
            rows=np.loadtxt(folder/'uv-rows.txt',skiprows=2);old_rows=np.loadtxt(old_folder/'uv-rows.txt',skiprows=2)
            m,n=run['m'],run['n'];mirror=np.array([(n-1-j)*m+k for j in range(n) for k in range(m)])
            gamma=rows[:,13]
            gamma_change=float(np.max(abs(gamma-old_rows[:,13])/np.maximum(abs(old_rows[:,13]),1.)))
            gamma_mirror=float(np.max(abs(gamma-gamma[mirror])/np.maximum(abs(gamma),1.)))
            coefs=next(c['coefficients'] for c in analysis['cases'] if all(c[k]==run[k] for k in ['variant','case','m','n','alpha']))
            comparisons.append(dict(case=run['case'],m=m,n=n,alpha=run['alpha'],gamma_change=gamma_change,gamma_mirror=gamma_mirror,
                side_coefficient=coefs['force'][1],yaw_coefficient=coefs['moment'][2],
                eligible=gamma_change<=1e-10 and gamma_mirror<=1e-9 and abs(coefs['force'][1])<=1e-9 and abs(coefs['moment'][2])<=1e-9))
            if run['alpha']==1:
                data=np.loadtxt(folder/'uv-force-points.txt',skiprows=2);p=data[:,1:4]
                a=rows[:,7:10];b=rows[:,10:13];size=m*n
                bound,_=kernels.stable_bound(p,a,b,np.eye(size,dtype=bool))
                native_trail=kernels.native_parts(p,a,b,np.eye(size,dtype=bool))[1]
                expected=bound+native_trail
                actual=np.loadtxt(folder/'uv-force-matrix.txt',skiprows=1)[:,2:].reshape(size,size,3)
                diff=np.linalg.norm(actual-expected,axis=-1)
                target=1e-10+1e-10*np.linalg.norm(expected,axis=-1)
                pair_comparisons.append(dict(case=run['case'],m=m,n=n,pairs=size**2,
                    maximum_absolute_difference=float(np.max(diff)),eligible=bool(np.all(diff<=target))))
    checks=dict(base_fourteen=analysis['all_criteria_pass'],nonzero_core_unchanged=all(r['outputs_equal'] for r in unchanged),
        zero_core_symmetry_and_gamma=all(r['eligible'] for r in comparisons),
        independent_pair_kernel=all(r['eligible'] for r in pair_comparisons))
    result=dict(protocol='rd/aero/avl-reference/records/avl-collinear-kernel-experiment.md#6',source_sha256=digest(Path(__file__)),
        dependency_sha256=digest(Path(kernels.__file__)),parent_runs_sha256=digest(path),
        analysis_sha256=digest(analysis_path),previous_runs_sha256=digest(previous_path),
        nonzero_core_comparisons=unchanged,zero_core_comparisons=comparisons,pair_comparisons=pair_comparisons,
        checks=checks,all_criteria_pass=all(checks.values()),
        status='limited corrected reference candidate; full grid, conservative strips and UVLAS comparison still pending')
    (FOLDER/'avl-stable-reference-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(checks=checks,max_side=max(abs(r['side_coefficient']) for r in comparisons),
        max_yaw=max(abs(r['yaw_coefficient']) for r in comparisons),max_gamma_change=max(r['gamma_change'] for r in comparisons))))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
