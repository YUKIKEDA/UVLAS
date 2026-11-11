import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib
import json
import numpy as np
import check_avl_linear_load_wrench as wrench

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    parent_path=FOLDER/'avl-band-aligned-grid-runs.json';parent=json.loads(parent_path.read_text());assert parent['completed']
    cases=[]
    for run in parent['runs']:
        folder=ROOT/run['directory'];rows=np.loadtxt(folder/'uv-rows.txt',skiprows=2);forces=np.loadtxt(folder/'uv-forces.txt',skiprows=2)
        a=rows[:,7:10];b=rows[:,10:13];finite=forces[forces[:,2]==0];other=forces[forces[:,2]!=0]
        p=finite[:,3:6];f=finite[:,6:9]
        couple=wrench.convert(a,b,p,f);r=wrench.observe(a,b,f,couple)
        trailing=wrench.old.aggregate(other[:,3:6],other[:,6:9]);point=wrench.old.aggregate(forces[:,3:6],forces[:,6:9])
        force_difference=float(np.max(abs(r['force']+trailing['force']-point['force']))/wrench.S)
        moment_difference=float(np.max(abs((r['moment']+trailing['moment']-point['moment'])/(wrench.S*np.array([wrench.B,wrench.C,wrench.B])))))
        low=np.minimum(a[:,1],b[:,1]);high=np.maximum(a[:,1],b[:,1]);bounds=wrench.BOUNDS[1:-1]
        crossing=np.any((low[:,None]<bounds[None,:]-1e-9)&(high[:,None]>bounds[None,:]+1e-9),axis=1)
        cases.append(dict(case=run['case'],m=run['m'],n=run['n'],alpha=run['alpha'],crossing_source_count=int(np.count_nonzero(crossing)),
            point_vs_reconstruction_force=force_difference,point_vs_reconstruction_moment=moment_difference))
    checks=dict(supports_do_not_cross=all(r['crossing_source_count']==0 for r in cases),
        point_reconstruction_equivalence=all(max(r['point_vs_reconstruction_force'],r['point_vs_reconstruction_moment'])<=1e-10 for r in cases))
    result=dict(protocol='rd/aero/avl-reference/records/avl-band-aligned-grid-experiment.md',source_sha256=digest(Path(__file__)),
        dependency_sha256=digest(Path(wrench.__file__)),parent_result_sha256=digest(parent_path),cases=cases,checks=checks,all_criteria_pass=all(checks.values()))
    (FOLDER/'avl-band-aligned-observation-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(checks=checks,max_force=max(r['point_vs_reconstruction_force'] for r in cases),max_moment=max(r['point_vs_reconstruction_moment'] for r in cases))))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
