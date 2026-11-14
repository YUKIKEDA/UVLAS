import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib
import json
import math
import time
import numpy as np
import check_q9_plate_convergence as memory

ROOT=rd_paths.ROOT;START=time.monotonic()
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def primitive(y):
    y=np.asarray(y);return .5*(y*np.sqrt(np.maximum(0.,1-(y/4)**2))+4*np.arcsin(np.clip(y/4,-1.,1.)))

def main():
    parent_path=ROOT/'rd/aero/avl-reference/records/avl-fixed-band-convergence-result.json';parent=json.loads(parent_path.read_text())
    runs_path=ROOT/'rd/aero/avl-reference/records/avl-fixed-band-grid-runs.json';runs=json.loads(runs_path.read_text())
    bounds=np.array(parent['bounds']);exact=np.diff(primitive(bounds));records=[]
    for ident in ['FW-V01','FW-V02']:
        for n in [12,24,48,96]:
            case=next(c for c in parent['cases'] if c['case']==ident and c['m']==32 and c['n']==n and c['alpha']==1)
            run=next(c for c in runs['runs'] if c['case']==ident and c['m']==32 and c['n']==n and c['alpha']==1)
            forces=np.loadtxt(ROOT/run['directory']/'uv-forces.txt',skiprows=2)
            points=forces[forces[:,2]==0,3:6].reshape(n,32,3)[:,0,:]
            intervals=np.array(case['input_strip_intervals'])
            source=primitive(intervals[:,1])-primitive(intervals[:,0])
            owner=np.searchsorted(bounds,points[:,1],side='right')-1
            concentrated=np.array([math.fsum(source[owner==j]) for j in range(20)])
            overlap=np.maximum(0.,np.minimum(bounds[1:,None],intervals[None,:,1])-np.maximum(bounds[:-1,None],intervals[None,:,0]))
            reconstructed=overlap@(source/(intervals[:,1]-intervals[:,0]))
            rel=lambda v:float(np.linalg.norm(v-exact)/np.linalg.norm(exact))
            records.append(dict(case=ident,n=n,source_loads=source.tolist(),source_points_y=points[:,1].tolist(),
                concentrated=concentrated.tolist(),reconstructed=reconstructed.tolist(),
                concentrated_relative=rel(concentrated),reconstructed_relative=rel(reconstructed),
                concentrated_conservation=abs(math.fsum(concentrated)-2*math.pi),
                reconstructed_conservation=abs(math.fsum(reconstructed)-2*math.pi)))
    checks=dict(conservative=all(max(c['concentrated_conservation'],c['reconstructed_conservation'])<=1e-12 for c in records),
        concentration_error_detected=any(c['concentrated_relative']>.01 for c in records),
        reconstruction_fine_accuracy=all(c['reconstructed_relative']<=.01 for c in records if c['n']==96))
    result=dict(protocol='rd/aero/avl-reference/records/avl-fixed-band-convergence-experiment.md#7',source_sha256=digest(Path(__file__)),
        dependency_sha256=digest(Path(memory.__file__)),parent_result_sha256=digest(parent_path),runs_sha256=digest(runs_path),
        exact_bands=exact.tolist(),cases=records,checks=checks,all_criteria_pass=all(checks.values()),
        seconds=time.monotonic()-START,peak_process_working_bytes=memory.memory_bytes(),
        limits='manufactured distribution only; reconstruction is not substituted into original AVL point-load observations')
    (ROOT/'rd/aero/avl-reference/records/fixed-band-manufactured-distribution-result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(checks=checks,fine=[dict(case=c['case'],point=c['concentrated_relative'],reconstruction=c['reconstructed_relative']) for c in records if c['n']==96])))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
