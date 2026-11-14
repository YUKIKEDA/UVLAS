"""Algebraic post-hoc decomposition of allocation-dependent moment observations."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,math,time
import numpy as np
import avl_reference_resources as provenance
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;digest=provenance.digest
def main():
    output=F/'avl-tip-couple-component-result.json'
    if output.exists():raise FileExistsError(output)
    start=time.monotonic();parent=F/'avl-tip-redistribution-result.json';r=json.loads(parent.read_text());assert r['methods_criteria_pass']
    records={(m['case'],m['allocation']):m for m in r['meshes']};pairs=[];maximum=0.
    for p in r['pairs']:
        a=records[(p['case'],p['from_allocation'])]['bands'];b=records[(p['case'],p['to_allocation'])]['bands']
        h=np.array([v['width'] for v in a]);couple=-np.array([y['Cy']-x['Cy'] for x,y in zip(a,b)])/(2*np.sqrt(h))
        line=np.array(p['channels']['j_line']['difference']);point=np.array(p['channels']['j_point']['difference'])
        error=float(np.max(abs(line+couple-point)));norm2=float(point@point);line2=float(line@line);couple2=float(couple@couple);cross=float(2*line@couple)
        maximum=max(maximum,error,abs(norm2-line2-couple2-cross))
        scalar=p['scalar_changes']['moment_point']['signed']['total'];line_scalar=p['scalar_changes']['moment_line']['signed']['total']
        # Existing S=8.4, Cref=1.1142857142857143.
        scalar_couple=math.fsum(y['Cy']-x['Cy'] for x,y in zip(a,b))/(8.4*1.1142857142857143)
        maximum=max(maximum,abs(scalar-line_scalar-scalar_couple))
        pairs.append(dict(case=p['case'],from_allocation=p['from_allocation'],to_allocation=p['to_allocation'],line_difference=line.tolist(),couple_difference=couple.tolist(),point_difference=point.tolist(),
            reconstruction_error=error,line_norm=math.sqrt(line2),couple_norm=math.sqrt(couple2),point_norm=math.sqrt(norm2),
            line_to_point_norm=math.sqrt(line2/norm2) if norm2 else None,couple_to_point_norm=math.sqrt(couple2/norm2) if norm2 else None,
            squared_terms=dict(line=line2,couple=couple2,cross=cross,total=norm2),
            scalar_line=line_scalar,scalar_couple=scalar_couple,scalar_point=scalar))
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),protocol_sha256=digest(F/'avl-tip-couple-diagnostic-addendum.md'),
        dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [provenance,memory]},pairs=pairs,maximum_reconstruction_error=maximum,
        seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),new_avl_runs=0,reference_acceptance=False,
        all_criteria_pass=len(pairs)==6 and maximum<=1e-12 and time.monotonic()-start<60 and memory.memory_bytes()<128*1024**2)
    output.write_text(json.dumps(result,indent=2)+'\n');assert output.stat().st_size<1024**2 and result['all_criteria_pass']
    print(json.dumps(dict(maximum_reconstruction_error=maximum,ratios=[{k:p[k] for k in ['case','from_allocation','to_allocation','line_to_point_norm','couple_to_point_norm']} for p in pairs]),indent=2))
if __name__=='__main__':main()
