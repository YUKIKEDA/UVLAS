"""Additional predeclared nonuniform physical bound panels, continuous wake."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import check_uvlm_trailing_edge as diagnostic

ROOT=rd_paths.ROOT
START=time.monotonic()


def install_grid(m,graded):
    grid=diagnostic.previous.geometry.make_grid(m,1,10000.,wake=False)
    bounds=(1-np.cos(np.pi*np.arange(m+1)/m))/2 if graded else np.linspace(0.,1.,m+1)
    widths=np.diff(bounds);xv=np.r_[bounds[:-1]+.25*widths,1.];xc=bounds[:-1]+.75*widths
    grid['points'][:,0]=np.repeat(xv,2);grid['material'][:,0]=np.repeat(xv,2)
    grid['collocation'][:,0]=xc;grid['coll_material'][:,0]=xc
    a=diagnostic.previous.geometry.aic(grid,0.)
    ends=np.array([[grid['points'][i],grid['points'][j]] for i,j in grid['edges']])
    dl=ends[:,1]-ends[:,0];ids=grid['body_ids']
    fc=1.2*8*(dl[ids,1]@grid['body'][ids])/10000.
    p=grid['points'][grid['triangles']]
    area=.5*np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])
    fu=np.bincount(grid['tri_owner'],weights=1.2*area[:,2],minlength=m)/10000.
    diagnostic.BASE.clear()
    diagnostic.BASE[m]=(a,fc,fu,xc)
    return dict(mesh='cosine-physical-panels-v1' if graded else 'uniform',m=m,
        physical_boundaries=bounds.tolist(),vortex_x=xv.tolist(),collocation_x=xc.tolist(),
        minimum_panel_width=float(widths.min()),maximum_panel_width=float(widths.max()),
        triangle_area_vectors=area.tolist(),positive_cells=bool(np.all(widths>0) and np.all(area[:,2]>0)))


def main():
    original=json.loads((ROOT/'rd/aero/wake-loads/records/uvlm-trailing-edge-result.json').read_text())
    rows=[];meshes=[];baseline=[];sensitivity=[];changes=[]
    for m in [16,32,64,128]:
        meshes.append(install_grid(m,False))
        for k in [.1,.5,1.]:
            kw,_=diagnostic.continuous(m,k,160.)
            row,_,_=diagnostic.solve(m,k,kw)
            saved=next(r for r in original['harmonic_cases'] if r['method']=='continuous' and r['derivative']=='exact' and r['m']==m and r['k']==k)
            baseline.append(dict(m=m,k=k,relative=diagnostic.previous.relative(row,saved)))
        meshes.append(install_grid(m,True))
        for k in [.1,.5,1.]:
            kw,error=diagnostic.continuous(m,k,160.)
            row,_,_=diagnostic.solve(m,k,kw);row['mesh']='cosine-physical-panels-v1'
            row['integration_error_estimate']=error;rows.append(row)
            print(json.dumps(row),flush=True)
            if m==128:
                tighter,estimate=diagnostic.continuous(m,k,160.,tight=True)
                tightrow,_,_=diagnostic.solve(m,k,tighter)
                sensitivity.append(dict(k=k,relative=diagnostic.previous.relative(row,tightrow),tight_error_estimate=estimate))
    for k in [.1,.5,1.]:
        coarse=next(r for r in rows if r['m']==64 and r['k']==k)
        fine=next(r for r in rows if r['m']==128 and r['k']==k)
        changes.append(dict(k=k,relative=diagnostic.previous.relative(coarse,fine)))
    fine=[r for r in rows if r['m']==128]
    checks=dict(positive_cells=all(m['positive_cells'] for m in meshes),
        uniform_baseline=all(r['relative']<=1e-12 for r in baseline),
        nonpenetration=all(r['nonpenetration_relative']<=1e-11 for r in rows),
        complex_accuracy=all(r['complex_relative']<=.02 for r in fine),
        amplitude_accuracy=all(r['amplitude_relative']<=.02 for r in fine),
        phase_accuracy=all(abs(r['phase_degrees'])<=1 for r in fine),
        spatial_change=all(c['relative']<=.005 for c in changes),
        integral_sensitivity=all(c['relative']<=1e-8 for c in sensitivity))
    source=Path(__file__);deps=['check_uvlm_trailing_edge.py',*original['dependency_sha256']]
    result=dict(protocol='rd/aero/wake-loads/records/uvlm-trailing-edge-diagnostic.md#6',geometry_id=diagnostic.previous.geometry.GEOMETRY,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={name:hashlib.sha256(rd_paths.sibling(name).read_bytes()).hexdigest() for name in deps},
        meshes=meshes,harmonic_cases=rows,uniform_baseline=baseline,refinement_changes=changes,
        integral_sensitivity=sensitivity,checks=checks,all_criteria_pass=bool(all(checks.values())),
        numpy=np.__version__,scipy=diagnostic.scipy.__version__,U=8.,rho=1.2,H=.001,c=1.,span=10000.,length=160.,
        seconds=time.monotonic()-START,peak_process_working_bytes=diagnostic.previous.memory.memory_bytes(),
        limits='continuous-wake frequency diagnosis only; no time stepping, Kutta pressure, moment or production mesh adoption')
    (ROOT/'rd/aero/wake-loads/records/uvlm-bound-grading-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=diagnostic.previous.core.json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,changes=changes,seconds=result['seconds'])),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':main()
