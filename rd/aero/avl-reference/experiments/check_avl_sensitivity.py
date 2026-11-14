"""Independent audits of the predeclared eight systems and both grid families."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from decimal import Decimal, localcontext
import hashlib
import json
import math
import time
import numpy as np
import scipy.linalg
import check_q9_plate_convergence as memory
import check_avl_fixed_bands as bands
import check_avl_linear_load_wrench as wrench

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(x):return Decimal.from_float(float(x))
def compensated(a,b):return math.fsum(float(x)*float(y) for x,y in zip(a,b))

def ledger():
    result=[]
    for family,name in [('original','avl-fixed-band-grid-runs.json'),('aligned','avl-band-aligned-grid-runs.json')]:
        for run in json.loads((FOLDER/name).read_text())['runs']:
            if run['alpha']!=1:continue
            p=ROOT/run['directory'];rows=np.loadtxt(p/'uv-rows.txt',skiprows=2)
            forces=np.loadtxt(p/'uv-forces.txt',skiprows=2);finite=forces[forces[:,2]==0]
            a=rows[:,7:10];b=rows[:,10:13];point=finite[:,3:6];g=rows[:,13]/math.sin(math.pi/180)
            dy=b[:,1]-a[:,1];f=2*g[:,None]*np.cross(np.array([1.,0.,0.]),b-a)
            point_values=bands.aggregate(point,f);couple=wrench.convert(a,b,point,f)
            line=wrench.observe(a,b,f,couple)
            w=np.column_stack([2*dy/bands.S,-2*(point[:,0]-.375)*dy/(bands.S*bands.C)])
            q=np.array([compensated(w[:,j],g) for j in range(2)])
            strip_records=[];strip_bands=np.zeros(20)
            for strip in sorted(set(finite[:,1].astype(int))):
                ix=finite[:,1].astype(int)==strip
                ys=float(point[ix,1][0]);owner=min(int(np.searchsorted(bands.BOUNDS,ys,side='right')-1),19)
                amount=2*math.fsum(g[ix]*dy[ix])/bands.S
                strip_bands[owner]+=amount
                strip_records.append(dict(strip=int(strip),y=ys,chord_gamma_sum=math.fsum(g[ix]),
                    weighted_gamma_dy_sum=math.fsum(g[ix]*dy[ix]),chord_count=int(ix.sum()),band=owner))
            point_bands=point_values['force'][:,2]/bands.S
            result.append(dict(family=family,case=run['case'],m=run['m'],n=run['n'],
                lift=float(q[0]),moment_y=float(q[1]),lift_bands=point_bands.tolist(),
                reconstructed_lift_bands=(line['force'][:,2]/bands.S).tolist(),strips=strip_records,
                scalar_force_difference=abs(float(q[0])-math.fsum(point_bands)),
                scalar_moment_difference=abs(float(q[1])-math.fsum(point_values['moment'][:,1])/(bands.S*bands.C)),
                strip_band_difference=float(np.max(abs(strip_bands-point_bands))),
                line_couple_moment_difference=abs(float(q[1])-math.fsum(line['moment'][:,1])/(bands.S*bands.C)),
                point_vs_line_band_difference=float(np.max(abs(point_bands-line['force'][:,2]/bands.S)))))
    return result

def analyze(run):
    folder=ROOT/run['directory']/'out'
    summary={k:float(v) for k,v in (line.split() for line in (folder/'us-summary.txt').read_text().splitlines())}
    n=int(summary['N']);vectors=np.loadtxt(folder/'us-vectors.txt',skiprows=1)
    g=vectors[:,1:5];r=vectors[:,5:9]
    responses=np.loadtxt(folder/'us-responses.txt',skiprows=1)
    adj=np.loadtxt(folder/'us-adjoints.txt',skiprows=1).reshape((22,n,5))
    cross=np.loadtxt(folder/'us-matrix-cross.txt',skiprows=1)
    rows=np.loadtxt(folder/'uv-rows.txt',skiprows=2);forces=np.loadtxt(folder/'uv-forces.txt',skiprows=2)
    p=forces[forces[:,2]==0,3:6];dy=rows[:,11]-rows[:,8]
    weights=np.zeros((22,n));weights[0]=2*dy/bands.S
    weights[1]=-2*(p[:,0]-.375)*dy/(bands.S*bands.C)
    owner=np.minimum(np.searchsorted(bands.BOUNDS,p[:,1],side='right')-1,19)
    weights[owner+2,np.arange(n)]=weights[0]
    sel=np.array([0,1,n//4-1,n//2-1,n//2,3*n//4-1,n-2,n-1])
    rowmat=np.full((8,n),np.nan);colmat=rowmat.copy()
    index={int(v):j for j,v in enumerate(sel)}
    for ii,jj,v in cross:
        i=int(ii)-1;j=int(jj)-1
        if i in index:rowmat[index[i],j]=v
        if j in index:colmat[index[j],i]=v
    assert np.all(np.isfinite(rowmat)) and np.all(np.isfinite(colmat))
    decimal_primal=[];decimal_adjoint=[]
    with localcontext() as ctx:
        ctx.prec=80
        for k,i in enumerate(sel):
            ar=[dec(v) for v in rowmat[k]];ac=[dec(v) for v in colmat[k]]
            for j in range(4):
                exact=-dec(rows[i,6])-sum((v*dec(x) for v,x in zip(ar,g[:,j])),Decimal(0))
                decimal_primal.append(dict(i=int(i+1),representation=j,residual=float(exact),difference=float(abs(exact-dec(r[i,j])))))
            for j in range(22):
                exact=dec(adj[j,i,2])-sum((v*dec(x) for v,x in zip(ac,adj[j,:,3])),Decimal(0))
                decimal_adjoint.append(dict(i=int(i+1),output=j+1,residual=float(exact),difference=float(abs(exact-dec(adj[j,i,4])))))
    small=None
    if n<=192:
        a=np.fromfile(folder/'us-matrix.bin',dtype='<f8').reshape((n,n),order='F')
        inverse=scipy.linalg.inv(a,check_finite=True)
        anorms=[float(np.linalg.norm(a,1)),float(np.linalg.norm(a,np.inf))]
        invnorms=[float(np.linalg.norm(inverse,1)),float(np.linalg.norm(inverse,np.inf))]
        estimates=[summary['AINV1_EST'],summary['AINVINF_EST']]
        norms=[summary['ANORM1'],summary['ANORMINF']]
        small=dict(matrix_binary_format='little-endian IEEE float64, column-major',a_norms=anorms,
            inverse_norms=invnorms,estimate_to_explicit=[x/y for x,y in zip(estimates,invnorms)],
            a_norm_relative=max(abs(x/y-1) for x,y in zip(norms,anorms)),
            independent_all_row_residual_difference=max(abs(-rows[i,6]-compensated(a[i],g[:,j])-r[i,j]) for i in range(n) for j in range(4)))
    q=responses[:,1:5];change=q[:,1:]-q[:,0,None]
    relative=[float(np.max(abs(change[0]))/abs(q[0,0])),float(np.max(abs(change[1]))/max(abs(q[1,0]),.1)),
        float(max(np.linalg.norm(change[2:,j]) for j in range(3))/np.linalg.norm(q[2:,0]))]
    backward=[float(np.max(abs(r[:,j]))/(summary['ANORMINF']*np.max(abs(g[:,j]))+np.max(abs(rows[:,6])))) for j in range(4)]
    perturb_abs=abs(responses[:,8]-responses[:,9])
    perturb_rel=perturb_abs/np.maximum(np.maximum(abs(responses[:,8]),abs(responses[:,9])),1e-300)
    recompute=np.array([[compensated(weights[k],g[:,j]) for j in range(4)] for k in range(22)])
    return dict(family=run['family'],case=run['case'],m=run['m'],n=run['n'],summary=summary,
        small_matrix=small,decimal_primal=decimal_primal,decimal_adjoint=decimal_adjoint,
        backward_errors=backward,primal_residual_maxima=np.max(abs(r),axis=0).tolist(),
        output_relative_changes=dict(zip(['lift','moment_y','lift_bands'],relative)),
        output_absolute_changes=dict(lift=float(max(abs(change[0]))),moment_y=float(max(abs(change[1]))),
            lift_bands=float(max(np.linalg.norm(change[2:,j]) for j in range(3)))),
        qbase=q[:,0].tolist(),max_weight_difference=float(np.max(abs(weights-adj[:,:,2]))),
        max_response_recompute_difference=float(np.max(abs(recompute-q))),
        max_adjoint_relative_residual=float(np.max(responses[:,7])),
        perturbation_abs_max=float(np.max(perturb_abs)),perturbation_relative_max=float(np.max(perturb_rel)),
        perturbation_check=bool(np.all((perturb_abs<=1e-13)|(perturb_rel<=1e-9))),
        rhs_perturbation_output_ranges=[float(min(responses[:,9])),float(max(responses[:,9]))],
        adjoint_base_correction_max=float(np.max(abs(responses[:,5]))),
        correction_vs_ir2_difference=float(np.max(abs(change[:,2]-responses[:,5]))))

def main():
    begin=time.monotonic();path=FOLDER/'avl-solve-sensitivity-runs.json';runs=json.loads(path.read_text())
    assert runs['completed']
    records=[analyze(r) for r in runs['runs']];grids=ledger()
    comparison=[]
    for rec in records:
        if rec['m'] not in [16,32,12,24] or rec['n'] not in [96,160]:continue
        here=next(r for r in grids if all(r[k]==rec[k] for k in ['family','case','m','n']))
        other=next(r for r in grids if r['family']==rec['family'] and r['case']==rec['case'] and r['n']==rec['n'] and r['m']==(rec['m']//2))
        diffs={k:abs(here[k]-other[k]) for k in ['lift','moment_y']}
        diffs['lift_bands']=float(np.linalg.norm(np.array(here['lift_bands'])-other['lift_bands']))
        comparison.append(dict(family=rec['family'],case=rec['case'],m=rec['m'],n=rec['n'],grid_changes=diffs,
            solve_change_to_grid_change={k:rec['output_absolute_changes'][k]/max(v,1e-12) for k,v in diffs.items()}))
    small=[r['small_matrix'] for r in records if r['small_matrix'] is not None]
    old=json.loads((FOLDER/'avl-linear-load-wrench-result.json').read_text())
    prior=max(max(abs(r[k]-next(x for x in old['linear'] if all(x[j]==r[j] for j in ['case','m','n']))[k]) for k in ['lift','moment_y']) for r in grids if r['family']=='original')
    checks=dict(eight_systems_unchanged=len(records)==8 and all(all(r['unchanged'].values()) for r in runs['runs']),
        norm_calibration=len(small)==2 and all(s['a_norm_relative']<=1e-12 and all(.1<=v<=1+1e-10 for v in s['estimate_to_explicit']) for s in small),
        decimal_residuals=all(max(v['difference'] for v in r['decimal_primal']+r['decimal_adjoint'])<=1e-13 for r in records),
        small_independent_residuals=all(s['independent_all_row_residual_difference']<=1e-13 for s in small),
        output_weights_and_recompute=all(max(r['max_weight_difference'],r['max_response_recompute_difference'])<=1e-11 for r in records),
        same_system_output_precision=all(max(r['output_relative_changes'].values())<=1e-10 for r in records),
        adjoint_residuals=all(r['max_adjoint_relative_residual']<=1e-11 for r in records),
        rhs_perturbation_identity=all(r['perturbation_check'] for r in records),
        circulation_load_ledger=len(grids)==40 and prior<=1e-11 and all(max(r[k] for k in ['scalar_force_difference','scalar_moment_difference','strip_band_difference','line_couple_moment_difference'])<=1e-11 for r in grids),
        aligned_observation=all(r['point_vs_line_band_difference']<=1e-11 for r in grids if r['family']=='aligned'),
        resource_budget=runs['seconds']<900 and runs['total_bytes']<256*1024**2 and max(r['peak_resident_bytes'] for r in runs['runs'])<2*1024**3 and time.monotonic()-begin<180 and memory.memory_bytes()<512*1024**2)
    result=dict(protocol='rd/aero/avl-reference/records/avl-solve-sensitivity-experiment.md',source_sha256=digest(Path(__file__)),
        dependency_sha256={name:digest(rd_paths.EXPERIMENTS/name) for name in ['check_avl_fixed_bands.py','check_avl_linear_load_wrench.py','check_q9_plate_convergence.py']},
        runs_sha256=digest(path),parent_linear_result_sha256=digest(FOLDER/'avl-linear-load-wrench-result.json'),
        records=records,circulation_ledger=grids,grid_vs_solve=comparison,previous_linear_difference=prior,
        checks=checks,all_criteria_pass=all(checks.values()),reference_convergence_accepted=False,
        qualification='Same stored float64 discrete systems only. Norm estimates are lower estimates, not certified bounds; matrix assembly and discretization error remain outside this result.',
        seconds=time.monotonic()-begin,peak_process_working_bytes=memory.memory_bytes())
    (FOLDER/'avl-solve-sensitivity-result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],peak=result['peak_process_working_bytes'],
        max_relative_output_change=max(max(r['output_relative_changes'].values()) for r in records),
        max_solve_grid_ratio=max(max(r['solve_change_to_grid_change'].values()) for r in comparison)),indent=2))

if __name__=='__main__':main()
