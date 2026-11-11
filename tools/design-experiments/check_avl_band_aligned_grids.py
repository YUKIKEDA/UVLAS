"""Conservative point-load bands and predeclared reference convergence."""
from pathlib import Path
import hashlib
import json
import math
import time
import numpy as np
import check_q9_plate_convergence as memory
from analyze_avl_source_measurement_corrected import fields,coefficients

ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification';START=time.monotonic()
BOUNDS=np.array([-4+.4*j for j in range(21)])
ORIGIN=np.array([.375,0.,0.]);S=8.4;C=1.1142857142857143;B=8.
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sum3(values):return np.array([math.fsum(values[:,j]) for j in range(3)])

def aggregate(points,force):
    p=np.asarray(points);f=np.asarray(force)
    if p.shape!=f.shape or p.ndim!=2 or p.shape[1]!=3:raise ValueError('invalid_shape')
    if not np.all(np.isfinite(p)) or not np.all(np.isfinite(f)):raise ValueError('nonfinite')
    y=p[:,1]
    if np.any(y<BOUNDS[0]-1e-9) or np.any(y>BOUNDS[-1]+1e-9):raise ValueError('outside_fixed_bands')
    # Clamp classification only; moments retain the actual input coordinates.
    classified=np.clip(y,BOUNDS[0],BOUNDS[-1])
    owner=np.searchsorted(BOUNDS,classified,side='right')-1;owner[owner==20]=19
    moments=np.cross(p-ORIGIN,f)
    bands_f=np.array([sum3(f[owner==j]) for j in range(20)])
    bands_m=np.array([sum3(moments[owner==j]) for j in range(20)])
    return dict(force=bands_f,moment=bands_m,owner=owner,counts=np.bincount(owner,minlength=20),
        clamped=int(np.count_nonzero(classified!=y)),maximum_outside=float(np.max(abs(classified-y),initial=0.)))

def manufactured():
    p=np.array([[-.3,-4.,0.],[.2,-3.6,.1],[.7,0.,-.2],[1.1,4.,.3]])
    f=np.array([[1.,2.,3.],[-2.,1.,4.],[3.,-1.,2.],[1.,3.,-2.]])
    before=digest_array(p)+digest_array(f);r=aggregate(p,f)
    moved=p.copy();moved[:,1]=(BOUNDS[r['owner']]+BOUNDS[r['owner']+1])/2
    negative=float(np.linalg.norm(sum3(np.cross(moved-ORIGIN,f))-sum3(r['moment'])))
    failures=[];buffer=np.full((20,3),17.);buffer_hash=digest_array(buffer)
    for label,pp,ff in [('outside',np.array([[0.,4+1e-6,0.]]),np.ones((1,3))),
        ('nonfinite_point',np.array([[0.,np.nan,0.]]),np.ones((1,3))),
        ('nonfinite_force',np.zeros((1,3)),np.array([[np.inf,1.,1.]]))]:
        try:buffer[:]=aggregate(pp,ff)['force'];refused=False
        except ValueError:refused=True
        failures.append(dict(kind=label,refused=refused,buffer_unchanged=digest_array(buffer)==buffer_hash))
    edge=aggregate(np.array([[0.,4+5e-10,0.]]),np.ones((1,3)))
    checks=dict(boundary_ownership=r['owner'].tolist()==[0,1,10,19],
        force=np.max(abs(sum3(r['force'])-[3.,5.,7.]))<=1e-13,
        moment=np.max(abs(sum3(r['moment'])-[-35.6,3.025,-6.875]))<=1e-13,
        negative_moved_points=negative>1e-3,failures=all(x['refused'] and x['buffer_unchanged'] for x in failures),
        inputs_unchanged=before==digest_array(p)+digest_array(f),classification_tolerance=edge['clamped']==1 and edge['owner'][0]==19)
    return dict(checks=checks,negative_moment_difference=negative,failures=failures,
        force=sum3(r['force']).tolist(),moment=sum3(r['moment']).tolist())

def digest_array(a):return hashlib.sha256(a.tobytes()).hexdigest()

def estimate(values,scale,vector=False):
    values=np.asarray(values);differences=np.diff(values,axis=0)
    magnitude=np.linalg.norm(differences,axis=1) if vector else abs(differences)
    if vector:
        # The band width is common to every component, so it cancels in relative metrics.
        magnitude=magnitude/np.sqrt(.4)
    result=dict(differences=magnitude.tolist(),status='not_asymptotic',estimate=None,orders=None)
    if np.max(magnitude)<=1e-10*scale:
        result.update(status='roundoff_candidate',estimate=float(4*np.max(magnitude)))
        return result
    if np.any(magnitude<=0):return result
    orders=np.log2(magnitude[:-1]/magnitude[1:]);result['orders']=orders.tolist()
    if vector:
        cosines=[float(np.dot(differences[j],differences[j+1])/(np.linalg.norm(differences[j])*np.linalg.norm(differences[j+1]))) for j in range(2)]
        result['difference_cosines']=cosines;direction=min(cosines)>=.9
    else:direction=bool(np.all(differences[:-1]*differences[1:]>0))
    if direction and np.all(np.diff(magnitude)<0) and .5<=orders[-1]<=6 and abs(orders[-1]-orders[-2])<=.5:
        result.update(status='asymptotic_candidate',estimate=float(1.25*magnitude[-1]/(2**orders[-1]-1)))
    return result

def main():
    parent_path=FOLDER/'avl-band-aligned-grid-runs.json';parent=json.loads(parent_path.read_text());assert parent['completed']
    cases=[];lookup={}
    for run in parent['runs']:
        directory=ROOT/run['directory'];rows=np.loadtxt(directory/'uv-rows.txt',skiprows=2)
        forces=np.loadtxt(directory/'uv-forces.txt',skiprows=2);total=fields(directory/'total.txt');coefs=coefficients(total)
        m,n=run['m'],run['n'];size=m*n
        assert rows.shape==(size,21) and forces.shape==(3*size,9)
        r=aggregate(forces[:,3:6],forces[:,6:9])
        fcoef=r['force']/S;mcoef=r['moment']/(S*np.array([B,C,B]))
        force_difference=float(np.max(abs(sum3(fcoef)-coefs['force'])))
        moment_difference=float(np.max(abs(sum3(mcoef)-coefs['moment'])))
        y=rows[:,[8,11]].reshape(n,m,2);intervals=y[:,0,:]
        interval_error=float(np.max(abs(y-intervals[:,None,:])))
        gap=float(np.max(abs(intervals[1:,0]-intervals[:-1,1]),initial=0.))
        interval_eligible=interval_error<=1e-12 and gap<=1e-12 and np.all(intervals[:,1]>intervals[:,0]) and max(abs(intervals[0,0]+4),abs(intervals[-1,1]-4))<=1e-9
        explicit_error=None
        if not run['reused']:
            strips=np.loadtxt(directory/'uv-strips.txt',skiprows=2)
            explicit_error=float(np.max(abs(strips[:,1:3]-intervals)))
            interval_eligible=bool(interval_eligible and explicit_error<=1e-13 and np.all(strips[:,3]==m))
        residual=float(np.max(abs(rows[:,20]))/max(np.max(np.linalg.norm(rows[:,14:17],axis=1)),1e-300))
        alpha=np.deg2rad(run['alpha']);lift=np.array([-np.sin(alpha),0.,np.cos(alpha)])
        lift_bands=fcoef@lift
        native_counts={(int(q[0]),int(q[2])) for q in forces}
        counts_ok=native_counts=={(i,k) for i in range(1,size+1) for k in range(3)}
        item=dict(case=run['case'],m=m,n=n,alpha=run['alpha'],reused=run['reused'],
            total=coefs,force_bands=fcoef.tolist(),moment_bands=mcoef.tolist(),lift_bands=lift_bands.tolist(),
            counts=r['counts'].tolist(),clamped=r['clamped'],maximum_outside=r['maximum_outside'],
            force_difference=force_difference,moment_difference=moment_difference,residual_relative=residual,
            count_eligible=counts_ok and int(sum(r['counts']))==3*size,
            input_strip_intervals=intervals.tolist(),interval_error=interval_error,gap_error=gap,
            explicit_interval_error=explicit_error,interval_eligible=bool(interval_eligible),
            export_unchanged=run.get('after_export_equal',True))
        cases.append(item);lookup[(run['case'],m,n,run['alpha'])]=item
    slopes=[];slopes_lookup={}
    for ident in ['FW-V01','FW-V02']:
        meshes=sorted({(r['m'],r['n']) for r in cases if r['case']==ident})
        for m,n in meshes:
            for delta in [1,2]:
                plus=lookup[(ident,m,n,delta)];minus=lookup[(ident,m,n,-delta)];denom=2*np.deg2rad(delta)
                value=dict(case=ident,m=m,n=n,delta=delta,
                    lift=(plus['total']['lift']-minus['total']['lift'])/denom,
                    moment_y=(plus['total']['moment'][1]-minus['total']['moment'][1])/denom,
                    lift_bands=((np.array(plus['lift_bands'])-minus['lift_bands'])/denom).tolist())
                slopes.append(value);slopes_lookup[(ident,m,n,delta)]=value
    convergence=[];angle=[]
    for ident in ['FW-V01','FW-V02']:
        for delta in [1,2]:
            finest=slopes_lookup[(ident,24,160,delta)]
            record=dict(case=ident,delta=delta,scalar={},distribution={})
            for name in ['lift','moment_y','lift_bands']:
                vector=name=='lift_bands'
                final=np.array(finest[name]);scale=float(np.linalg.norm(final)/np.sqrt(.4)) if vector else max(abs(float(final)),.1 if name=='moment_y' else 1e-300)
                chord=[slopes_lookup[(ident,m,160,delta)][name] for m in [3,6,12,24]]
                span=[slopes_lookup[(ident,24,n,delta)][name] for n in [20,40,80,160]]
                ce=estimate(chord,scale,vector);se=estimate(span,scale,vector)
                relative=None if ce['estimate'] is None or se['estimate'] is None else (ce['estimate']+se['estimate'])/scale
                details=dict(scale=scale,chord=ce,span=se,combined_relative_estimate=relative,
                    eligible=relative is not None and relative<=(.005 if vector else .003))
                (record['distribution'] if vector else record['scalar'])[name]=details
            convergence.append(record)
        one=slopes_lookup[(ident,24,160,1)];two=slopes_lookup[(ident,24,160,2)]
        differences={name:abs(one[name]-two[name])/max(abs(one[name]),.1 if name=='moment_y' else 1e-300) for name in ['lift','moment_y']}
        angle.append(dict(case=ident,differences=differences,eligible=max(differences.values())<=.001))
    fixture=manufactured()
    checks=dict(conservative_force=all(r['force_difference']<=1e-10 for r in cases),
        conservative_moment=all(r['moment_difference']<=1e-10 for r in cases),
        complete_contributions=all(r['count_eligible'] for r in cases),strip_tessellation=all(r['interval_eligible'] for r in cases),
        residual=all(r['residual_relative']<=1e-11 for r in cases),export_regression=all(r['export_unchanged'] for r in cases),
        manufactured=all(fixture['checks'].values()),
        scalar_convergence=all(v['eligible'] for r in convergence for v in r['scalar'].values()),
        angle_sensitivity=all(r['eligible'] for r in angle),
        distribution_convergence=all(v['eligible'] for r in convergence for v in r['distribution'].values()))
    result=dict(protocol='docs/verification/avl-band-aligned-grid-experiment.md',source_sha256=digest(Path(__file__)),
        dependency_sha256={name:digest(ROOT/'tools/design-experiments'/name) for name in ['check_q9_plate_convergence.py','analyze_avl_source_measurement_corrected.py']},
        parent_result_sha256=digest(parent_path),bounds=BOUNDS.tolist(),observation='point-force measure; moment kept about original O',
        cases=cases,slopes=slopes,convergence=convergence,angle=angle,manufactured=fixture,checks=checks,
        all_criteria_pass=all(checks.values()),seconds=time.monotonic()-START,peak_process_working_bytes=memory.memory_bytes(),
        limits='reference side only; no UVLAS load comparison; band distribution is a declared point-load observation')
    (FOLDER/'avl-band-aligned-convergence-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else str(v))+'\n')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],peak_bytes=result['peak_process_working_bytes'],
        max_force=max(r['force_difference'] for r in cases),max_moment=max(r['moment_difference'] for r in cases))),flush=True)
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
