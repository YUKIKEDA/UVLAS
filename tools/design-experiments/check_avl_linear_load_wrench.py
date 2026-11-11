"""Predeclared exact linear load and conservative observation experiment."""
from pathlib import Path
import hashlib
import json
import math
import time
import numpy as np
import check_q9_plate_convergence as memory
import check_avl_fixed_bands as old

ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification';START=time.monotonic()
BOUNDS=old.BOUNDS;O=old.ORIGIN;S=old.S;C=old.C;B=old.B
PARTIAL=[]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(v):return np.linalg.norm(v,axis=-1)
def sum3(v):return np.array([math.fsum(v[:,j]) for j in range(3)])

def convert(a,b,p,f):
    arrays=[np.asarray(x,dtype=float) for x in [a,b,p,f]];a,b,p,f=arrays
    if any(x.ndim!=2 or x.shape[1]!=3 or x.shape!=a.shape for x in arrays):raise ValueError('invalid_shape')
    if not all(np.all(np.isfinite(x)) for x in arrays):raise ValueError('nonfinite')
    edge=b-a;length2=np.sum(edge*edge,axis=1)
    if np.any(length2<=1e-24):raise ValueError('degenerate_support')
    t=np.sum((p-a)*edge,axis=1)/length2
    distance=norm(p-a-t[:,None]*edge)
    if np.any((t<-1e-9)|(t>1+1e-9)|(distance>1e-9)):raise ValueError('point_off_support')
    return np.cross(p-(a+b)/2,f)

def observe(a,b,f,couple,origin=O,axis=np.array([0.,1.,0.]),span_origin=np.zeros(3),bounds=BOUNDS):
    arrays=[np.asarray(x,dtype=float) for x in [a,b,f,couple]];a,b,f,couple=arrays
    axis=np.asarray(axis);origin=np.asarray(origin);span_origin=np.asarray(span_origin);bounds=np.asarray(bounds)
    if any(x.ndim!=2 or x.shape[1]!=3 or x.shape!=a.shape for x in arrays):raise ValueError('invalid_shape')
    if not all(np.all(np.isfinite(x)) for x in arrays+[axis,origin,span_origin,bounds]):raise ValueError('nonfinite')
    if abs(norm(axis)-1)>1e-12 or len(bounds)<2 or np.any(np.diff(bounds)<=0):raise ValueError('invalid_observer')
    sa=(a-span_origin)@axis;sb=(b-span_origin)@axis
    if np.any(sa<bounds[0]-1e-9)|np.any(sa>bounds[-1]+1e-9)|np.any(sb<bounds[0]-1e-9)|np.any(sb>bounds[-1]+1e-9):raise ValueError('support_outside')
    sa=np.clip(sa,bounds[0],bounds[-1]);sb=np.clip(sb,bounds[0],bounds[-1])
    if np.any(abs(sb-sa)<=1e-12*np.maximum(norm(b-a),1.)):raise ValueError('zero_span_width')
    swap=sa>sb
    aa=np.where(swap[:,None],b,a);bb=np.where(swap[:,None],a,b)
    low=np.minimum(sa,sb);high=np.maximum(sa,sb);width=high-low
    left=np.clip((bounds[:-1][None,:]-low[:,None])/width[:,None],0.,1.)
    right=np.clip((bounds[1:][None,:]-low[:,None])/width[:,None],0.,1.)
    weight=np.maximum(0.,right-left)
    if np.any(abs(weight.sum(axis=1)-1)>1e-12):raise ValueError('incomplete_support')
    centers=aa[:,None,:]+((left+right)/2)[:,:,None]*(bb-aa)[:,None,:]
    ff=weight[:,:,None]*f[:,None,:]
    mm=weight[:,:,None]*(np.cross(centers-origin,f[:,None,:])+couple[:,None,:])
    return dict(force=np.array([sum3(ff[:,j]) for j in range(len(bounds)-1)]),
        moment=np.array([sum3(mm[:,j]) for j in range(len(bounds)-1)]),
        weight_error=float(np.max(abs(weight.sum(axis=1)-1))),
        couple_sum=sum3(couple),couple_norm_sum=float(math.fsum(norm(couple))))

def quadrature(a,b,f,couple,origin=O,axis=np.array([0.,1.,0.]),span_origin=np.zeros(3),bounds=BOUNDS):
    q,w=np.polynomial.legendre.leggauss(8);q=(q+1)/2;w=w/2
    sa=(a-span_origin)@axis;sb=(b-span_origin)@axis;ds=sb-sa
    force=[];moment=[]
    for lo,hi in zip(bounds[:-1],bounds[1:]):
        t0=(lo-sa)/ds;t1=(hi-sa)/ds
        left=np.maximum(0.,np.minimum(t0,t1));right=np.minimum(1.,np.maximum(t0,t1));length=np.maximum(0.,right-left)
        nodes=left[:,None]+length[:,None]*q[None,:]
        positions=a[:,None,:]+nodes[:,:,None]*(b-a)[:,None,:]
        ww=length[:,None]*w[None,:]
        integrand=np.cross(positions-origin,f[:,None,:])+couple[:,None,:]
        force.append(sum3(np.sum(ww[:,:,None]*f[:,None,:],axis=1)))
        moment.append(sum3(np.sum(ww[:,:,None]*integrand,axis=1)))
    return dict(force=np.array(force),moment=np.array(moment))

def manufactured():
    a=np.array([[0.,-1.,0.]]);b=np.array([[2.,1.,1.]]);p=a+.25*(b-a);f=np.array([[1.,2.,3.]])
    couple=convert(a,b,p,f);r=observe(a,b,f,couple);ind=quadrature(a,b,f,couple)
    expected_m=sum3(np.cross(p-O,f));d=np.array([.3,-.2,.4]);shift=observe(a,b,f,couple,origin=O+d)
    rz=np.array([[np.cos(.4),-np.sin(.4),0],[np.sin(.4),np.cos(.4),0],[0,0,1.]])
    ry=np.array([[np.cos(-.3),0,np.sin(-.3)],[0,1.,0],[-np.sin(-.3),0,np.cos(-.3)]])
    rotation=rz@ry;translation=np.array([.3,-.2,.5])
    transformed=observe(a@rotation.T+translation,b@rotation.T+translation,f@rotation.T,couple@rotation.T,
        origin=O@rotation.T+translation,axis=rotation@np.array([0.,1.,0.]),span_origin=translation)
    middle=a+.3*(b-a)
    split=observe(np.vstack([a,middle]),np.vstack([middle,b]),np.vstack([.3*f,.7*f]),np.vstack([.3*couple,.7*couple]))
    reverse=observe(b,a,f,couple);without=observe(a,b,f,np.zeros_like(f))
    v=np.array([.1,.2,-.3]);omega=np.array([.4,-.2,.5])
    original_work=float(np.sum(f*(v+np.cross(omega,p-O))))
    rigid_work=float(sum3(r['force'])@v+sum3(r['moment'])@omega)
    nonrigid_original=float(f[0,2]*p[0,1]**2)
    q,w=np.polynomial.legendre.leggauss(8);yy=a[0,1]+(q+1)/2*(b[0,1]-a[0,1])
    nonrigid_spread=float(f[0,2]*np.dot(w/2,yy**2))
    failures=[];buffer=np.full((20,3),13.);before=hashlib.sha256(buffer.tobytes()).hexdigest()
    for name,aa,bb,ff,axis,bounds in [
        ('outside',a+np.array([0.,10.,0.]),b+np.array([0.,10.,0.]),f,np.array([0.,1.,0.]),BOUNDS),
        ('zero_span',a,a+np.array([1.,0.,0.]),f,np.array([0.,1.,0.]),BOUNDS),
        ('nonfinite',a,b,np.array([[np.nan,1.,2.]]),np.array([0.,1.,0.]),BOUNDS),
        ('axis',a,b,f,np.zeros(3),BOUNDS),('bounds',a,b,f,np.array([0.,1.,0.]),BOUNDS[::-1])]:
        try:buffer[:]=observe(aa,bb,ff,couple,axis=axis,bounds=bounds)['force'];refused=False
        except ValueError:refused=True
        failures.append(dict(kind=name,refused=refused,buffer_unchanged=hashlib.sha256(buffer.tobytes()).hexdigest()==before))
    try:convert(a,b,p+np.array([0.,0.,.1]),f);off_refused=False
    except ValueError:off_refused=True
    checks=dict(independent_line_quadrature=max(np.max(abs(r[k]-ind[k])) for k in ['force','moment'])<=1e-12,
        global_wrench=max(np.max(abs(sum3(r['force'])-f[0])),np.max(abs(sum3(r['moment'])-expected_m)))<=1e-12,
        origin_shift=np.max(abs(shift['moment']-(r['moment']-np.cross(d,r['force']))))<=1e-12,
        rigid_covariance=max(np.max(abs(transformed[k]-r[k]@rotation.T)) for k in ['force','moment'])<=1e-12,
        parameter_reversal=max(np.max(abs(reverse[k]-r[k])) for k in ['force','moment'])<=1e-12,
        subdivision=max(np.max(abs(split[k]-r[k])) for k in ['force','moment'])<=1e-12,
        rigid_work=abs(original_work-rigid_work)<=1e-12,
        negative_omit_couple=norm(sum3(without['moment'])-expected_m)>1e-3,
        nonrigid_force_work_changes=abs(nonrigid_original-nonrigid_spread)>1e-3,
        failures=all(x['refused'] and x['buffer_unchanged'] for x in failures) and off_refused)
    return dict(checks=checks,source_couple=couple.tolist(),negative_moment_difference=float(norm(sum3(without['moment'])-expected_m)),
        rigid_work_difference=abs(original_work-rigid_work),nonrigid_force_work_difference=abs(nonrigid_original-nonrigid_spread),failures=failures)

def analytic_distribution(parent):
    records=[]
    primitive=lambda y:.5*(y*np.sqrt(np.maximum(0.,1-(y/4)**2))+4*np.arcsin(np.clip(y/4,-1.,1.)))
    first=lambda y:-16/3*np.maximum(0.,1-(y/4)**2)**1.5
    exact_f=np.diff(primitive(BOUNDS));exact_m=np.diff(first(BOUNDS))
    for ident in ['FW-V01','FW-V02']:
        for n in [12,24,48,96]:
            c=next(c for c in parent['cases'] if c['case']==ident and c['m']==32 and c['n']==n and c['alpha']==1)
            intervals=np.array(c['input_strip_intervals']);amount=np.diff(primitive(intervals),axis=1)[:,0];moment=np.diff(first(intervals),axis=1)[:,0]
            a=np.zeros((n,3));b=a.copy();a[:,1]=intervals[:,0];b[:,1]=intervals[:,1]
            p=np.zeros_like(a);p[:,1]=moment/amount;f=np.zeros_like(a);f[:,2]=amount
            couple=convert(a,b,p,f);value=observe(a,b,f,couple)
            error_f=float(norm(value['force'][:,2]-exact_f)/norm(exact_f));error_m=float(norm(value['moment'][:,0]-exact_m)/norm(exact_m))
            conservation=max(abs(math.fsum(value['force'][:,2])-2*math.pi),abs(math.fsum(value['moment'][:,0])))
            records.append(dict(case=ident,n=n,force_relative=error_f,moment_relative=error_m,conservation=conservation,
                force_bands=value['force'][:,2].tolist(),moment_bands=value['moment'][:,0].tolist()))
    return records

def convergence(records):
    lookup={(r['case'],r['m'],r['n'],r.get('delta',0)):r for r in records};result=[]
    deltas=sorted({r.get('delta',0) for r in records})
    for ident in ['FW-V01','FW-V02']:
        for delta in deltas:
            finest=lookup[(ident,32,96,delta)];item=dict(case=ident,delta=delta,scalar={},distribution={})
            for name in ['lift','moment_y','lift_bands']:
                vector=name=='lift_bands';value=np.array(finest[name])
                scale=float(norm(value)/np.sqrt(.4)) if vector else max(abs(float(value)),.1 if name=='moment_y' else 1e-300)
                chord=[lookup[(ident,m,96,delta)][name] for m in [4,8,16,32]]
                span=[lookup[(ident,32,n,delta)][name] for n in [12,24,48,96]]
                ce=old.estimate(chord,scale,vector);se=old.estimate(span,scale,vector)
                ratio=None if ce['estimate'] is None or se['estimate'] is None else (ce['estimate']+se['estimate'])/scale
                (item['distribution'] if vector else item['scalar'])[name]=dict(chord=ce,span=se,combined_relative=ratio,eligible=ratio is not None and ratio<=(.005 if vector else .003))
            result.append(item)
    return result

def main():
    parent_path=FOLDER/'avl-fixed-band-convergence-result.json';parent=json.loads(parent_path.read_text())
    runs_path=FOLDER/'avl-fixed-band-grid-runs.json';raw=json.loads(runs_path.read_text())
    previous={(r['case'],r['m'],r['n'],r['alpha']):r for r in parent['cases']}
    numeric={};observations=[];quad_cases=[]
    for run in raw['runs']:
        folder=ROOT/run['directory'];rows=np.loadtxt(folder/'uv-rows.txt',skiprows=2);forces=np.loadtxt(folder/'uv-forces.txt',skiprows=2)
        a=rows[:,7:10];b=rows[:,10:13];finite=forces[forces[:,2]==0];other=forces[forces[:,2]!=0]
        p=finite[:,3:6];f=finite[:,6:9];couple=convert(a,b,p,f);recon=observe(a,b,f,couple)
        points=old.aggregate(other[:,3:6],other[:,6:9]);ff=(recon['force']+points['force'])/S;mm=(recon['moment']+points['moment'])/(S*np.array([B,C,B]))
        key=(run['case'],run['m'],run['n'],run['alpha']);native=previous[key]
        difference_f=float(np.max(abs(sum3(ff)-native['total']['force'])));difference_m=float(np.max(abs(sum3(mm)-native['total']['moment'])))
        angle=np.deg2rad(run['alpha']);lift=ff@np.array([-np.sin(angle),0.,np.cos(angle)])
        observations.append(dict(case=run['case'],m=run['m'],n=run['n'],alpha=run['alpha'],force_difference=difference_f,moment_difference=difference_m,
            force_bands=ff.tolist(),moment_bands=mm.tolist(),lift_bands=lift.tolist(),
            couple_sum_over_q=recon['couple_sum'].tolist(),couple_norm_sum_over_q=recon['couple_norm_sum']))
        numeric[key]=dict(gamma=rows[:,13],normal_x_max=float(np.max(abs(rows[:,4]))),a=a,b=b,p=p)
        if (run['m'],run['n']) in [(4,12),(8,24)] and run['alpha']==1:
            q=quadrature(a,b,f,couple)
            quad_cases.append(dict(case=run['case'],m=run['m'],n=run['n'],force_difference=float(np.max(abs(q['force']-recon['force']))/S),
                moment_difference=float(np.max(abs(q['moment']-recon['moment']))/(S*C))))
        PARTIAL.append(dict(case=run['case'],m=run['m'],n=run['n'],alpha=run['alpha']))
        if time.monotonic()-START>180:raise TimeoutError('postprocessing_limit')
    linear=[];decomposition=[];lookup={(r['case'],r['m'],r['n'],r['alpha']):r for r in observations}
    finite_slopes=[]
    for ident in ['FW-V01','FW-V02']:
        meshes=sorted({(r['m'],r['n']) for r in observations if r['case']==ident})
        for m,n in meshes:
            zero=numeric[(ident,m,n,0)];gz=[]
            for delta in [1,2]:
                plus=numeric[(ident,m,n,delta)];minus=numeric[(ident,m,n,-delta)]
                gz.append((plus['gamma']-minus['gamma'])/(2*np.sin(np.deg2rad(delta))))
            g=gz[0];a=zero['a'];b=zero['b'];p=zero['p'];f=2*g[:,None]*np.cross(np.array([1.,0.,0.]),b-a)
            couple=convert(a,b,p,f);r=observe(a,b,f,couple);pure_m=sum3(np.cross((a+b)/2-O,f))
            f0=sum3(f)/S;m0=sum3(np.cross(p-O,f))/(S*np.array([B,C,B]));q=quadrature(a,b,f,np.zeros_like(f))
            linear.append(dict(case=ident,m=m,n=n,delta=0,lift=float(f0[2]),moment_y=float(m0[1]),
                line_moment_y=float(pure_m[1]/(S*C)),couple_moment_y=float(sum3(couple)[1]/(S*C)),
                lift_bands=(r['force'][:,2]/S).tolist(),
                gamma_zero=float(np.max(abs(zero['gamma']))),normal_x_max=zero['normal_x_max'],
                gamma_sine_difference=float(np.max(abs(gz[0]-gz[1])/np.maximum(abs(gz[0]),1.))),
                line_quadrature_force=float(np.max(abs(sum3(q['force'])/S-f0))),
                line_quadrature_moment=float(np.max(abs(sum3(q['moment'])-pure_m))/(S*C))))
            for delta in [1,2]:
                d=np.deg2rad(delta);plus=previous[(ident,m,n,delta)];minus=previous[(ident,m,n,-delta)]
                native_l=(plus['total']['lift']-minus['total']['lift'])/(2*d)
                native_m=(plus['total']['moment'][1]-minus['total']['moment'][1])/(2*d)
                x2=(plus['total']['force'][0]+minus['total']['force'][0])/(2*np.sin(d)**2)
                predict_l=np.sin(d)*np.cos(d)**2/d*f0[2]-np.sin(d)**3/d*x2
                predict_m=np.sin(d)*np.cos(d)/d*m0[1]
                decomposition.append(dict(case=ident,m=m,n=n,delta=delta,native_lift=float(native_l),native_moment=float(native_m),
                    predicted_lift=float(predict_l),predicted_moment=float(predict_m),lift_difference=float(abs(predict_l-native_l)),moment_difference=float(abs(predict_m-native_m)),
                    relative_finite_lift=float(abs(native_l-f0[2])/abs(f0[2])),relative_finite_moment=float(abs(native_m-m0[1])/max(abs(m0[1]),.1))))
                pp=lookup[(ident,m,n,delta)];nn=lookup[(ident,m,n,-delta)]
                finite_slopes.append(dict(case=ident,m=m,n=n,delta=delta,lift=float(native_l),moment_y=float(native_m),
                    lift_bands=((np.array(pp['lift_bands'])-nn['lift_bands'])/(2*d)).tolist()))
    unit=manufactured();analytic=analytic_distribution(parent)
    linear_conv=convergence(linear);finite_conv=convergence(finite_slopes)
    checks=dict(linear_input=all(r['gamma_zero']<=1e-12 and r['normal_x_max']<=1e-12 and r['gamma_sine_difference']<=1e-10 for r in linear),
        finite_angle_decomposition=all(max(r['lift_difference'],r['moment_difference'])<=1e-10 for r in decomposition),
        independent_linear_line=all(max(r['line_quadrature_force'],r['line_quadrature_moment'])<=1e-11 for r in linear),
        manufactured_wrench=all(unit['checks'].values()),
        analytic_distribution=all(r['force_relative']<=.01 and r['moment_relative']<=.01 for r in analytic if r['n']==96) and all(r['conservation']<=1e-12 for r in analytic),
        aerodynamic_wrench_conservation=all(max(r['force_difference'],r['moment_difference'])<=1e-10 for r in observations),
        independent_aerodynamic_line=all(max(r['force_difference'],r['moment_difference'])<=1e-11 for r in quad_cases),
        linear_scalar_convergence=all(v['eligible'] for r in linear_conv for v in r['scalar'].values()),
        linear_distribution_convergence=all(v['eligible'] for r in linear_conv for v in r['distribution'].values()),
        finite_distribution_convergence=all(v['eligible'] for r in finite_conv for v in r['distribution'].values()))
    source=Path(__file__).resolve()
    dependencies=['check_q9_plate_convergence.py','check_avl_fixed_bands.py','analyze_avl_source_measurement_corrected.py']
    result=dict(protocol='docs/verification/avl-linear-load-wrench-experiment.md',source_sha256=digest(source),
        dependency_sha256={name:digest(source.parent/name) for name in dependencies},parent_result_sha256=digest(parent_path),runs_sha256=digest(runs_path),
        observation_id='span-wrench-P0-v1',original_point_observation_unchanged=True,linear=linear,decomposition=decomposition,
        observations=observations,finite_slopes=finite_slopes,linear_convergence=linear_conv,finite_convergence=finite_conv,
        manufactured=unit,analytic_distribution=analytic,independent_line_cases=quad_cases,
        checks=checks,all_criteria_pass=all(checks.values()),seconds=time.monotonic()-START,peak_process_working_bytes=memory.memory_bytes(),
        limits='observation-only reconstruction; not a continuous-pressure or general structural work-preserving transfer')
    (FOLDER/'avl-linear-load-wrench-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else str(v))+'\n')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],peak_bytes=result['peak_process_working_bytes'])),flush=True)
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':
    try:main()
    except Exception as exc:
        (ROOT/'build/analysis/avl-linear-load-wrench-stopped.json').write_text(json.dumps(dict(error=repr(exc),completed=PARTIAL),indent=2)+'\n')
        raise
