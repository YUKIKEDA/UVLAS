"""Predeclared TE/wake/differentiation diagnosis; not production numerics."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad_vec
from scipy.linalg import lu_factor,lu_solve
import check_uvlm_harmonic_heave as previous

ROOT=rd_paths.ROOT
START=time.monotonic()
ROWS=[]
BASE={}


def budget():
    if time.monotonic()-START>180 or previous.memory.memory_bytes()>512*1024**2:
        raise RuntimeError('experiment_budget_exceeded')


def horseshoe(s,x,b):
    r=np.asarray(s)-np.asarray(x)
    return b/(2*np.pi*r*(np.sqrt(r*r+b*b)+r))


def density(s,x,b):
    r=np.asarray(s)-np.asarray(x)
    return b/(2*np.pi*r*r*np.sqrt(r*r+b*b))


def weights(s,d,x,b):
    r=np.asarray(s)-np.asarray(x);root=np.sqrt(r*r+b*b)
    endroot=np.sqrt((r+d)**2+b*b)
    diff=d*(2*r+d)/(root+endroot)
    average=((diff-d)/b+np.log1p(d/r)-np.log1p(diff/(b+root)))/(2*np.pi*d)
    left=horseshoe(s,x,b)-average
    right=average-horseshoe(s+d,x,b)
    return left,right


def algebraic_tests():
    q,w=leggauss(64);q=(q+1)/2;w=w/2;rows=[]
    for x in [.25,.9,.998046875]:
        for s in [1.,1.01,2.,81.]:
            for d in [.00078125,.025,.1]:
                p=np.array([[s,-5000.,0.],[s,5000.,0.],[s+d,5000.,0.],[s+d,-5000.,0.]])
                direct=previous.core.kernel([[x,0.,0.]],p,np.roll(p,-1,axis=0),0.)[0,:,2].sum()
                expected=horseshoe(s,x,5000.)-horseshoe(s+d,x,5000.)
                a,b=weights(s,d,x,5000.)
                f=density(s+d*q,x,5000.)
                aq=d*np.dot(w,f*(1-q));bq=d*np.dot(w,f*q)
                eligible=lambda error,ref:bool(abs(error)<=1e-11 or abs(error)/max(abs(ref),1e-300)<=1e-10)
                rows.append(dict(x=x,s=s,d=d,ring_error=float(direct-expected),
                    left_error=float(a-aq),right_error=float(b-bq),constant_error=float(a+b-expected),
                    eligible=eligible(direct-expected,expected) and eligible(a-aq,aq) and eligible(b-bq,bq)
                    and eligible(a+b-expected,expected)))
    return rows


def base(m):
    if m not in BASE:
        grid=previous.geometry.make_grid(m,1,10000.,wake=False)
        a=previous.geometry.aic(grid,0.)
        ends=np.array([[grid['points'][i],grid['points'][j]] for i,j in grid['edges']])
        dl=ends[:,1]-ends[:,0];ids=grid['body_ids']
        fc=1.2*8*(dl[ids,1]@grid['body'][ids])/10000.
        p=grid['points'][grid['triangles']]
        area=.5*np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])[:,2]
        fu=np.bincount(grid['tri_owner'],weights=1.2*area,minlength=m)/10000.
        BASE[m]=(a,fc,fu,grid['collocation'][:,0])
    return BASE[m]


def continuous(m,k,length,tight=False):
    budget();x=base(m)[3];kappa=2*k
    def f(s):return density(1+s,x,5000.)*np.exp(-1j*kappa*s)
    value,error=quad_vec(f,0.,length,epsrel=1e-12 if tight else 1e-10,
        epsabs=1e-13 if tight else 1e-11,limit=4096)
    return value,float(error)


def discrete(m,k,delta,length,keep_nodes=False):
    budget();x=base(m)[3][:,None];n=int(np.ceil(length/delta))
    old=np.zeros(m,complex);p1=old.copy();wrong=old.copy()
    nodes=np.zeros((m,n+1)) if keep_nodes else None
    z=np.exp(2j*k*delta)
    for start in range(0,n,512):
        budget();j=np.arange(start,min(start+512,n));s=1.+delta*j[None,:]
        a,b=weights(s,delta,x,5000.)
        cell=horseshoe(s,x,5000.)-horseshoe(s+delta,x,5000.)
        phase=z**(-j)
        old+=cell@(phase/z);p1+=a@phase+b@(phase/z);wrong+=cell@phase
        if keep_nodes:
            nodes[:,j]+=a;nodes[:,j+1]+=b
    left0=weights(1.,delta,x[:,0],5000.)[0]
    return dict(old=old,p1=p1,unconverted=wrong,left0=left0,nodes=nodes,n=n,effective_length=n*delta)


def solve(m,k,kw,delta=None,method='continuous',derivative='exact',length=160.):
    a,fc,fu,_=base(m);matrix=a.astype(complex);matrix[:,-1]+=kw
    rhs=np.full(m,1j*(16*k)*.001)
    gamma=np.linalg.solve(matrix,rhs)
    d=1j*16*k if derivative=='exact' else (1-np.exp(-2j*k*delta))/(delta/8)
    circ=fc@gamma;noncirc=d*(fu@gamma);force=circ+noncirc
    ref=previous.reference(k)[0]
    row=dict(m=m,k=k,delta=delta,length=length,method=method,derivative=derivative,
        force=previous.pair(force),circulating=previous.pair(circ),unsteady=previous.pair(noncirc),
        reference=previous.pair(ref),complex_relative=float(abs(force/ref-1)),
        amplitude_relative=float(abs(abs(force/ref)-1)),phase_degrees=float(np.angle(force/ref,deg=True)),
        nonpenetration_relative=float(np.linalg.norm(matrix@gamma-rhs)/np.linalg.norm(rhs)))
    return row,gamma,force


def march():
    m=32;k=.5;delta=np.pi/(k*512);steps=1024
    operators=discrete(m,k,delta,40.,keep_nodes=True)
    row,gamma,force=solve(m,k,operators['p1'],delta,'W-P1','BDF1',operators['effective_length'])
    a,fc,fu,_=base(m);matrix=a.copy();matrix[:,-1]+=operators['left0']
    factor=lu_factor(matrix);g=operators['nodes'];z=np.exp(2j*k*delta)
    accepted=(z**(-np.arange(operators['n']+1))*gamma[-1]).real
    old=gamma.real.copy();maximum=0.;residual=0.;trace=[]
    for n in range(1,steps+1):
        if n%128==0:budget()
        trial=np.empty_like(accepted);trial[1:]=accepted[:-1]
        velocity=float(np.real(1j*16*k*.001*z**n))
        rhs=np.full(m,velocity)-g[:,1:]@trial[1:]
        current=lu_solve(factor,rhs);trial[0]=current[-1]
        load=float(fc@current+fu@((current-old)/(delta/8)))
        target=float(np.real(force*z**n))
        maximum=max(maximum,abs(load-target)/abs(force))
        full_residual=a@current+g@trial-velocity
        residual=max(residual,float(np.linalg.norm(full_residual)/max(np.linalg.norm(rhs),1e-12)))
        if n%64==0:trace.append(dict(step=n,force=load,target=target))
        accepted=trial;old=current
    return dict(m=m,k=k,delta=delta,steps=steps,steps_per_period=512,length=40.,
        effective_length=operators['effective_length'],nwake=operators['n'],initialization='periodic nodal history',
        max_periodic_difference=maximum,max_nonpenetration_relative=residual,trace=trace,frequency_case=row)


def main():
    algebra=algebraic_tests();sensitivity=[];convergence=[];negatives=[];baseline=[];component=[]
    saved=json.loads((ROOT/'rd/aero/unsteady-reference/records/uvlm-harmonic-heave-result.json').read_text())['harmonic_cases']
    for k in [.1,.5,1.]:
        cont={}
        for m in [16,32,64,128]:
            kw,error=continuous(m,k,160.);cont[m]=kw
            row,_,_=solve(m,k,kw);row['integration_error_estimate']=error;ROWS.append(row)
            print(json.dumps(row),flush=True)
        tighter,error=continuous(128,k,160.,tight=True)
        row,_,_=solve(128,k,tighter);normal=next(r for r in ROWS if r['k']==k and r['m']==128 and r['method']=='continuous')
        sensitivity.append(dict(k=k,force_change=previous.relative(row,normal),tight_error_estimate=error))
        finest={}
        for m in [32,64,128]:
            for delta in [.05,.025,.0125]+([.00625] if m==128 else []):
                ops=discrete(m,k,delta,160.)
                current={}
                for method in ['old','p1']:
                    for differentiation in ['BDF1','exact']:
                        row,_,_=solve(m,k,ops[method],delta,method,differentiation)
                        ROWS.append(row);current[(method,differentiation)]=row
                cdf,_,_=solve(m,k,cont[m],delta,'continuous','BDF1')
                ROWS.append(cdf)
                component.append(dict(m=m,k=k,delta=delta,
                    old_to_continuous_exact=previous.relative(current[('old','exact')],
                        next(r for r in ROWS if r['k']==k and r['m']==m and r['method']=='continuous' and r['derivative']=='exact')),
                    p1_to_continuous_exact=previous.relative(current[('p1','exact')],
                        next(r for r in ROWS if r['k']==k and r['m']==m and r['method']=='continuous' and r['derivative']=='exact')),
                    BDF1_to_exact_old=previous.relative(current[('old','BDF1')],current[('old','exact')]),
                    BDF1_to_exact_p1=previous.relative(current[('p1','BDF1')],current[('p1','exact')])))
                for old in saved:
                    if old['m']==m and old['k']==k and old['delta']==delta and old['length']==160. and old['span']==10000. and old['eps']==0:
                        baseline.append(dict(m=m,k=k,delta=delta,relative=previous.relative(current[('old','BDF1')],old)))
                if m==128 and delta in [.0125,.00625]:
                    finest[delta]=float(np.linalg.norm(ops['p1']-cont[m])/np.linalg.norm(cont[m]))
                if m==128 and delta==.00625:
                    correct=current[('p1','BDF1')]
                    omitted,_,_=solve(m,k,ops['p1']-ops['left0'],delta,'omit-left0','BDF1')
                    unconverted,_,_=solve(m,k,ops['unconverted'],delta,'unconverted','BDF1')
                    negatives.append(dict(k=k,omitted=previous.relative(omitted,correct),
                        unconverted=previous.relative(unconverted,correct)))
            print(json.dumps(dict(stage='discrete_done',k=k,m=m)),flush=True)
        convergence.append(dict(k=k,coarse=finest[.0125],fine=finest[.00625]))
    time_result=march()
    contfine=[r for r in ROWS if r['method']=='continuous' and r['m']==128 and r['derivative']=='exact']
    p1fine=[r for r in ROWS if r['method']=='p1' and r['m']==128 and r['delta']==.00625 and r['derivative']=='BDF1']
    checks=dict(analytic_coefficients=all(r['eligible'] for r in algebra),
        continuous_integral_sensitivity=all(r['force_change']<=1e-8 for r in sensitivity),
        saved_B1_match=len(baseline)==3 and all(r['relative']<=1e-9 for r in baseline),
        p1_operator_convergence=all(r['fine']<r['coarse'] for r in convergence),
        continuous_complex=all(r['complex_relative']<=.02 for r in contfine),
        continuous_amplitude=all(r['amplitude_relative']<=.02 for r in contfine),
        continuous_phase=all(abs(r['phase_degrees'])<=1 for r in contfine),
        p1_complex=all(r['complex_relative']<=.02 for r in p1fine),
        p1_amplitude=all(r['amplitude_relative']<=.02 for r in p1fine),
        p1_phase=all(abs(r['phase_degrees'])<=1 for r in p1fine),
        nonpenetration=all(r['nonpenetration_relative']<=1e-11 for r in ROWS),
        periodic_time_match=time_result['max_periodic_difference']<=1e-10,
        periodic_nonpenetration=time_result['max_nonpenetration_relative']<=1e-11,
        negative_omitted=max(r['omitted'] for r in negatives)>.001,
        negative_unconverted=max(r['unconverted'] for r in negatives)>.001)
    source=Path(__file__)
    deps=['check_uvlm_harmonic_heave.py',*json.loads((ROOT/'rd/aero/unsteady-reference/records/uvlm-harmonic-heave-result.json').read_text())['dependency_sha256']]
    result=dict(protocol='rd/aero/wake-loads/records/uvlm-trailing-edge-diagnostic.md',geometry_id=previous.geometry.GEOMETRY,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={name:hashlib.sha256(rd_paths.sibling(name).read_bytes()).hexdigest() for name in deps},
        U=8.,rho=1.2,H=.001,c=1.,span=10000.,core=0.,numpy=np.__version__,scipy=scipy.__version__,
        algebraic_cases=algebra,harmonic_cases=ROWS,integral_sensitivity=sensitivity,
        baseline_saved_comparison=baseline,p1_operator_convergence=convergence,component_differences=component,
        negative_controls=negatives,time_march=time_result,checks=checks,all_criteria_pass=bool(all(checks.values())),
        seconds=time.monotonic()-START,peak_process_working_bytes=previous.memory.memory_bytes(),
        limits='planar centerline one-span scalar diagnostic; no direct TE pressure/Kutta verification or 3D production adoption')
    (ROOT/'rd/aero/wake-loads/records/uvlm-trailing-edge-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=previous.core.json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],sensitivity=sensitivity,convergence=convergence,
        baseline=baseline,time_error=time_result['max_periodic_difference']),default=previous.core.json_scalar),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':
    try:main()
    except Exception as error:
        (ROOT/'build/analysis/uvlm-trailing-edge-stopped.json').write_text(json.dumps(dict(
            error=repr(error),partial_harmonic_cases=ROWS,seconds=time.monotonic()-START),indent=2,allow_nan=False)+'\n',encoding='utf-8')
        raise
