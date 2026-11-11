"""Predeclared scalar, cancellation-aware refinement of the same P0 problem."""
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.integrate import quad_vec
from scipy.linalg import lu_factor,lu_solve
import check_uvlm_local_pressure_pitch as original

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()
ROWS=[]


def budget():
    if time.monotonic()-START>180 or original.P.memory.memory_bytes()>512*1024**2:
        raise RuntimeError('experiment_budget_exceeded')


def signed_h(s,x):
    r=np.asarray(s)-np.asarray(x)
    if np.any(r==0):raise ValueError('observation_on_vortex')
    return (np.sqrt(r*r+5000.**2)/r-1)/(2*np.pi*5000.)


def mesh_operators(m):
    bounds=(1-np.cos(np.pi*np.arange(m+1)/m))/2;dx=np.diff(bounds)
    xv=np.r_[bounds[:-1]+.25*dx,1.];xc=bounds[:-1]+.75*dx
    matrix=np.empty((m,m))
    for start in range(0,m,128):
        budget();x=xc[start:start+128,None]
        matrix[start:start+128]=signed_h(xv[:-1][None,:],x)-signed_h(xv[1:][None,:],x)
    fc=np.zeros(m);fc[-1]=1.2*8
    fu=1.2*np.diff(xv)
    mc0=8*fu;mc0[-1]-=1.2*8
    mu0=-.6*(xv[1:]**2-xv[:-1]**2)
    return dict(m=m,xv=xv,xc=xc,A=matrix,fc=fc,fu=fu,mc0=mc0,mu0=mu0)


def solve(op,k):
    budget();x=op['xc']
    def integrand(t):
        return original.D.density(1+t,x,5000.)*np.expm1(-2j*k*t)
    correction,error=quad_vec(integrand,0.,160.,epsrel=1e-11,epsabs=1e-12,limit=4096)
    matrix=op['A'].astype(complex)
    matrix[:,-1]=signed_h(op['xv'][-2],x)-signed_h(161.,x)+correction
    rhs=np.column_stack([1j*16*k*mode['H']-(8+1j*16*k*(x-mode['xp']))*mode['alpha'] for mode in original.MODES])
    factor=lu_factor(matrix);gamma=lu_solve(factor,rhs)
    rows=[];edges=[.05,.15,.3,.5,.7,.85,.95,1.]
    intervals=list(zip(edges[:-1],edges[1:]))+[(1-w,1.) for w in [.1,.05,.025,.0125]]
    for j,mode in enumerate(original.MODES):
        g=gamma[:,j];xp=mode['xp'];derivative=1j*16*k
        force=op['fc']@g+derivative*(op['fu']@g)
        moment=(op['mc0']+xp*op['fc'])@g+derivative*((op['mu0']+xp*op['fu'])@g)
        profile=[];ps=original.scales(mode)[0]
        for lower,upper in intervals:
            xv=op['xv'];overlap=np.maximum(0.,np.minimum(xv[1:],upper)-np.maximum(xv[:-1],lower))
            li=np.searchsorted(xv,lower,side='right')-1;ui=np.searchsorted(xv,upper,side='left')-1
            value=1.2*(derivative*(overlap@g)+8*(g[ui]-g[li]))/(upper-lower)
            ref=original.pressure_integral(k,mode,lower,upper)/(upper-lower)
            profile.append(dict(lower=lower,upper=upper,pressure=original.P.pair(value),reference=original.P.pair(ref),
                normalized_error=float(abs(value-ref)/max(ps,abs(ref)))))
        width=1/np.sqrt(op['m']);lower=1-width;xv=op['xv']
        overlap=np.maximum(0.,xv[1:]-np.maximum(xv[:-1],lower));li=np.searchsorted(xv,lower,side='right')-1
        pmean=1.2*(derivative*(overlap@g)+8*(g[-1]-g[li]))/width
        b,_=original.coefficients(k,mode);limit=(2/3)*1.2*8**2*(b[0]+4*b[1]+8*b[2])
        rows.append(dict(m=op['m'],k=k,mode=mode['name'],loads=[original.P.pair(force),original.P.pair(moment)],
            pressure_profile=profile,nonpenetration_relative=float(np.linalg.norm(matrix@g-rhs[:,j])/np.linalg.norm(rhs[:,j])),
            integral_error_estimate=float(error),trailing_width=float(width),
            trailing_divided_by_sqrt=original.P.pair(pmean/np.sqrt(width)),trailing_asymptotic=original.P.pair(limit),
            normalized_asymptotic_difference=float(abs(pmean/np.sqrt(width)-limit)/max(ps,abs(limit)))))
    return rows


def main():
    op=mesh_operators(128);_,(a,fc,fu,_),(mc0,mu0)=original.model(128)
    operator_errors=dict(A=float(np.linalg.norm(op['A']-a)/np.linalg.norm(a)),
        fc=float(np.linalg.norm(op['fc']-fc)),fu=float(np.linalg.norm(op['fu']-fu)),
        mc0=float(np.linalg.norm(op['mc0']-mc0)),mu0=float(np.linalg.norm(op['mu0']-mu0)))
    saved=json.loads((ROOT/'docs/verification/uvlm-local-pressure-pitch-result.json').read_text())
    baseline=[]
    for k in [.1,.5,1.]:
        rows=solve(op,k)
        for row,mode in zip(rows,original.MODES):
            old=next(r for r in saved['continuous_cases'] if r['m']==128 and r['k']==k and r['mode']==mode['name'])
            ps,fs,ms=original.scales(mode)
            load_error=max(abs(complex(*v)-complex(*r))/max(s,abs(complex(*r))) for v,r,s in zip(row['loads'],old['loads'],[fs,ms]))
            pressure_error=max(abs(complex(*v['pressure'])-complex(*r['pressure']))/max(ps,abs(complex(*r['reference']))) for v,r in zip(row['pressure_profile'],old['pressure_profile']))
            baseline.append(dict(k=k,mode=mode['name'],loads=float(load_error),pressure=float(pressure_error)))
    del op
    for m in [1024,2048]:
        op=mesh_operators(m)
        for k in [.1,.5,1.]:
            rows=solve(op,k);ROWS.extend(rows)
            print(json.dumps(dict(m=m,k=k,max_pressure_error=max(p['normalized_error'] for r in rows for p in r['pressure_profile']))),flush=True)
        del op
    changes=[]
    for k in [.1,.5,1.]:
        for mode in original.MODES:
            coarse=next(r for r in ROWS if r['m']==1024 and r['k']==k and r['mode']==mode['name'])
            fine=next(r for r in ROWS if r['m']==2048 and r['k']==k and r['mode']==mode['name'])
            ps=original.scales(mode)[0]
            value=max(abs(complex(*c['pressure'])-complex(*f['pressure']))/max(ps,abs(complex(*f['reference']))) for c,f in zip(coarse['pressure_profile'],fine['pressure_profile']))
            changes.append(dict(k=k,mode=mode['name'],maximum=float(value)))
    fine=[r for r in ROWS if r['m']==2048]
    checks=dict(operator_equivalence=max(operator_errors.values())<=1e-11,
        original_solution_equivalence=all(max(r['loads'],r['pressure'])<=1e-8 for r in baseline),
        local_pressure=all(p['normalized_error']<=.05 for r in fine for p in r['pressure_profile']),
        local_space_change=all(r['maximum']<=.02 for r in changes),
        nonpenetration=all(r['nonpenetration_relative']<=1e-11 for r in ROWS))
    source=Path(__file__);deps=['check_uvlm_local_pressure_pitch.py',*saved['dependency_sha256']]
    result=dict(protocol='docs/verification/uvlm-local-pressure-pitch-experiment.md#7',
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={name:hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in deps},
        baseline_result_sha256=hashlib.sha256((ROOT/'docs/verification/uvlm-local-pressure-pitch-result.json').read_bytes()).hexdigest(),
        U=8.,rho=1.2,c=1.,span=10000.,length=160.,modes=original.MODES,
        operator_errors=operator_errors,original_solution_comparison=baseline,continuous_cases=ROWS,
        pressure_refinement=changes,checks=checks,all_criteria_pass=bool(all(checks.values())),
        seconds=time.monotonic()-START,peak_process_working_bytes=original.P.memory.memory_bytes(),
        numpy=np.__version__,scipy=original.D.scipy.__version__,
        limits='same planar P0 pressure means with finer bound mesh; no P0 point pressure or 3D production adoption')
    (ROOT/'docs/verification/uvlm-local-pressure-refinement-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=original.P.core.json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,changes=changes,seconds=result['seconds']),default=original.P.core.json_scalar),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':
    try:main()
    except Exception as error:
        (ROOT/'build/analysis/uvlm-local-pressure-refinement-stopped.json').write_text(json.dumps(dict(
            error=repr(error),partial_cases=ROWS,seconds=time.monotonic()-START),indent=2,allow_nan=False,default=original.P.core.json_scalar)+'\n',encoding='utf-8')
        raise
