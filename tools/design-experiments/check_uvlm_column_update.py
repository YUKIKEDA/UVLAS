"""Predeclared real-LU/column-update pressure convergence; no model changes."""
import hashlib
import json
import math
import time
from decimal import Decimal,localcontext
from pathlib import Path
import numpy as np
from scipy.integrate import quad_vec
from scipy.linalg import lu_factor,lu_solve
import check_uvlm_residual_trace as previous

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()
O=previous.O
ROWS=[]
SYSTEMS=[]


def budget():
    if time.monotonic()-START>300 or O.P.memory.memory_bytes()>512*1024**2:
        raise RuntimeError('experiment_budget_exceeded')


# Only the resource callback is redirected to this predeclared 300 s budget.
previous.budget=budget


def real_operator(m):
    bounds=(1-np.cos(np.pi*np.arange(m+1)/m))/2;dx=np.diff(bounds)
    xv=np.r_[bounds[:-1]+.25*dx,1.];xc=bounds[:-1]+.75*dx
    a=np.empty((m,m))
    for start in range(0,m,64):
        budget();x=xc[start:start+64,None]
        a[start:start+64]=previous.baseline.signed_h(xv[:-1][None,:],x)-previous.baseline.signed_h(xv[1:][None,:],x)
    a[:,-1]=previous.baseline.signed_h(xv[-2],xc)-previous.baseline.signed_h(161.,xc)
    fc=np.zeros(m);fc[-1]=1.2*8;fu=1.2*np.diff(xv)
    mc0=8*fu;mc0[-1]-=1.2*8;mu0=-.6*(xv[1:]**2-xv[:-1]**2)
    return dict(m=m,xv=xv,xc=xc,A=a,fc=fc,fu=fu,mc0=mc0,mu0=mu0)


def block(op,v,start,end):
    a=op['A'][start:end].astype(complex);a[:,-1]+=v[start:end]
    return a


def effective_hash(op,v):
    h=hashlib.sha256()
    for start in range(0,op['m'],64):h.update(block(op,v,start,start+64).view(np.uint8))
    return h.hexdigest()


def base_solver(factor,b):
    return lu_solve(factor,np.asfortranarray(b.real),check_finite=False)+1j*lu_solve(factor,np.asfortranarray(b.imag),check_finite=False)


def update_solver(factor,v):
    q=base_solver(factor,v);den=1+q[-1]
    ratio=abs(den)/(1+abs(q[-1]))
    if not np.isfinite(ratio) or ratio<=1e-12:raise ValueError('unsafe_column_update_denominator')
    def solve(b):
        x=base_solver(factor,b)
        return x-q[:,None]*(x[-1][None,:]/den) if x.ndim==2 else x-q*x[-1]/den
    return solve,den,float(ratio)


def residual(op,v,hi,b,lo=None):
    result=np.empty_like(b);standard=np.empty_like(b);backward=np.zeros(3)
    for start in range(0,op['m'],64):
        budget();end=min(start+64,op['m']);a=block(op,v,start,end)
        r,_=previous.compensated(a,hi,b[start:end],lo)
        result[start:end]=r
        standard[start:end]=b[start:end]-a@hi
        if lo is not None:standard[start:end]-=a@lo
        absolute=np.abs(hi) if lo is None else np.abs(hi)+np.abs(lo)
        denom=np.abs(b[start:end])+np.abs(a)@absolute
        backward=np.maximum(backward,np.max(np.divide(np.abs(r),denom,out=np.zeros_like(denom),where=denom>0),axis=0))
    norms=np.linalg.norm(b,axis=0)
    return result,dict(compensated_relative=(np.linalg.norm(result,axis=0)/norms).tolist(),
        blocked_BLAS_relative=(np.linalg.norm(standard,axis=0)/norms).tolist(),componentwise_backward=backward.tolist())


def decimal_audit(op,v,hi,b,r,lo=None):
    rows=[]
    with localcontext() as ctx:
        ctx.prec=80
        hr=[[Decimal.from_float(float(x)) for x in hi[:,j].real] for j in range(3)]
        im=[[Decimal.from_float(float(x)) for x in hi[:,j].imag] for j in range(3)]
        if lo is not None:
            for j in range(3):
                hr[j]=[h+Decimal.from_float(float(x)) for h,x in zip(hr[j],lo[:,j].real)]
                im[j]=[h+Decimal.from_float(float(x)) for h,x in zip(im[j],lo[:,j].imag)]
        for i in [0,1,op['m']//2,op['m']-2,op['m']-1]:
            a=block(op,v,i,i+1)[0]
            ar=[Decimal.from_float(float(x)) for x in a.real];ai=Decimal.from_float(float(a[-1].imag))
            for j in range(3):
                real=Decimal.from_float(float(b[i,j].real))-sum((x*y for x,y in zip(ar,hr[j])),Decimal(0))+ai*im[j][-1]
                imag=Decimal.from_float(float(b[i,j].imag))-sum((x*y for x,y in zip(ar,im[j])),Decimal(0))-ai*hr[j][-1]
                value=complex(float(real),float(imag))
                rows.append(dict(row=i,column=j,normalized_difference=float(abs(value-r[i,j])/np.linalg.norm(b[:,j]))))
    return rows


def correction(op,v,b,solve):
    initial=solve(b);hi=initial.copy();lo=None;history=[];low_history=[]
    for step in range(5):
        r,stat=residual(op,v,hi,b);history.append(dict(corrections=step,**stat))
        print(json.dumps(dict(stage='refinement',m=op['m'],corrections=step,max_residual=max(stat['compensated_relative']))),flush=True)
        if max(stat['compensated_relative'])<=1e-11:break
        if step==4:break
        updated=hi+solve(r)
        if np.array_equal(updated,hi):break
        hi=updated
    kind='float64'
    if max(history[-1]['compensated_relative'])>1e-11:
        kind='high_low_float64';hi=initial.copy();lo=np.zeros_like(hi)
        for step in range(5):
            r,stat=residual(op,v,hi,b,lo);low_history.append(dict(corrections=step,**stat))
            print(json.dumps(dict(stage='high_low_refinement',m=op['m'],corrections=step,max_residual=max(stat['compensated_relative']))),flush=True)
            if max(stat['compensated_relative'])<=1e-11:break
            if step==4:break
            d=solve(r);sr,er=previous.two_sum(hi.real,d.real);si,ei=previous.two_sum(hi.imag,d.imag)
            hr,lr=previous.two_sum(sr,lo.real+er);hh,ll=previous.two_sum(si,lo.imag+ei)
            hi=hr+1j*hh;lo=lr+1j*ll
    accuracy=(low_history if lo is not None else history)[-1]
    audit=decimal_audit(op,v,hi,b,r,lo)
    if lo is None:projection=accuracy['compensated_relative']
    else:projection=residual(op,v,hi+lo,b)[1]['compensated_relative']
    return hi,lo,dict(selected_representation=kind,float64_history=history,high_low_history=low_history,
        selected_accuracy=accuracy,float64_projection_compensated_relative=projection,Decimal_rows=audit)


def qualification():
    rng=np.random.default_rng(613);a=rng.normal(size=(12,12));a+=20*np.eye(12)
    v=.1*(rng.normal(size=12)+1j*rng.normal(size=12));b=rng.normal(size=(12,3))+1j*rng.normal(size=(12,3))
    factor=lu_factor(a);solve,_,_=update_solver(factor,v)
    dense=a.astype(complex);dense[:,-1]+=v
    direct=np.linalg.solve(dense,b);small=float(np.linalg.norm(solve(b)-direct)/np.linalg.norm(direct))
    refused=False
    try:update_solver(lu_factor(np.eye(3)),np.array([0.,0.,-1.],complex))
    except ValueError:refused=True
    h,l=previous.two_sum(np.array([1.]),np.array([2.**-54]))
    with localcontext() as ctx:
        ctx.prec=100
        retained=Decimal.from_float(float(h[0]))+Decimal.from_float(float(l[0]))==Decimal(1)+Decimal(2)**-54
    op=real_operator(128);k=.5;x=op['xc']
    v,_=quad_vec(lambda t:O.D.density(1+t,x,5000.)*np.expm1(-2j*k*t),0.,160.,epsrel=1e-11,epsabs=1e-12,limit=4096)
    dense=block(op,v,0,128)
    b=np.column_stack([1j*16*k*m['H']-(8+1j*16*k*(x-m['xp']))*m['alpha'] for m in O.MODES])
    solve,_,_=update_solver(lu_factor(op['A']),v);g=solve(b);target=np.linalg.solve(dense,b)
    block_r,_=residual(op,v,g,b);full_r,_=previous.compensated(dense,g,b)
    first=previous.evaluate(op,k,g);last=previous.evaluate(op,k,target);change=0.
    for mode,aa,bb in zip(O.MODES,first,last):
        ps,fs,ms=O.scales(mode)
        change=max(change,max(abs(complex(*u)-complex(*w))/max(s,abs(complex(*w))) for u,w,s in zip(aa['loads'],bb['loads'],[fs,ms])),
            max(abs(complex(*u['pressure'])-complex(*w['pressure']))/max(ps,abs(complex(*w['reference']))) for u,w in zip(aa['pressure_profile'],bb['pressure_profile'])))
    return dict(synthetic_direct_difference=small,unsafe_denominator_refused=refused,
        high_low_addition_exact=retained,float64_projection_loses_low=bool(float(h[0]+l[0])==1.),
        aero_matrix_hash_equal=effective_hash(op,v)==previous.hash_array(dense),
        aero_direct_solution_difference=float(np.linalg.norm(g-target)/np.linalg.norm(target)),
        aero_observation_difference=float(change),blocked_compensated_difference=float(np.max(np.linalg.norm(block_r-full_r,axis=0)/np.linalg.norm(b,axis=0))))


def main():
    tests=qualification()
    old=json.loads((ROOT/'docs/verification/uvlm-residual-trace-result.json').read_text());comparisons=[]
    for m in [2048,4096]:
        op=real_operator(m);base_hash=previous.hash_array(op['A'])
        # Fortran copy is the only second square real array; LU overwrites that copy.
        factor=lu_factor(np.array(op['A'],order='F'),overwrite_a=True,check_finite=False)
        for k in [.1,.5,1.]:
            budget();x=op['xc']
            v,error=quad_vec(lambda t:O.D.density(1+t,x,5000.)*np.expm1(-2j*k*t),0.,160.,epsrel=1e-11,epsabs=1e-12,limit=4096)
            b=np.column_stack([1j*16*k*mode['H']-(8+1j*16*k*(x-mode['xp']))*mode['alpha'] for mode in O.MODES])
            hashes=[effective_hash(op,v),previous.hash_array(b)]
            solve,den,ratio=update_solver(factor,v)
            hi,lo,report=correction(op,v,b,solve)
            rows=previous.evaluate(op,k,hi,lo);ROWS.extend(rows)
            if m==2048:
                saved=next(s for s in old['systems'] if s['m']==m and s['k']==k)
                differences=[]
                for mode,row in zip(O.MODES,rows):
                    ref=next(r for r in old['continuous_cases'] if r['m']==m and r['k']==k and r['mode']==mode['name'])
                    ps,fs,ms=O.scales(mode)
                    difference=max(max(abs(complex(*u)-complex(*w))/max(s,abs(complex(*w))) for u,w,s in zip(row['loads'],ref['loads'],[fs,ms])),
                        max(abs(complex(*u['pressure'])-complex(*w['pressure']))/max(ps,abs(complex(*w['reference']))) for u,w in zip(row['pressure_profile'],ref['pressure_profile'])))
                    differences.append(float(difference))
                comparisons.append(dict(k=k,input_hashes_equal=hashes==saved['input_sha256'],observation_differences=differences))
            report.update(m=m,k=k,input_sha256=hashes,base_sha256=base_hash,
                inputs_unchanged=hashes==[effective_hash(op,v),previous.hash_array(b)] and base_hash==previous.hash_array(op['A']),
                denominator=O.P.pair(den),relative_denominator=ratio,integral_error_estimate=float(error))
            SYSTEMS.append(report)
            print(json.dumps(dict(stage='system_done',m=m,k=k,representation=report['selected_representation'],accuracy=report['selected_accuracy'])),flush=True)
        del solve,factor,op,hi,lo
    changes=[]
    for k in [.1,.5,1.]:
        for mode in O.MODES:
            coarse=next(r for r in ROWS if r['m']==2048 and r['k']==k and r['mode']==mode['name'])
            fine=next(r for r in ROWS if r['m']==4096 and r['k']==k and r['mode']==mode['name']);ps=O.scales(mode)[0]
            differences=[abs(complex(*a['pressure'])-complex(*b['pressure']))/max(ps,abs(complex(*b['reference']))) for a,b in zip(coarse['pressure_profile'],fine['pressure_profile'])]
            index=int(np.argmax(differences));p=fine['pressure_profile'][index]
            changes.append(dict(k=k,mode=mode['name'],maximum=float(differences[index]),lower=p['lower'],upper=p['upper']))
    fine=[r for r in ROWS if r['m']==4096]
    checks=dict(synthetic_update=tests['synthetic_direct_difference']<=1e-12,
        unsafe_denominator=tests['unsafe_denominator_refused'],high_low_addition=tests['high_low_addition_exact'] and tests['float64_projection_loses_low'],
        aero_update_equivalence=tests['aero_matrix_hash_equal'] and tests['aero_direct_solution_difference']<=1e-10 and tests['aero_observation_difference']<=1e-8,
        blocked_residual_equivalence=tests['blocked_compensated_difference']<=1e-18,
        saved_input_equivalence=all(r['input_hashes_equal'] for r in comparisons),
        saved_observation_equivalence=all(max(r['observation_differences'])<=1e-8 for r in comparisons),
        immutable_inputs=all(s['inputs_unchanged'] for s in SYSTEMS),
        compensated_residual=all(max(s['selected_accuracy']['compensated_relative'])<=1e-11 for s in SYSTEMS),
        backward_error=all(max(s['selected_accuracy']['componentwise_backward'])<=1e-12 for s in SYSTEMS),
        Decimal_residual=all(r['normalized_difference']<=1e-18 for s in SYSTEMS for r in s['Decimal_rows']),
        local_pressure=all(p['normalized_error']<=.05 for r in fine for p in r['pressure_profile']),
        local_space_change=all(r['maximum']<=.02 for r in changes))
    source=Path(__file__);deps=['check_uvlm_residual_trace.py',*old['dependency_sha256']]
    result=dict(protocol='docs/verification/uvlm-column-update-convergence.md',
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={name:hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in deps},
        baseline_result_sha256=hashlib.sha256((ROOT/'docs/verification/uvlm-residual-trace-result.json').read_bytes()).hexdigest(),
        U=8.,rho=1.2,c=1.,span=10000.,length=160.,modes=O.MODES,
        qualification=tests,saved_comparisons=comparisons,systems=SYSTEMS,continuous_cases=ROWS,pressure_refinement=changes,
        checks=checks,all_criteria_pass=bool(all(checks.values())),numpy=np.__version__,scipy=O.D.scipy.__version__,
        seconds=time.monotonic()-START,peak_process_working_bytes=O.P.memory.memory_bytes(),
        resource_callback_override='previous.budget=budget; numerical helper unchanged; declared 300s/512MiB',
        limits='same P0 local means and finite planar diagnostic; not point pressure, general 3D, solver or product precision adoption')
    (ROOT/'docs/verification/uvlm-column-update-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=O.P.core.json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,space_changes=changes,seconds=result['seconds']),default=O.P.core.json_scalar),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':
    try:main()
    except Exception as error:
        (ROOT/'build/analysis/uvlm-column-update-stopped.json').write_text(json.dumps(dict(error=repr(error),systems=SYSTEMS,continuous_cases=ROWS,seconds=time.monotonic()-START),indent=2,allow_nan=False,default=O.P.core.json_scalar)+'\n',encoding='utf-8')
        raise
