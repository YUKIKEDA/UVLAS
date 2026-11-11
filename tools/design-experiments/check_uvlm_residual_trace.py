"""Predeclared compensated residual/refinement and independent trace diagnosis."""
import hashlib
import itertools
import json
import math
import time
from decimal import Decimal,localcontext
from pathlib import Path
import numpy as np
from scipy.integrate import quad_vec
from scipy.linalg import lu_factor,lu_solve
from numpy.polynomial.legendre import leggauss
import check_uvlm_local_pressure_refinement as baseline

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()
O=baseline.original
ROWS=[]
SYSTEMS=[]


def budget():
    if time.monotonic()-START>180 or O.P.memory.memory_bytes()>512*1024**2:
        raise RuntimeError('experiment_budget_exceeded')


def two_product(a,b):
    for value in [a,b]:
        ab=np.abs(value)
        if np.any(~np.isfinite(ab)) or np.any(ab>2.**200) or np.any((ab!=0)&(ab<2.**-200)):
            raise ValueError('two_product_outside_tested_range')
    p=a*b;ca=(2.**27+1)*a;cb=(2.**27+1)*b
    ah=ca-(ca-a);al=a-ah;bh=cb-(cb-b);bl=b-bh
    e1=p-ah*bh;e2=e1-al*bh;e3=e2-ah*bl
    return p,al*bl-e3


def two_sum(a,b):
    s=a+b;bb=s-a;error=(a-(s-bb))+(b-bb)
    return s,error


def compensated(a,hi,rhs,lo=None,rounded=False):
    if np.any(a.imag[:,:-1]!=0):raise ValueError('unsupported_complex_pattern')
    out=np.empty_like(rhs);ordinary=np.empty_like(rhs) if rounded else None
    for start in range(0,len(a),16):
        budget();ar=a.real[start:start+16];ai=a.imag[start:start+16,-1]
        for j in range(hi.shape[1]):
            values=[]
            for vector in ([hi] if lo is None else [hi,lo]):
                pr,er=two_product(ar,vector[:,j].real)
                pi,ei=two_product(ar,vector[:,j].imag)
                qr,fr=two_product(ai,vector[-1,j].imag)
                qi,fi=two_product(ai,vector[-1,j].real)
                values.append((pr,er,pi,ei,qr,fr,qi,fi))
            for row in range(len(ar)):
                real=[[rhs[start+row,j].real]];imag=[[rhs[start+row,j].imag]]
                rawreal=real.copy();rawimag=imag.copy()
                for pr,er,pi,ei,qr,fr,qi,fi in values:
                    real.extend([-pr[row],-er[row],[qr[row],fr[row]]])
                    imag.extend([-pi[row],-ei[row],[-qi[row],-fi[row]]])
                    if rounded:
                        rawreal.extend([-pr[row],[qr[row]]]);rawimag.extend([-pi[row],[-qi[row]]])
                out[start+row,j]=complex(math.fsum(itertools.chain.from_iterable(real)),math.fsum(itertools.chain.from_iterable(imag)))
                if rounded:ordinary[start+row,j]=complex(math.fsum(itertools.chain.from_iterable(rawreal)),math.fsum(itertools.chain.from_iterable(rawimag)))
    return out,ordinary


def qualification():
    rng=np.random.default_rng(61030);a=np.ldexp(rng.uniform(-1,1,50),rng.integers(-20,21,50));b=np.ldexp(rng.uniform(-1,1,50),rng.integers(-20,21,50))
    hi,lo=two_product(a,b)
    with localcontext() as ctx:
        ctx.prec=100
        exact=all(Decimal.from_float(float(x))*Decimal.from_float(float(y))==Decimal.from_float(float(h))+Decimal.from_float(float(l)) for x,y,h,l in zip(a,b,hi,lo))
    matrix=np.array([[1+2.**-27,1]],complex);vector=np.array([[1-2.**-27],[-1]],complex)
    r,ordinary=compensated(matrix,vector,np.zeros((1,1),complex),rounded=True)
    ma=np.array([[1.e16,1.,-1.e16]],complex);xx=np.ones((3,1),complex)
    sr,_=compensated(ma,xx,np.zeros((1,1),complex))
    return dict(product_pairs=50,Decimal_product_identity=exact,
        product_cancellation=float(-r[0,0].real),product_expected=-2.**-54,
        rounded_products_result=float(-ordinary[0,0].real),sum_cancellation=float(-sr[0,0].real),sum_expected=1.)


def decimal_rows(a,hi,rhs,residual,lo=None):
    indexes=[0,1,len(a)//2,len(a)-2,len(a)-1];rows=[]
    with localcontext() as ctx:
        ctx.prec=80
        hr=[[Decimal.from_float(float(v)) for v in hi[:,j].real] for j in range(hi.shape[1])]
        im=[[Decimal.from_float(float(v)) for v in hi[:,j].imag] for j in range(hi.shape[1])]
        if lo is not None:
            for j in range(hi.shape[1]):
                hr[j]=[h+Decimal.from_float(float(v)) for h,v in zip(hr[j],lo[:,j].real)]
                im[j]=[h+Decimal.from_float(float(v)) for h,v in zip(im[j],lo[:,j].imag)]
        for i in indexes:
            ar=[Decimal.from_float(float(v)) for v in a[i].real];ai=Decimal.from_float(float(a[i,-1].imag))
            for j in range(hi.shape[1]):
                rr=Decimal.from_float(float(rhs[i,j].real))-sum((u*v for u,v in zip(ar,hr[j])),Decimal(0))+ai*im[j][-1]
                ri=Decimal.from_float(float(rhs[i,j].imag))-sum((u*v for u,v in zip(ar,im[j])),Decimal(0))-ai*hr[j][-1]
                value=complex(float(rr),float(ri));difference=abs(value-residual[i,j])/np.linalg.norm(rhs[:,j])
                rows.append(dict(row=i,column=j,decimal_residual=O.P.pair(value),normalized_difference=float(difference)))
    return rows


def evaluate(op,k,hi,lo=None):
    out=[];xv=op['xv'];edges=[.05,.15,.3,.5,.7,.85,.95,1.]
    intervals=list(zip(edges[:-1],edges[1:]))+[(1-w,1.) for w in [.1,.05,.025,.0125]]
    for j,mode in enumerate(O.MODES):
        vectors=[hi[:,j]] if lo is None else [hi[:,j],lo[:,j]]
        def add(values):return complex(math.fsum(v.real for v in values),math.fsum(v.imag for v in values))
        xp=mode['xp'];force=add([op['fc']@v+1j*16*k*(op['fu']@v) for v in vectors])
        moment=add([(op['mc0']+xp*op['fc'])@v+1j*16*k*((op['mu0']+xp*op['fu'])@v) for v in vectors])
        profile=[];ps=O.scales(mode)[0]
        for lower,upper in intervals:
            overlap=np.maximum(0.,np.minimum(xv[1:],upper)-np.maximum(xv[:-1],lower))
            li=np.searchsorted(xv,lower,side='right')-1;ui=np.searchsorted(xv,upper,side='left')-1
            integral=add([overlap@v for v in vectors]);left=add([v[li] for v in vectors]);right=add([v[ui] for v in vectors])
            value=1.2*(1j*16*k*integral+8*(right-left))/(upper-lower)
            ref=O.pressure_integral(k,mode,lower,upper)/(upper-lower)
            profile.append(dict(lower=lower,upper=upper,pressure=O.P.pair(value),reference=O.P.pair(ref),
                normalized_error=float(abs(value-ref)/max(ps,abs(ref))),integral=O.P.pair(integral),
                left_trace=O.P.pair(left),right_trace=O.P.pair(right)))
        out.append(dict(m=op['m'],k=k,mode=mode['name'],loads=[O.P.pair(force),O.P.pair(moment)],pressure_profile=profile))
    return out


def mu_reference(k,mode,points,n):
    q,w=leggauss(n);q=(q+1)/2;w=w/2;points=np.asarray(points)
    lower=np.arccos(2*points-1);theta=lower[:,None]+(np.pi-lower[:,None])*q[None,:]
    xx=(1+np.cos(theta))/2;b,_=O.coefficients(k,mode)
    p=b[0]*(1-np.cos(theta))+2*b[1]*np.sin(theta)**2+4*b[2]*np.cos(theta)*np.sin(theta)**2
    return 4*np.exp(-2j*k*points)*(np.pi-lower)*np.sum(w[None,:]*np.exp(2j*k*xx)*p,axis=1)


def trace_diagnosis(op,k,rows):
    result=[];validation=[];xv=op['xv']
    for mode,row in zip(O.MODES,rows):
        points=set()
        for p in row['pressure_profile']:
            for x,side in [(p['lower'],'right'),(p['upper'],'left')]:
                i=np.searchsorted(xv,x,side=side)-1
                points.update([x,(xv[i]+xv[i+1])/2])
        points=sorted(points);r128=mu_reference(k,mode,points,128);r256=mu_reference(k,mode,points,256)
        sensitivity=float(np.max(np.abs(r128-r256))/(8*(abs(mode['H'])+abs(mode['alpha']))))
        validation.append(dict(m=op['m'],k=k,mode=mode['name'],normalized_mu_sensitivity=sensitivity))
        refs=dict(zip(points,r256));ps=O.scales(mode)[0]
        for p in row['pressure_profile']:
            lower=p['lower'];upper=p['upper'];width=upper-lower
            li=np.searchsorted(xv,lower,side='right')-1;ui=np.searchsorted(xv,upper,side='left')-1
            pref=complex(*p['reference']);integral=complex(*p['integral'])
            iref=(pref*width/1.2-8*(refs[upper]-refs[lower]))/(1j*16*k)
            bulk=1.2*1j*16*k*(integral-iref)/width
            left=-1.2*8*(complex(*p['left_trace'])-refs[lower])/width
            right=1.2*8*(complex(*p['right_trace'])-refs[upper])/width
            total=complex(*p['pressure'])-pref;den=max(ps,abs(pref))
            center=(xv[li]+xv[li+1])/2
            position=-1.2*8*(refs[center]-refs[lower])/width
            result.append(dict(m=op['m'],k=k,mode=mode['name'],lower=lower,upper=upper,
                normalized_bulk=float(abs(bulk)/den),normalized_combined_trace=float(abs(left+right)/den),
                normalized_left=float(abs(left)/den),normalized_right=float(abs(right)/den),
                bulk=O.P.pair(bulk),left=O.P.pair(left),right=O.P.pair(right),
                closure=float(abs(bulk+left+right-total)/den),
                lower_cell_fraction=float((lower-xv[li])/(xv[li+1]-xv[li])),
                left_reference_position_term=O.P.pair(position),normalized_left_position=float(abs(position)/den)))
    return result,validation


def build(op,k):
    x=op['xc']
    correction,error=quad_vec(lambda t:O.D.density(1+t,x,5000.)*np.expm1(-2j*k*t),0.,160.,epsrel=1e-11,epsabs=1e-12,limit=4096)
    matrix=op['A'].astype(complex);matrix[:,-1]=baseline.signed_h(op['xv'][-2],x)-baseline.signed_h(161.,x)+correction
    rhs=np.column_stack([1j*16*k*mode['H']-(8+1j*16*k*(x-mode['xp']))*mode['alpha'] for mode in O.MODES])
    return matrix,rhs,float(error)


def hash_array(a):return hashlib.sha256(np.ascontiguousarray(a).view(np.uint8)).hexdigest()


def system(op,k):
    budget();a,b,error=build(op,k);hashes=[hash_array(a),hash_array(b)]
    factor=lu_factor(a);initial=lu_solve(factor,b);hi=initial.copy();lo=None
    comp,rounded=compensated(a,initial,b,rounded=True)
    norms=np.linalg.norm(b,axis=0)
    audits_initial=decimal_rows(a,initial,b,comp)
    history=[];low_history=[]
    def record(h,l,r,step):
        projected=h if l is None else h+l
        standard=np.array([np.linalg.norm(b[:,j]-a@projected[:,j])/norms[j] for j in range(3)])
        denominator=np.abs(b)+np.abs(a)@(np.abs(h) if l is None else np.abs(h)+np.abs(l))
        backward=np.max(np.divide(np.abs(r),denominator,out=np.zeros_like(denominator),where=denominator>0),axis=0)
        return dict(corrections=step,compensated_relative=(np.linalg.norm(r,axis=0)/norms).tolist(),
            BLAS_projected_relative=standard.tolist(),componentwise_backward=backward.tolist())
    for step in range(5):
        entry=record(hi,None,comp,step);history.append(entry)
        if max(entry['compensated_relative'])<=1e-11:break
        if step==4:break
        delta=lu_solve(factor,comp);updated=hi+delta
        if np.array_equal(updated,hi):break
        hi=updated;comp,_=compensated(a,hi,b)
    selected='float64'
    if max(history[-1]['compensated_relative'])>1e-11:
        hi=initial.copy();lo=np.zeros_like(initial);selected='high_low_float64'
        comp,_=compensated(a,hi,b,lo)
        for step in range(5):
            entry=record(hi,lo,comp,step);low_history.append(entry)
            if max(entry['compensated_relative'])<=1e-11:break
            if step==4:break
            delta=lu_solve(factor,comp)
            sr,er=two_sum(hi.real,delta.real);si,ei=two_sum(hi.imag,delta.imag)
            hr,lr=two_sum(sr,lo.real+er);hh,ll=two_sum(si,lo.imag+ei)
            hi=hr+1j*hh;lo=lr+1j*ll;comp,_=compensated(a,hi,b,lo)
    chosen=(low_history if lo is not None else history)[-1]
    audits_final=decimal_rows(a,hi,b,comp,lo)
    projected=hi if lo is None else hi+lo
    projection,_=compensated(a,projected,b)
    before=evaluate(op,k,initial);after=evaluate(op,k,hi,lo)
    changes=[]
    for mode,first,last in zip(O.MODES,before,after):
        ps,fs,ms=O.scales(mode)
        force=max(abs(complex(*u)-complex(*v))/max(s,abs(complex(*v))) for u,v,s in zip(first['loads'],last['loads'],[fs,ms]))
        pressure=max(abs(complex(*u['pressure'])-complex(*v['pressure']))/max(ps,abs(complex(*v['reference']))) for u,v in zip(first['pressure_profile'],last['pressure_profile']))
        changes.append(dict(mode=mode['name'],loads=float(force),pressure=float(pressure)))
    trace,mu_tests=trace_diagnosis(op,k,after)
    result=dict(m=op['m'],k=k,input_sha256=hashes,inputs_unchanged=hashes==[hash_array(a),hash_array(b)],
        initial_rounded_products_fsum_relative=(np.linalg.norm(rounded,axis=0)/norms).tolist(),
        float64_history=history,high_low_history=low_history,selected_representation=selected,
        selected_accuracy=chosen,float64_projection_compensated_relative=(np.linalg.norm(projection,axis=0)/norms).tolist(),
        Decimal_initial=audits_initial,Decimal_final=audits_final,observation_changes=changes,
        trace_diagnosis=trace,mu_reference_checks=mu_tests,continuous_integral_error_estimate=error)
    return result,after


def main():
    tests=qualification()
    for m in [1024,2048]:
        op=baseline.mesh_operators(m)
        for k in [.1,.5,1.]:
            report,rows=system(op,k);SYSTEMS.append(report);ROWS.extend(rows)
            print(json.dumps(dict(m=m,k=k,representation=report['selected_representation'],accuracy=report['selected_accuracy'],changes=report['observation_changes'])),flush=True)
        del op
    changes=[]
    for k in [.1,.5,1.]:
        for mode in O.MODES:
            coarse=next(r for r in ROWS if r['m']==1024 and r['k']==k and r['mode']==mode['name'])
            fine=next(r for r in ROWS if r['m']==2048 and r['k']==k and r['mode']==mode['name'])
            ps=O.scales(mode)[0]
            value=max(abs(complex(*a['pressure'])-complex(*b['pressure']))/max(ps,abs(complex(*b['reference']))) for a,b in zip(coarse['pressure_profile'],fine['pressure_profile']))
            changes.append(dict(k=k,mode=mode['name'],maximum=float(value)))
    fine=[r for r in ROWS if r['m']==2048]
    checks=dict(EFT_qualification=tests['Decimal_product_identity'] and tests['product_cancellation']==tests['product_expected'] and tests['sum_cancellation']==1.,
        Decimal_residual=all(r['normalized_difference']<=1e-18 for s in SYSTEMS for r in s['Decimal_initial']+s['Decimal_final']),
        compensated_residual=all(max(s['selected_accuracy']['compensated_relative'])<=1e-11 for s in SYSTEMS),
        backward_error=all(max(s['selected_accuracy']['componentwise_backward'])<=1e-12 for s in SYSTEMS),
        immutable_inputs=all(s['inputs_unchanged'] for s in SYSTEMS),
        small_observation_change=all(max(r['loads'],r['pressure'])<=1e-6 for s in SYSTEMS for r in s['observation_changes']),
        local_pressure=all(p['normalized_error']<=.05 for r in fine for p in r['pressure_profile']),
        local_space_change=all(r['maximum']<=.02 for r in changes),
        mu_reference=all(r['normalized_mu_sensitivity']<=1e-10 for s in SYSTEMS for r in s['mu_reference_checks']),
        error_decomposition=all(r['closure']<=1e-11 for s in SYSTEMS for r in s['trace_diagnosis']))
    source=Path(__file__);saved=json.loads((ROOT/'docs/verification/uvlm-local-pressure-refinement-result.json').read_text())
    deps=['check_uvlm_local_pressure_refinement.py',*saved['dependency_sha256']]
    result=dict(protocol='docs/verification/uvlm-residual-trace-experiment.md',
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={name:hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in deps},
        baseline_result_sha256=hashlib.sha256((ROOT/'docs/verification/uvlm-local-pressure-refinement-result.json').read_bytes()).hexdigest(),
        U=8.,rho=1.2,c=1.,span=10000.,length=160.,modes=O.MODES,
        float64_epsilon=float(np.finfo(np.float64).eps),longdouble_epsilon=float(np.finfo(np.longdouble).eps),
        qualification=tests,systems=SYSTEMS,continuous_cases=ROWS,pressure_refinement=changes,
        checks=checks,all_criteria_pass=bool(all(checks.values())),numpy=np.__version__,scipy=O.D.scipy.__version__,
        seconds=time.monotonic()-START,peak_process_working_bytes=O.P.memory.memory_bytes(),
        limits='represented float64 systems and tested compensated evaluator; not interval bounds, exact assembly, point pressure, or product type selection')
    (ROOT/'docs/verification/uvlm-residual-trace-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=O.P.core.json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,space_changes=changes,seconds=result['seconds']),default=O.P.core.json_scalar),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':
    try:main()
    except Exception as error:
        (ROOT/'build/analysis/uvlm-residual-trace-stopped.json').write_text(json.dumps(dict(error=repr(error),systems=SYSTEMS,continuous_cases=ROWS,seconds=time.monotonic()-START),indent=2,allow_nan=False,default=O.P.core.json_scalar)+'\n',encoding='utf-8')
        raise
