"""Prewritten protocol: docs/verification/q9-transient-accuracy-diagnostic.md.
All 69 coordinates retained; structural time diagnostic, not product code.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.integrate import solve_ivp
import check_q9_nonlinear_restart as previous
import check_q9_material_tangent as material
import check_q9_plate_convergence as plate

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()


def budget():
    if time.monotonic()-START>240 or plate.memory_bytes()>4*1024**3:
        raise TimeoutError('240 second / 4 GiB experiment limit')


class Model:
    def __init__(self):
        f=previous.Fixture();self.f=f
        _,k,_=f.evaluate(np.zeros(69))
        values,self.v=eigh(k,f.m)
        self.omega=np.sqrt(values);self.w1=self.omega[0]
        self.basis=f.p@self.v
        self.d=np.einsum('iak,piJ->paJk',self.basis.reshape(27,3,69),f.ctx['g'])
        self.fstar=np.einsum('ia,piJ->paJ',f.qstar.reshape(27,3),f.ctx['g'])
        self.e0=f.ctx['e0'];self.c=f.ctx['c'];self.w=f.ctx['w']
        e=.5*(np.einsum('paI,paJ->pIJ',self.fstar,self.fstar)-np.eye(3))-self.e0
        self.sstar=material.stress_tensor(material.engineering(e)@self.c.T)
        self.rstar=self.v.T@f.evaluate(np.zeros(69))[0]
        self.y0=self.v.T@f.m@f.x0
        self.amplitude=np.linalg.norm(self.y0)
        self.initial_speed_scale=self.w1*self.amplitude

    def evaluate(self,y,tangent=True):
        h=np.einsum('paJk,k->paJ',self.d,y)
        f=self.fstar+h
        de=.5*(np.einsum('paI,paJ->pIJ',self.fstar,h)+
            np.einsum('paI,paJ->pIJ',h,self.fstar)+np.einsum('paI,paJ->pIJ',h,h))
        ds=material.stress_tensor(material.engineering(de)@self.c.T)
        dp=np.einsum('paI,pIJ->paJ',h,self.sstar)+np.einsum('paI,pIJ->paJ',self.fstar,ds)+np.einsum('paI,pIJ->paJ',h,ds)
        force=self.rstar+np.einsum('p,paJk,paJ->k',self.w,self.d,dp)
        if not tangent:return force
        b=np.zeros((len(self.w),6,69),dtype=np.result_type(y))
        for row,(i,j) in enumerate(material.PAIRS):
            b[:,row,:]=np.einsum('pa,pak->pk',f[:,:,i],self.d[:,:,j,:])
            if i!=j:b[:,row,:]+=np.einsum('pa,pak->pk',f[:,:,j],self.d[:,:,i,:])
        kmat=np.einsum('p,pik,ij,pjl->kl',self.w,b,self.c,b,optimize=True)
        kgeo=np.einsum('p,paIk,pIJ,paJl->kl',self.w,self.d,self.sstar+ds,self.d,optimize=True)
        return force,kmat+kgeo,float(np.linalg.det(f).min())

    def physical(self,y):
        return self.f.s@(self.basis@y).reshape(27,3)


def equivalence(model):
    rng=np.random.default_rng(20261030);rows=[]
    for case in range(6):
        y=model.y0.copy()
        if case:
            perturb=rng.standard_normal(69)/model.omega
            perturb*=2e-6/np.max(np.abs(model.physical(perturb)))
            y+=perturb
        r,k,det=model.evaluate(y)
        oldr,oldk,_=model.f.evaluate(model.v@y)
        # Natural oscillator energy coordinates prevent mixed-unit matrix norms.
        invw=1/model.omega
        force_error=np.linalg.norm(invw*(r-model.v.T@oldr))/max(np.linalg.norm(invw*r),model.initial_speed_scale)
        transformed=model.v.T@oldk@model.v
        tangent_error=np.linalg.norm(invw[:,None]*(k-transformed)*invw[None,:])/np.linalg.norm(invw[:,None]*k*invw[None,:])
        direction=rng.standard_normal(69)/model.omega
        direction/=np.linalg.norm(direction)
        cs=np.imag(model.evaluate(y.astype(complex)+1j*1e-25*direction,False))/1e-25
        derivative_error=np.linalg.norm(invw*(cs-k@direction))/np.linalg.norm(invw*cs)
        rows.append(dict(case=case,force_scaled_relative=float(force_error),tangent_scaled_relative=float(tangent_error),complex_step_relative=float(derivative_error),sampled_det_f_min=det))
    mass_error=float(np.linalg.norm(model.v.T@model.f.m@model.v-np.eye(69)))
    return dict(cases=rows,mass_orthogonality=mass_error,eligible=bool(mass_error<=1e-9 and all(
        r['force_scaled_relative']<=1e-7 and r['tangent_scaled_relative']<=1e-7 and r['complex_step_relative']<=1e-9 for r in rows)))


def midpoint(model,duration,n,samples,tolerance=1e-11):
    assert n%samples==0
    dt=duration/n;y=model.y0.copy();v=np.zeros(69)
    ys=[y.copy()];vs=[v.copy()];maxres=0.;mindet=1e30;maxiteration=0
    for step in range(n):
        budget();z=v.copy()
        for iteration in range(13):
            r,k,det=model.evaluate(y+.5*dt*z)
            residual=z-v+.5*dt*r
            correction=np.linalg.solve(np.eye(69)+.25*dt**2*k,-residual)
            relative=float(np.linalg.norm(residual)/model.initial_speed_scale)
            physical=float(np.max(np.abs(model.physical(dt*correction))))
            if relative<=tolerance and physical<=1e-12:break
            if iteration==12:raise ValueError(dict(reason='newton_limit',n=n,step=step,relative=relative,tolerance=tolerance))
            z+=correction
        y=y+dt*z;v=2*z-v
        # Check endpoints too, not only midpoint geometry.
        _,_,enddet=model.evaluate(y)
        mindet=min(mindet,det,enddet);maxres=max(maxres,relative);maxiteration=max(maxiteration,iteration)
        if (step+1)%(n//samples)==0:ys.append(y.copy());vs.append(v.copy())
    y=np.array(ys);v=np.array(vs)
    modal_energy=.5*(v*v+y*y*model.omega**2)
    initial=.5*np.sum(model.y0**2*model.omega**2)
    return dict(y=y,v=v,metadata=dict(n=n,tolerance=tolerance,residual_max=maxres,sampled_det_f_min=mindet,max_iterations=maxiteration,
        modal_energy_max_fraction=(np.max(modal_energy,axis=0)/initial).tolist(),
        high65_linear_energy_fraction_max=float(np.max(np.sum(modal_energy[:,4:],axis=1))/initial)))


def errors(a,b,model):
    dy=(a['y']-b['y'])/model.amplitude
    dv=(a['v']-b['v'])/model.initial_speed_scale
    return dict(displacement_mass_max=float(np.linalg.norm(dy,axis=1).max()),
        velocity_scaled_mass_max=float(np.linalg.norm(dv,axis=1).max()),
        phase_space_max=float(np.sqrt(np.sum(dy*dy+dv*dv,axis=1)).max()),
        low4_displacement_max=float(np.linalg.norm(dy[:,:4],axis=1).max()),
        high65_displacement_max=float(np.linalg.norm(dy[:,4:],axis=1).max()),
        low4_velocity_max=float(np.linalg.norm(dv[:,:4],axis=1).max()),
        high65_velocity_max=float(np.linalg.norm(dv[:,4:],axis=1).max()))


def oracle(model,duration,rtol,atol):
    # Dimensionless time/coordinates give both tolerances a stated meaning.
    def rhs(t,z):
        budget()
        return np.r_[z[69:],-model.evaluate(model.amplitude*z[:69],False)/(model.w1**2*model.amplitude)]
    z0=np.r_[model.y0/model.amplitude,np.zeros(69)]
    ts=np.linspace(0.,model.w1*duration,161)
    start=time.monotonic()
    solution=solve_ivp(rhs,(0.,ts[-1]),z0,method='DOP853',rtol=rtol,atol=atol,t_eval=ts)
    if not solution.success:raise ValueError(solution.message)
    return dict(y=solution.y[:69].T*model.amplitude,v=solution.y[69:].T*model.initial_speed_scale,
        metadata=dict(method='DOP853',rtol=rtol,atol=atol,nfev=solution.nfev,seconds=time.monotonic()-start))


def original_path(model,duration):
    state=previous.initial_state(model.f);ys=[model.y0.copy()];vs=[np.zeros(69)]
    for step in range(320):
        trial=previous.try_step(model.f,state,duration/320)
        state=previous.commit(state,trial)
        if (step+1)%4==0:
            ys.append(model.v.T@model.f.m@np.array(state['x']))
            vs.append(model.v.T@model.f.m@np.array(state['v']))
    return dict(y=np.array(ys),v=np.array(vs))


def main():
    refined='--refine' in sys.argv
    model=Model();same=equivalence(model)
    print(json.dumps(dict(stage='equivalence',eligible=same['eligible'])),flush=True)
    if not same['eligible']:raise ValueError('incremental formulation equivalence failed')
    period=2*np.pi/model.w1;long=[];short=[]
    for n in [80,160,320,640]:
        case=midpoint(model,period/2,n,80);long.append(case)
        print(json.dumps(dict(stage='half_period',**{k:v for k,v in case['metadata'].items() if k!='modal_energy_max_fraction'})),flush=True)
    strict=midpoint(model,period/2,320,80,1e-13)
    longdiff=[errors(a,b,model) for a,b in zip(long[:-1],long[1:])]
    iteration_difference=errors(long[2],strict,model)
    original_difference=errors(long[2],original_path(model,period/2),model) if refined else None
    duration=period/200
    ref=oracle(model,duration,1e-10,1e-12)
    ref2=oracle(model,duration,1e-12,1e-14)
    sensitivity=errors(ref,ref2,model)
    for n in ([160,320,640,1280,2560] if refined else [160,320,640,1280]):
        case=midpoint(model,duration,n,160);short.append(case)
        print(json.dumps(dict(stage='short_window',n=n,error=errors(case,ref2,model)['phase_space_max'])),flush=True)
    strict_short=midpoint(model,duration,640,160,1e-13)
    short_errors=[errors(case,ref2,model) for case in short]
    short_iteration=errors(short[2],strict_short,model)
    order=float(np.log2(short_errors[-2]['phase_space_max']/short_errors[-1]['phase_space_max']))
    nonlinear=model.evaluate(model.y0,False)-model.omega**2*model.y0
    correction=-nonlinear/model.omega**2
    nonlinear_estimate=dict(modal_acceleration=nonlinear.tolist(),linear_static_compatibility_correction=correction.tolist(),
        physical_correction_max_m=float(np.max(np.abs(model.physical(correction)))),
        correction_mass_relative=float(np.linalg.norm(correction)/model.amplitude))
    phase=[]
    for n in [80,160,320,640,2560]:
        dt=(period/2)/n
        lag=(period/2)*(model.omega-2*np.arctan(.5*dt*model.omega)/dt)
        phase.append(dict(half_period_steps=n,highest_h_omega=float(dt*model.omega[-1]),
            lag_lowest_rad=float(lag[0]),lag_highest_rad=float(lag[-1]),all_lags_rad=lag.tolist()))
    checks=dict(formulation_equivalence=same['eligible'],reference_sensitivity=sensitivity['phase_space_max']<=1e-7,
        short_window_error=short_errors[-1]['phase_space_max']<=1e-3,short_window_order=1.7<=order<=2.3,
        short_iteration_budget=short_iteration['phase_space_max']<=.1*short_errors[2]['phase_space_max'],
        half_period_iteration_relative=iteration_difference['phase_space_max']<=.1*longdiff[-1]['phase_space_max'])
    if refined:checks['original_formulation_difference']=original_difference['phase_space_max']<=.1*longdiff[-1]['phase_space_max']
    output=ROOT/'build/analysis/q9-transient-accuracy-trajectories.npz'
    np.savez_compressed(output,omega=model.omega,y0=model.y0,**{f'long_y_{c["metadata"]["n"]}':c['y'] for c in long},
        **{f'long_v_{c["metadata"]["n"]}':c['v'] for c in long},reference_y=ref2['y'],reference_v=ref2['v'],
        **{f'short_y_{c["metadata"]["n"]}':c['y'] for c in short},**{f'short_v_{c["metadata"]["n"]}':c['v'] for c in short})
    source=Path(__file__)
    dependencies=['check_q9_nonlinear_restart.py','check_q9_full_space_equilibrium.py','check_q9_material_tangent.py']
    report=dict(protocol='docs/verification/q9-transient-accuracy-diagnostic.md',refined=refined,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={name:hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in dependencies},
        numpy=np.__version__,scipy=scipy.__version__,frequency_rad_s=model.omega.tolist(),frequency_ratio=float(model.omega[-1]/model.w1),
        equivalence=same,nonlinear_initial_acceleration=nonlinear_estimate,
        half_period_cases=[c['metadata'] for c in long],half_period_successive_differences=longdiff,
        half_period_iteration_difference=iteration_difference,half_period_strict=strict['metadata'],
        original_formulation_difference=original_difference,
        short_window_seconds=duration,oracles=[ref['metadata'],ref2['metadata']],reference_sensitivity=sensitivity,
        short_cases=[c['metadata'] for c in short],short_errors=short_errors,short_order=order,
        short_iteration_difference=short_iteration,short_strict=strict_short['metadata'],linear_midpoint_phase=phase,
        trajectory_file=str(output.relative_to(ROOT)),trajectory_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        checks=checks,all_criteria_pass=all(checks.values()),seconds=time.monotonic()-START,peak_process_working_bytes=plate.memory_bytes(),
        limitations='fixed 69-DOF isotropic one-element model; short-window oracle uses same internal force; no general long-time accuracy, damping or model truncation adopted')
    (ROOT/'docs/verification/q9-transient-accuracy-result.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,short_order=order,seconds=report['seconds'],sensitivity=sensitivity,short_error=short_errors[-1],long_iteration=iteration_difference)),flush=True)
    if not report['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':main()
