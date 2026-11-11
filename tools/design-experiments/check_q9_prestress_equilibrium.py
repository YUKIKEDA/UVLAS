"""Prestress equilibrium -> modes -> small response, limited 4-coordinate fixture.
Not a general UVLAS nonlinear/transient implementation.
Protocol: docs/verification/q9-prestress-equilibrium-experiment.md.
"""
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path
import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.integrate import solve_ivp
from scipy.optimize import linear_sum_assignment
import check_q9_material_tangent as material
import check_q9_curved_load as curved
import check_q9_interface_reparameterization as interface
import check_q9_independent_load as common
import check_q9_plate_convergence as plate

YOUNG=70e9;NU=.3;DENSITY=2700.;THICKNESS=.02;VOLUME=THICKNESS;MASS=DENSITY*VOLUME
E0=np.diag([-.001,.0003,.0003])
OMEGA=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,0.]])


def reference_parameters():
    lam=1.002;elastic=(lam*lam-1)/2-E0[0,0]
    lateral=np.sqrt(1+2*(E0[1,1]-NU*elastic))
    stretch=np.array([lam,lateral,lateral]);stress=np.array([YOUNG*elastic,0.,0.])
    mu=YOUNG/(2*(1+NU));lame=YOUNG*NU/((1+NU)*(1-2*NU))
    normal=np.eye(3)*2*mu+lame
    nominal=lam*stress[0]
    return dict(stretch=stretch,stress=stress,c=normal,p=nominal,force_scale=VOLUME*nominal)


REF=reference_parameters()
Q0=curved.geometry('RECT')


def map_state(x):
    rotation=material.rz(x[3]);f=rotation@np.diag(x[:3]);base=Q0.reshape(27,3)
    q=(base@f.T).reshape(81);g=np.zeros((81,4));h=np.zeros((81,4,4))
    for i in range(3):
        g[:,i]=np.outer(base[:,i],rotation[:,i]).reshape(81)
        cross=(g[:,i].reshape(27,3)@OMEGA.T).reshape(81)
        h[:,i,3]=h[:,3,i]=cross
    g[:,3]=(q.reshape(27,3)@OMEGA.T).reshape(81)
    h[:,3,3]=(g[:,3].reshape(27,3)@OMEGA.T).reshape(81)
    return q,g,h


def element_context():
    p,w=leggauss(5);indexes=list(itertools.product(range(5),repeat=3))
    coords=np.array([[p[i],p[j],p[k]] for i,j,k in indexes])
    weights=np.array([w[i]*w[j]*w[k] for i,j,k in indexes])
    s,g=interface.shape_tables(coords,0.,THICKNESS/2)
    jac=np.einsum('ia,pib->pab',Q0.reshape(27,3),g);det=np.linalg.det(jac)
    grad=np.linalg.solve(jac.transpose(0,2,1),g.transpose(0,2,1)).transpose(0,2,1)
    measure=weights*det;c=material.isotropic(YOUNG,NU)
    ctx=dict(w=measure,g=grad,axes=np.broadcast_to(np.eye(3),(len(s),3,3)).copy(),
        e0=np.broadcast_to(E0,(len(s),3,3)).copy(),c=c,c4=material.fourth(c))
    mass=np.kron(DENSITY*(s.T@(measure[:,None]*s)),np.eye(3))
    return ctx,mass,s


def nominal_load():
    p,w=leggauss(5);result=np.zeros(81);nominal=np.diag([REF['p'],0.,0.])
    for sign in [-1,1]:
        coords=np.array([[sign,a,b] for a,b in itertools.product(p,repeat=2)])
        weights=np.array([a*b for a,b in itertools.product(w,repeat=2)])
        s,g=interface.shape_tables(coords,0.,THICKNESS/2)
        jac=np.einsum('ia,pib->pab',Q0.reshape(27,3),g)
        area=sign*np.cross(jac[:,:,1],jac[:,:,2])
        traction=area@nominal.T*weights[:,None]
        result+=(s.T@traction).reshape(81)
    return result


def evaluate(x,ctx,mass,load):
    if not np.all(np.isfinite(x)) or np.any(x[:3]<=0):raise ValueError('outside_positive_stretch_branch')
    q,g,h=map_state(x);state=material.analytic_evaluation(q,ctx)
    residual=state['force']-load
    rx=g.T@residual;basic=g.T@state['k']@g
    curvature=np.einsum('i,iab->ab',residual,h)
    return dict(q=q,g=g,h=h,element=state,residual_q=residual,residual=rx,
        stiffness=basic+curvature,basic_stiffness=basic,curvature_stiffness=curvature,
        mass=g.T@mass@g,potential=float(state['energy']-load@(q-Q0)))


def reference_matrices(x):
    stretches=x[:3];theta=x[3]
    elastic=(stretches*stretches-1)/2-np.diag(E0);stress=REF['c']@elastic
    residual=np.zeros(4);residual[:3]=VOLUME*stretches*stress
    residual[0]-=REF['force_scale']*np.cos(theta)
    residual[3]=REF['force_scale']*stretches[0]*np.sin(theta)
    k=np.zeros((4,4));k[:3,:3]=VOLUME*(np.diag(stress)+stretches[:,None]*REF['c']*stretches[None,:])
    k[0,3]=k[3,0]=REF['force_scale']*np.sin(theta)
    k[3,3]=REF['force_scale']*stretches[0]*np.cos(theta)
    md=MASS*np.array([1/12,1/12,THICKNESS**2/12,(stretches[0]**2+stretches[1]**2)/12])
    return residual,k,np.diag(md)


def potential_increment(before,after,ctx,load):
    dq=after-before
    f=np.einsum('ia,pib->pab',before.reshape(27,3),ctx['g'])
    df=np.einsum('ia,pib->pab',dq.reshape(27,3),ctx['g'])
    e=.5*(np.einsum('pji,pjk->pik',f,f)-np.eye(3))
    de=.5*(np.einsum('pji,pjk->pik',f,df)+np.einsum('pji,pjk->pik',df,f)+np.einsum('pji,pjk->pik',df,df))
    axes=ctx['axes']
    em=np.einsum('pji,pjk,pkl->pil',axes,e,axes)-ctx['e0']
    dem=np.einsum('pji,pjk,pkl->pil',axes,de,axes)
    ev=material.engineering(em);dev=material.engineering(dem)
    internal=np.einsum('p,pa,ab,pb->',ctx['w'],ev,ctx['c'],dev)
    internal+=.5*np.einsum('p,pa,ab,pb->',ctx['w'],dev,ctx['c'],dev)
    return float(internal-load@dq)


def solve_equilibrium(ctx,mass,load,s,deadline,stable_increment=False):
    x=np.array([1.,1.,1.,.05]);history=[]
    for iteration in range(21):
        if time.monotonic()>deadline:raise TimeoutError('experiment budget')
        state=evaluate(x,ctx,mass,load)
        correction=np.linalg.solve(state['stiffness'],-state['residual'])
        dq=state['g']@correction
        physical=float(np.max(np.abs(s@dq.reshape(27,3))))
        residual=float(np.linalg.norm(state['residual'])/REF['force_scale'])
        row=dict(iteration=iteration,residual=residual,candidate_physical_step=physical,
                 potential=state['potential'],x=x.tolist())
        history.append(row)
        if residual<=1e-10 and physical<=1e-10:
            state=evaluate(x,ctx,mass,load)
            return x,state,dict(eligible=True,history=history)
        if iteration==20:break
        slope=float(state['residual']@correction)
        if slope>=0:return x,state,dict(eligible=False,reason='non_descent_Newton_direction',history=history)
        accepted=False
        for reduction in range(13):
            step=2.**(-reduction);trial=x+step*correction
            try:next_state=evaluate(trial,ctx,mass,load)
            except ValueError:continue
            change=(potential_increment(state['q'],next_state['q'],ctx,load) if stable_increment
                    else next_state['potential']-state['potential'])
            accepted_test=(change<=1e-4*step*slope if stable_increment else
                           next_state['potential']<=state['potential']+1e-4*step*slope)
            if accepted_test:
                x=trial;row['accepted_step']=step;row['potential_change']=change
                row['total_value_difference']=next_state['potential']-state['potential'];accepted=True;break
        if not accepted:return x,state,dict(eligible=False,reason='line_search_failed',history=history)
    return x,state,dict(eligible=False,reason='Newton_iteration_limit',history=history)


def modal_checks(x,state):
    _,kr,mr=reference_matrices(x);k=state['stiffness'];m=state['mass']
    eigenvalues,vectors=eigh(k,m);ref_values,ref_vectors=eigh(kr,mr)
    residuals=[]
    for value,v in zip(eigenvalues,vectors.T):
        residuals.append(float(np.linalg.norm(k@v-value*(m@v))/((np.linalg.norm(k)+abs(value)*np.linalg.norm(m))*np.linalg.norm(v))))
    orthogonal=float(np.max(np.abs(vectors.T@m@vectors-np.eye(4))))
    frequency_error=float(np.max(np.abs(np.sqrt(eigenvalues/ref_values)-1)))
    ref_theta=np.array([0.,0.,0.,1/np.sqrt(mr[3,3])])
    mac=np.abs(vectors.T@m@ref_theta)**2;theta_index=int(np.argmax(mac))
    expected=REF['force_scale']*REF['stretch'][0]/(MASS*(REF['stretch'][0]**2+REF['stretch'][1]**2)/12)
    theta_error=float(abs(eigenvalues[theta_index]/expected-1))
    differences=dict(stiffness=common.matrix_error(kr,k),mass=common.matrix_error(mr,m))
    passed=np.all(eigenvalues>0) and max(residuals)<=1e-10 and orthogonal<=1e-10 and frequency_error<=1e-8 and theta_error<=1e-8 and mac[theta_index]>=1-1e-8
    passed=passed and all(v<=1e-9 for field in differences.values() for v in field.values())
    return dict(eligible=bool(passed),squared_frequencies=eigenvalues.tolist(),frequencies_hz=(np.sqrt(eigenvalues)/(2*np.pi)).tolist(),
        reference_frequencies_hz=(np.sqrt(ref_values)/(2*np.pi)).tolist(),reference_frequency_relative=frequency_error,
        eigenpair_residuals=residuals,mass_orthogonality=orthogonal,theta_mode_index=theta_index,theta_mac=float(mac[theta_index]),
        theta_squared_frequency_relative=theta_error,matrix_comparison=differences),float(np.sqrt(expected))


def midpoint_checks(k,m,omega):
    records=[];amplitude=.001;period=2*np.pi/omega
    frequencies=np.sqrt(eigh(k,m,eigvals_only=True))
    for resolution in [40,80,160,320]:
        step=period/resolution;eye=np.eye(4)
        e=np.block([[eye,-step/2*eye],[step/2*k,m]])
        a=np.block([[eye,step/2*eye],[-step/2*k,m]])
        update=np.linalg.solve(e,a)
        mu=np.linalg.eigvals(update)
        plus=(1+1j*frequencies*step/2)/(1-1j*frequencies*step/2)
        expected=np.concatenate([plus,np.conjugate(plus)])
        cost=np.abs(mu[:,None]-expected[None,:]);rows,cols=linear_sum_assignment(cost)
        eigen_error=float(np.max(cost[rows,cols]));unit_error=float(np.max(np.abs(np.abs(mu)-1)))
        theta_expected=(1+1j*omega*step/2)/(1-1j*omega*step/2)
        selected=mu[int(np.argmin(np.abs(mu-theta_expected)))]
        phase_error=float(abs(np.angle(selected)-2*np.arctan(omega*step/2)))
        z=np.zeros(8);z[3]=amplitude
        initial_energy=float(.5*z[:4]@k@z[:4]+.5*z[4:]@m@z[4:])
        error=0.;energy_error=0.
        for n in range(1,4*resolution+1):
            z=update@z
            error=max(error,abs(z[3]-amplitude*np.cos(omega*n*step))/amplitude)
            energy=float(.5*z[:4]@k@z[:4]+.5*z[4:]@m@z[4:])
            energy_error=max(energy_error,abs(energy/initial_energy-1))
        records.append(dict(steps_per_period=resolution,step=step,max_theta_relative=error,
            max_energy_relative=energy_error,discrete_eigenvalue_distance=eigen_error,unit_circle_error=unit_error,
            theta_phase_error=phase_error,continuous_discrete_phase_difference=float(abs(2*np.arctan(omega*step/2)-omega*step))))
    order=float(np.log2(records[-2]['max_theta_relative']/records[-1]['max_theta_relative']))
    passed=records[-1]['max_theta_relative']<=1e-3 and 1.8<=order<=2.2
    passed=passed and all(r['max_energy_relative']<=1e-10 and r['discrete_eigenvalue_distance']<=1e-9 and r['unit_circle_error']<=1e-9 and r['theta_phase_error']<=1e-9 for r in records)
    return dict(eligible=bool(passed),cases=records,last_order=order)


def stable_reference_gradient(displacement):
    u=displacement[:3];theta=displacement[3];lam=REF['stretch']+u
    de=REF['stretch']*u+.5*u*u;ds=REF['c']@de
    r=np.zeros(4);r[:3]=VOLUME*(lam*ds+REF['stress']*u)
    r[0]+=REF['force_scale']*2*np.sin(theta/2)**2
    r[3]=REF['force_scale']*lam[0]*np.sin(theta)
    return r


def reference_inertia(displacement,velocity):
    lam=REF['stretch']+displacement[:3]
    md=MASS*np.array([1/12,1/12,THICKNESS**2/12,(lam[0]**2+lam[1]**2)/12])
    h=np.zeros(4);h[:2]=-MASS/12*lam[:2]*velocity[3]**2
    h[3]=2*MASS/12*np.dot(lam[:2],velocity[:2])*velocity[3]
    return md,h


def nonlinear_oracle(ctx,mass,load,omega,deadline):
    period=2*np.pi/omega;times=np.linspace(0,period,401);records=[]
    def rhs(_,z):
        if time.monotonic()>deadline:raise TimeoutError('nonlinear reference budget')
        md,h=reference_inertia(z[:4],z[4:])
        return np.concatenate([z[4:],-(h+stable_reference_gradient(z[:4]))/md])
    for amplitude in [.01,.005,.0025]:
        initial=np.zeros(8);initial[3]=amplitude
        sol=solve_ivp(rhs,(0,period),initial,method='DOP853',t_eval=times,rtol=1e-11,atol=1e-13)
        if not sol.success or len(sol.t)!=len(times):raise RuntimeError('nonlinear_reference_failed')
        error=float(np.max(np.abs(sol.y[3]-amplitude*np.cos(omega*times)))/amplitude)
        residual=0.;mass_difference=0.;gradient_difference=0.;convective_difference=0.
        for idx in np.linspace(0,400,9,dtype=int):
            state=sol.y[:,idx];x=np.concatenate([REF['stretch']+state[:3],[state[3]]])
            fe=evaluate(x,ctx,mass,load);velocity=state[4:];acceleration=rhs(sol.t[idx],state)[4:]
            curvature_acceleration=np.einsum('iab,a,b->i',fe['h'],velocity,velocity)
            convective=fe['g'].T@mass@curvature_acceleration
            md,hr=reference_inertia(state[:4],velocity)
            residual=max(residual,float(np.linalg.norm(fe['mass']@acceleration+convective+fe['residual'])/REF['force_scale']))
            mass_difference=max(mass_difference,float(np.linalg.norm(fe['mass']-np.diag(md))/np.linalg.norm(md)))
            gradient_difference=max(gradient_difference,float(np.linalg.norm(fe['residual']-stable_reference_gradient(state[:4]))/REF['force_scale']))
            convective_difference=max(convective_difference,float(np.linalg.norm(convective-hr)/REF['force_scale']))
        passed=max(residual,gradient_difference,convective_difference)<=1e-9 and mass_difference<=1e-10
        records.append(dict(amplitude=amplitude,linear_limit_relative=error,nfev=int(sol.nfev),
            projected_ode_residual=residual,mass_reference_relative=mass_difference,
            gradient_reference_scaled=gradient_difference,convective_reference_scaled=convective_difference,eligible=bool(passed)))
    orders=[float(np.log2(records[i]['linear_limit_relative']/records[i+1]['linear_limit_relative'])) for i in range(2)]
    passed=all(r['eligible'] for r in records) and all(1.8<=o<=2.2 for o in orders)
    return dict(eligible=bool(passed),cases=records,amplitude_orders=orders)


def main():
    start=time.monotonic();deadline=start+120;ctx,mass,s=element_context();load=nominal_load()
    stable_increment='--stable-increment' in sys.argv
    filename='q9-prestress-equilibrium-refined-result.json' if stable_increment else 'q9-prestress-equilibrium-result.json'
    def model_fingerprint():
        fields=[Q0,ctx['c'],ctx['axes'],ctx['e0'],mass,load]
        return hashlib.sha256(b'affine-stretches-zrotation-v1'+b''.join(v.tobytes() for v in fields)).hexdigest()
    fingerprint=model_fingerprint()
    x,eq,newton=solve_equilibrium(ctx,mass,load,s,deadline,stable_increment=stable_increment)
    if not newton['eligible']:
        failure=dict(protocol='docs/verification/q9-prestress-equilibrium-experiment.md',all_criteria_pass=False,
            stopped='equilibrium_failed',stable_increment=stable_increment,equilibrium=newton,
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        path=Path(__file__).resolve().parents[2]/'build/analysis'/filename
        path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(failure,indent=2)+'\n',encoding='utf8')
        print(json.dumps(failure),flush=True)
        raise SystemExit(1)
    expected=np.concatenate([REF['stretch'],[0.]])
    scale=interface.natural_scale(THICKNESS,9)
    full_residual=float(np.linalg.norm(scale*eq['residual_q'])/np.linalg.norm(scale*load))
    equilibrium=dict(newton=newton,x=x.tolist(),expected=expected.tolist(),coordinate_error=float(np.max(np.abs(x-expected))),full_residual=full_residual)
    equilibrium['eligible']=bool(newton['eligible'] and equilibrium['coordinate_error']<=1e-10 and full_residual<=1e-8)
    print(json.dumps(dict(stage='equilibrium',eligible=equilibrium['eligible'],iterations=len(newton['history'])-1)),flush=True)
    modes,omega=modal_checks(x,eq)
    linear=midpoint_checks(eq['stiffness'],eq['mass'],omega)
    print(json.dumps(dict(stage='modes_and_linear_response',modes=modes['eligible'],linear=linear['eligible'])),flush=True)
    unbalanced_x=np.array([1.,1.,1.,0.]);unbalanced=evaluate(unbalanced_x,ctx,mass,load)
    accel=np.linalg.solve(unbalanced['mass'],-unbalanced['residual'])
    accel_res=float(np.linalg.norm(unbalanced['mass']@accel+unbalanced['residual'])/np.linalg.norm(unbalanced['residual']))
    initialization=dict(x=unbalanced_x.tolist(),unchanged=bool(np.array_equal(unbalanced_x,[1.,1.,1.,0.])),
        acceleration=accel.tolist(),equilibrium_residual=float(np.linalg.norm(unbalanced['residual'])/REF['force_scale']),
        acceleration_residual=accel_res,eligible=bool(accel_res<=1e-10 and np.linalg.norm(accel)>0))
    nonlinear=nonlinear_oracle(ctx,mass,load,omega,deadline)
    print(json.dumps(dict(stage='nonlinear_reference',eligible=nonlinear['eligible'],orders=nonlinear['amplitude_orders'])),flush=True)
    negatives=[]
    no_geo=eq['g'].T@eq['element']['kmat']@eq['g']+eq['curvature_stiffness']
    missing=abs(no_geo[3,3]/eq['mass'][3,3]/omega**2-1)
    negatives.append(dict(name='omitted_geometric_tangent',mismatch=float(missing),detected=bool(missing>.9)))
    zeroctx=dict(ctx,e0=np.zeros_like(ctx['e0']));zero_state=evaluate(x,zeroctx,mass,load)
    mismatch=float(np.linalg.norm(zero_state['residual'])/REF['force_scale'])
    negatives.append(dict(name='reset_initial_strain',mismatch=mismatch,detected=mismatch>.1))
    off=evaluate(np.array([1.,1.,1.,.2]),ctx,mass,load);_,ref_k,_=reference_matrices(np.array([1.,1.,1.,.2]))
    mismatch=float(abs(off['basic_stiffness'][3,3]/ref_k[3,3]-1))
    negatives.append(dict(name='omitted_constraint_curvature',mismatch=mismatch,detected=mismatch>1e-3))
    negatives.append(dict(name='unbalanced_reference_as_equilibrium',mismatch=initialization['equilibrium_residual'],detected=initialization['equilibrium_residual']>.1))
    phase=linear['cases'][0]['continuous_discrete_phase_difference']
    negatives.append(dict(name='continuous_discrete_frequency_confusion',mismatch=phase,detected=phase>1e-4))
    test_x=expected.copy();test_x[3]=.01;v=np.array([.001,-.002,0.,5.]);fe=evaluate(test_x,ctx,mass,load)
    hvec=fe['g'].T@mass@np.einsum('iab,a,b->i',fe['h'],v,v)
    mismatch=float(np.linalg.norm(hvec)/REF['force_scale'])
    negatives.append(dict(name='omitted_velocity_inertia',mismatch=mismatch,detected=mismatch>1e-6))
    final_fingerprint=model_fingerprint()
    elapsed=time.monotonic()-start;memory=plate.memory_bytes()
    passed=equilibrium['eligible'] and modes['eligible'] and linear['eligible'] and initialization['eligible'] and nonlinear['eligible'] and all(n['detected'] for n in negatives)
    passed=passed and fingerprint==final_fingerprint and elapsed<=120 and memory<=4*1024**3
    source=Path(__file__).resolve();modules=[material,curved,interface,common,plate]
    report=dict(protocol='docs/verification/q9-prestress-equilibrium-experiment.md',stable_increment=stable_increment,scipy=scipy.__version__,numpy=np.__version__,
        equilibrium=equilibrium,modes=modes,linear_response=linear,nonequilibrium_initialization=initialization,
        nonlinear_reference=nonlinear,negative_controls=negatives,unchanged_model=fingerprint==final_fingerprint,
        model_fingerprint=fingerprint,fingerprint_fields=['reference_geometry','elasticity','material_axes','initial_strain','mass','nominal_load','kinematic_map_id'],
        seconds=elapsed,peak_process_working_bytes=memory,all_criteria_pass=bool(passed),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependencies=[dict(path=Path(m.__file__).name,sha256=hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest()) for m in modules],
        limitation='4-coordinate affine kinematic fixture embedded in Q9; nonlinear orbit is an independent analytic reduced model; not full-DOF equilibrium/transient, wing modes, aeroelasticity or restart')
    output=source.parents[2]/'build/analysis'/filename;output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(all_criteria_pass=passed,seconds=elapsed,source_sha256=report['source_sha256'])),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
