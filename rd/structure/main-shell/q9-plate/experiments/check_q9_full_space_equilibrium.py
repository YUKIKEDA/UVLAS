"""Full 81-coordinate Q9 with 12 physical support conditions, not product code.
Protocol: rd/structure/main-shell/q9-plate/records/q9-full-space-equilibrium-experiment.md.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import hashlib
import itertools
import json
import time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import scipy
import skfem
from scipy.linalg import eigh
from numpy.polynomial.legendre import leggauss
from skfem import MeshHex,ElementHex2,ElementVector,Basis,BilinearForm,asm
import check_q9_prestress_equilibrium as prestress
import check_q9_material_tangent as material
import check_q9_curved_load as curved
import check_q9_interface_reparameterization as interface
import check_q9_independent_load as common
import check_q9_plate_convergence as plate

H=prestress.THICKNESS;RHO=prestress.DENSITY;YOUNG=prestress.YOUNG;NU=prestress.NU
Q0=prestress.Q0;PARAM=prestress.REF
W=interface.natural_scale(H,9)
LEFT_NODES=np.flatnonzero(common.nodes[:,0]==-1)
FIXED=np.array(sorted([9*i+3*k for i in LEFT_NODES for k in range(3)]+[9*7+1,9*7+2,9*3+2]))
FREE=np.setdiff1d(np.arange(81),FIXED)


def determinant_certificate(q):
    polys=curved.geometry_polynomials(q,H)
    jac=[[curved.derivative(p,a) for a in range(3)] for p in polys];det={}
    for perm in itertools.permutations(range(3)):
        inversions=sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        term=curved.multiply(curved.multiply(jac[0][perm[0]],jac[1][perm[1]]),jac[2][perm[2]])
        det=curved.add(det,curved.times(term,F((-1)**inversions)))
    return curved.bernstein_bounds(det)


def make_context(q0,rotation):
    pts,weights=leggauss(5);ijk=list(itertools.product(range(5),repeat=3))
    coords=np.array([[pts[i],pts[j],pts[k]] for i,j,k in ijk])
    w=np.array([weights[i]*weights[j]*weights[k] for i,j,k in ijk])
    s,g=interface.shape_tables(coords,0.,H/2)
    jac=np.einsum('ia,pib->pab',q0.reshape(27,3),g);det=np.linalg.det(jac)
    grad=np.linalg.solve(jac.transpose(0,2,1),g.transpose(0,2,1)).transpose(0,2,1)
    c=material.isotropic(YOUNG,NU);measure=w*det
    ctx=dict(w=measure,g=grad,axes=np.broadcast_to(rotation,(len(s),3,3)).copy(),
        e0=np.broadcast_to(prestress.E0,(len(s),3,3)).copy(),c=c,c4=material.fourth(c))
    mass=np.kron(RHO*(s.T@(measure[:,None]*s)),np.eye(3))
    return ctx,mass,s


def face_load(q0,rotation,sign):
    pts,weights=leggauss(5)
    coords=np.array([[sign,a,b] for a,b in itertools.product(pts,repeat=2)])
    w=np.array([a*b for a,b in itertools.product(weights,repeat=2)])
    s,g=interface.shape_tables(coords,0.,H/2)
    jac=np.einsum('ia,pib->pab',q0.reshape(27,3),g)
    area=sign*np.cross(jac[:,:,1],jac[:,:,2])
    nominal=rotation@np.diag([PARAM['p'],0.,0.])@rotation.T
    return (s.T@(area@nominal.T*w[:,None])).reshape(81)


def exact_q(rotation,translation):
    base=(Q0.reshape(27,3)@np.diag(PARAM['stretch'])).reshape(81)
    base=curved.transform(base,np.eye(3),np.array([.5*(PARAM['stretch'][0]-1),0.,0.]))
    return curved.transform(base,rotation,translation)


def solve(q0,p,ctx,s,load,deadline):
    rng=np.random.default_rng(20261030);direction=rng.standard_normal(len(FREE))
    dq=p@direction
    df=np.einsum('ia,pib->pab',dq.reshape(27,3),ctx['g'])
    direction/=np.max(np.linalg.norm(df,axis=(1,2)))
    x=2e-6*direction;history=[];initial_cert=determinant_certificate(q0+p@x)
    if not initial_cert['exact_positive']:return None,dict(eligible=False,reason='initial_geometry',initial_certificate=initial_cert)
    accepted_cert=initial_cert
    force_scale=np.linalg.norm(p.T@load)
    for iteration in range(26):
        if time.monotonic()>deadline:raise TimeoutError('case budget')
        q=q0+p@x;state=material.analytic_evaluation(q,ctx);r=p.T@(state['force']-load)
        k=p.T@state['k']@p;dx=np.linalg.solve(k,-r)
        residual=float(np.linalg.norm(r)/force_scale)
        correction=float(np.max(np.abs(s@(p@dx).reshape(27,3))))
        row=dict(iteration=iteration,residual=residual,candidate_physical_step=correction)
        history.append(row)
        if residual<=1e-10 and correction<=1e-10:
            return q,dict(eligible=True,history=history,initial_certificate=initial_cert,current_certificate=accepted_cert)
        if iteration==25:break
        slope=float(r@dx)
        if slope>=0:return None,dict(eligible=False,reason='non_descent',history=history)
        accepted=False
        for reduction in range(13):
            step=2.**(-reduction);candidate=q0+p@(x+step*dx)
            try:material.analytic_evaluation(candidate,ctx)
            except ValueError:continue
            change=prestress.potential_increment(q,candidate,ctx,load)
            if change<=1e-4*step*slope:
                cert=determinant_certificate(candidate)
                if not cert['exact_positive']:continue
                x+=step*dx;accepted_cert=cert;accepted=True
                row.update(step=step,potential_change=change);break
        if not accepted:return None,dict(eligible=False,reason='line_search_or_geometry',history=history)
    return None,dict(eligible=False,reason='iteration_limit',history=history)


def external_reference():
    mesh=MeshHex.init_tensor(np.array([-.5,.5]),np.array([-.5,.5]),np.array([-H/2,H/2]))
    basis=Basis(mesh,ElementVector(ElementHex2()),intorder=9)
    fstar=np.diag(PARAM['stretch']);sstar=np.diag(PARAM['stress'])
    mu=YOUNG/(2*(1+NU));lam=YOUNG*NU/((1+NU)*(1-2*NU))
    @BilinearForm
    def tangent(u,v,w):
        gu=np.einsum('ai,aj...->ij...',fstar,u.grad);gv=np.einsum('ai,aj...->ij...',fstar,v.grad)
        eu=.5*(gu+np.swapaxes(gu,0,1));ev=.5*(gv+np.swapaxes(gv,0,1))
        traceu=np.einsum('ii...->...',eu);tracev=np.einsum('ii...->...',ev)
        return lam*traceu*tracev+2*mu*np.sum(eu*ev,axis=(0,1))+np.einsum('ai...,ij,aj...->...',u.grad,sstar,v.grad)
    @BilinearForm
    def mass_form(u,v,w):return RHO*np.einsum('i...,i...->...',u,v)
    kp=asm(tangent,basis).toarray();mp=asm(mass_form,basis).toarray()
    transform=np.zeros((81,81));fixedp=[]
    for component,ids in enumerate(basis.split_indices()):
        for row in ids:
            x,y,z=basis.doflocs[:,row];z=H*round(2*z/H)/2
            matches=np.flatnonzero(np.max(np.abs(common.nodes-[2*x,2*y]),axis=1)<1e-12)
            if len(matches)!=1:raise ValueError('external_node_correspondence')
            node=int(matches[0]);transform[row,[9*node+component,9*node+3+component,9*node+6+component]]=[1,z,.5*z*z]
            atleft=abs(x+.5)<1e-12
            if atleft and (component==0 or (abs(y)<1e-12 and abs(z)<1e-12 and component in (1,2)) or (abs(y-.5)<1e-12 and abs(z)<1e-12 and component==2)):
                fixedp.append(row)
    freep=np.setdiff1d(np.arange(81),fixedp);p=np.eye(81)[:,FREE]*W[FREE][None,:]
    free_transform=transform[np.ix_(freep,np.arange(81))]@p
    mfree=mp[np.ix_(freep,freep)];kfree=kp[np.ix_(freep,freep)]
    d=1/np.sqrt(np.diag(mfree));vals,native=eigh(d[:,None]*kfree*d[None,:],d[:,None]*mfree*d[None,:])
    vectors=np.linalg.solve(free_transform,d[:,None]*native)
    mr=p.T@transform.T@mp@transform@p;kr=p.T@transform.T@kp@transform@p
    constraint=(transform[fixedp,:]*W[None,:])
    result=dict(k=kr,m=mr,values=vals,vectors=vectors,
        metadata=dict(native_dofs=int(basis.N),free_native_dofs=len(freep),volume_points=basis.X.shape[1],support_count=len(fixedp),
            support_ranks=common.rank(constraint),support_free_error=float(np.max(np.abs(transform[fixedp,:]@p)))))
    return result


def modes(k,m,external):
    values,vectors=eigh(k,m);frequency=np.sqrt(values)/(2*np.pi)
    reference_frequency=np.sqrt(external['values'])/(2*np.pi)
    errors=np.abs(frequency/reference_frequency-1)
    residuals=[float(np.linalg.norm(k@v-val*m@v)/((np.linalg.norm(k)+abs(val)*np.linalg.norm(m))*np.linalg.norm(v))) for val,v in zip(values,vectors.T)]
    orth=float(np.max(np.abs(vectors.T@m@vectors-np.eye(69))))
    groups=[];begin=0
    for end in range(1,70):
        if end<69 and (reference_frequency[end]-reference_frequency[end-1])/reference_frequency[end-1]<=1e-6:continue
        a=vectors[:,begin:end];b=external['vectors'][:,begin:end]
        gram=b.T@m@b;b=b@np.linalg.inv(np.linalg.cholesky(gram).T)
        singular=np.linalg.svd(a.T@m@b,compute_uv=False)
        sine=float(np.sqrt(max(0.,1-float(singular.min())**2)))
        groups.append(dict(first=begin,last=end-1,max_sin_angle=sine));begin=end
    comparisons=dict(stiffness=common.matrix_error(external['k'],k),mass=common.matrix_error(external['m'],m))
    passed=np.all(values>0) and max(residuals)<=1e-9 and orth<=1e-9 and float(errors.max())<=1e-6 and float(errors[:8].max())<=1e-7
    passed=passed and max(g['max_sin_angle'] for g in groups)<=1e-5 and all(v<=1e-8 for x in comparisons.values() for v in x.values())
    return dict(eligible=bool(passed),mode_count=len(values),frequencies_hz=frequency.tolist(),reference_frequencies_hz=reference_frequency.tolist(),
        maximum_frequency_relative=float(errors.max()),low8_frequency_relative=float(errors[:8].max()),
        eigenpair_residual_max=max(residuals),mass_orthogonality=orth,groups=groups,matrix_comparison=comparisons),values,vectors


def run_case(rotation,translation,external):
    start=time.monotonic();q0=curved.transform(Q0,rotation,translation);t=np.kron(np.eye(27),rotation)
    ctx,mass,s=make_context(q0,rotation)
    c=np.eye(81)[FIXED,:]@t.T;p=t[:,FREE]*W[FREE][None,:]
    refcert=determinant_certificate(q0)
    if not refcert['exact_positive']:raise ValueError('reference_geometry_unverified')
    load=face_load(q0,rotation,1);left=face_load(q0,rotation,-1)
    q,newton=solve(q0,p,ctx,s,load,start+120)
    record=dict(transformed=bool(np.linalg.norm(translation)>0),support_count=len(FIXED),free_coordinates=len(FREE),reference_certificate=refcert,newton=newton)
    if q is None:return dict(**record,eligible=False),None
    state=material.analytic_evaluation(q,ctx);residual=state['force']-load
    multiplier=-c@residual;reaction=-c.T@multiplier
    exact=exact_q(rotation,translation)
    expected_u=s@(exact-q0).reshape(27,3);error_u=s@(q-exact).reshape(27,3)
    f=np.einsum('ia,pib->pab',q.reshape(27,3),ctx['g']);fe=np.einsum('ia,pib->pab',exact.reshape(27,3),ctx['g'])
    green=.5*(np.einsum('pji,pjk->pik',f,f)-np.eye(3));expected_green=.5*(np.einsum('pji,pjk->pik',fe,fe)-np.eye(3))
    energy_expected=.5*YOUNG*((PARAM['stretch'][0]**2-1)/2-prestress.E0[0,0])**2*H
    rigid=curved.rigid_modes(q,translation);wrench_r=rigid.T@reaction;wrench_ext=rigid.T@load
    force_scale=PARAM['p']*H
    expected_force=-force_scale*rotation[:,0]
    reaction_error=common.force_error(state['k'],left,reaction)
    support_error=float(np.max(np.abs((t.T@(q-q0))[FIXED]/W[FIXED])))
    generalized_balance=float(np.linalg.norm(W*(residual-reaction))/np.linalg.norm(W*load))
    reduced_k=p.T@state['k']@p;reduced_m=p.T@mass@p
    modal,values,vectors=modes(reduced_k,reduced_m,external)
    coordinates_work=np.random.default_rng(1029).standard_normal(69)
    work=float(abs((p@coordinates_work)@reaction)/(np.linalg.norm(W*reaction)*np.linalg.norm(coordinates_work)))
    metrics=dict(displacement_relative=float(np.max(np.abs(error_u))/np.max(np.abs(expected_u))),
        green_relative=float(np.max(np.abs(green-expected_green))/np.max(np.abs(expected_green))),
        energy_relative=float(abs(state['energy']/energy_expected-1)),support_scaled_error=support_error,
        reaction_covector_relative=reaction_error,generalized_balance=generalized_balance,
        reaction_force_scaled=float(np.linalg.norm(wrench_r[:3]-expected_force)/force_scale),
        reaction_moment_scaled=float(np.linalg.norm(wrench_r[3:])/force_scale),
        combined_wrench_scaled=float(np.linalg.norm(wrench_r+wrench_ext)/force_scale),reaction_virtual_work=work)
    record.update(metrics=metrics,moment_reference_point=translation.tolist(),reaction_wrench=wrench_r.tolist(),external_wrench=wrench_ext.tolist(),modes=modal,
        seconds=time.monotonic()-start,peak_process_working_bytes=plate.memory_bytes())
    passed=all(metrics[v]<=1e-7 for v in ['displacement_relative','green_relative','energy_relative']) and support_error<=1e-12
    passed=passed and all(metrics[v]<=1e-8 for v in ['reaction_covector_relative','generalized_balance','reaction_force_scaled','reaction_moment_scaled','combined_wrench_scaled']) and work<=1e-10
    record['eligible']=bool(passed and modal['eligible'] and record['seconds']<=120 and record['peak_process_working_bytes']<=4*1024**3)
    return record,dict(q=q,q0=q0,ctx=ctx,mass=mass,p=p,c=c,k=reduced_k,m=reduced_m,state=state,load=load,reaction=reaction,values=values,vectors=vectors,
        rotation=t,wrench_r=wrench_r,wrench_ext=wrench_ext)


def transient_check(aux):
    k=aux['k'];m=aux['m'];omega=np.sqrt(aux['values'][0]);phi=aux['vectors'][:,0]
    # Normalize by the actual physical point displacement, not a mixed q norm.
    _,_,s=make_context(aux['q0'],np.eye(3))
    physical=s@(aux['p']@phi).reshape(27,3);amplitude=1e-5/np.max(np.abs(physical))
    period=2*np.pi/omega;records=[]
    for n in [80,160,320]:
        dt=period/n;a=m+.25*dt*dt*k
        # Equivalent constant midpoint update avoids rebuilding a factor each step.
        aa=np.linalg.solve(a,m-.25*dt*dt*k);bb=np.linalg.solve(a,dt*m)
        u=amplitude*phi;v=np.zeros(69);initial=.5*u@k@u;error=0.;energy=0.
        for step in range(1,2*n+1):
            unew=aa@u+bb@v;v=2*(unew-u)/dt-v;u=unew
            difference=u-amplitude*phi*np.cos(omega*step*dt)
            error=max(error,float(np.sqrt(max(0.,difference@m@difference))/amplitude))
            energy=max(energy,float(abs((.5*u@k@u+.5*v@m@v)/initial-1)))
        records.append(dict(steps_per_period=n,mass_norm_relative_error=error,energy_relative=energy))
    order=float(np.log2(records[-2]['mass_norm_relative_error']/records[-1]['mass_norm_relative_error']))
    return dict(eligible=bool(records[-1]['mass_norm_relative_error']<=1e-3 and 1.8<=order<=2.2 and all(r['energy_relative']<=1e-9 for r in records)),cases=records,last_order=order)


def negative_controls(aux):
    p=aux['p'];q=aux['q'];state=aux['state'];load=aux['load'];force_scale=np.linalg.norm(p.T@load)
    no_geo=eigh(p.T@state['kmat']@p,aux['m'],eigvals_only=True)
    missing=float(np.max(np.abs(np.sqrt(no_geo[:8]/aux['values'][:8])-1)))
    released=np.array([9*i+3*k for i in LEFT_NODES for k in [1,2]])
    mismatch=float(np.linalg.norm(W[released]*(state['force']-load)[released])/np.linalg.norm(W*load))
    wrong_reaction=float(np.linalg.norm(aux['wrench_ext']-aux['wrench_r'])/(PARAM['p']*H))
    zeroctx=dict(aux['ctx'],e0=np.zeros_like(aux['ctx']['e0']))
    zero=material.analytic_evaluation(q,zeroctx)
    reset=float(np.linalg.norm(p.T@(zero['force']-load))/force_scale)
    original=material.analytic_evaluation(aux['q0'],aux['ctx'])
    unbalanced=float(np.linalg.norm(p.T@(original['force']-load))/force_scale)
    all_left=np.concatenate([np.arange(9*i,9*i+9) for i in LEFT_NODES])
    clamp_gap=float(np.max(np.abs((q-aux['q0'])[all_left]/W[all_left])))
    return [dict(name='omitted_geometric_tangent',mismatch=missing,detected=missing>1e-3),
        dict(name='position_only_face_support',mismatch=mismatch,detected=mismatch>1e-4),
        dict(name='reaction_sign_error',mismatch=wrong_reaction,detected=wrong_reaction>1),
        dict(name='reset_initial_strain',mismatch=reset,detected=reset>.1),
        dict(name='unbalanced_reference_as_equilibrium',mismatch=unbalanced,detected=unbalanced>.1),
        dict(name='roller_replaced_by_full_clamp',mismatch=clamp_gap,detected=clamp_gap>1e-5)]


def main():
    start=time.monotonic();cases=[];auxiliaries=[];stopped=None;external=None
    try:
        external=external_reference()
        for r,t in [(np.eye(3),np.zeros(3)),(curved.ROTATION,curved.SHIFT)]:
            case,aux=run_case(r,t,external);cases.append(case);auxiliaries.append(aux)
            print(json.dumps(dict(transformed=case['transformed'],eligible=case['eligible'],seconds=case.get('seconds'))),flush=True)
            if not case['eligible']:stopped='case_failed';break
    except Exception as exc:stopped=str(exc)
    dynamic=None;negative=[];covariance=None
    if stopped is None and len(auxiliaries)==2:
        try:
            base,moved=auxiliaries
            covariance=dict(displacement=float(np.max(np.abs((moved['rotation'].T@(moved['q']-moved['q0']))-(base['q']-base['q0'])))/np.max(np.abs(base['q']-base['q0']))),
                reaction=common.force_error(base['state']['k'],base['reaction'],moved['rotation'].T@moved['reaction']),
                stiffness=common.matrix_error(base['k'],moved['k']),mass=common.matrix_error(base['m'],moved['m']))
            dynamic=transient_check(base);negative=negative_controls(base)
        except Exception as exc:stopped=str(exc)
    cov_ok=bool(covariance and covariance['displacement']<=1e-7 and covariance['reaction']<=1e-8 and all(v<=1e-8 for key in ['stiffness','mass'] for v in covariance[key].values()))
    elapsed=time.monotonic()-start
    passed=stopped is None and len(cases)==2 and all(c['eligible'] for c in cases) and bool(dynamic and dynamic['eligible']) and cov_ok and len(negative)==6 and all(c['detected'] for c in negative) and elapsed<=240
    source=Path(__file__).resolve();modules=[prestress,material,curved,interface,common,plate]
    report=dict(protocol='rd/structure/main-shell/q9-plate/records/q9-full-space-equilibrium-experiment.md',numpy=np.__version__,scipy=scipy.__version__,skfem=skfem.__version__,
        external_reference_metadata=external['metadata'] if external else None,cases=cases,coordinate_covariance=covariance,
        linear_response=dynamic,negative_controls=negative,seconds=elapsed,stopped=stopped,all_criteria_pass=bool(passed),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),dependencies=[dict(path=Path(m.__file__).name,sha256=hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest()) for m in modules],
        limitation='One isotropic element, 12 linear supports and manufactured homogeneous equilibrium; external reference supplies basis/assembly of analytic prestress second variation, not an independent nonlinear material solver')
    output=rd_paths.ROOT/'build/analysis/q9-full-space-equilibrium-result.json';output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(all_criteria_pass=bool(passed),stopped=stopped,seconds=elapsed)),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
