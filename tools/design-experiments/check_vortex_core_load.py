"""Limited core/self/physical steady-limit trial; see prewritten verification protocol."""
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad,quad_vec
import check_q9_plate_convergence as memory

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()


def json_scalar(value):
    if isinstance(value,np.generic):return value.item()
    raise TypeError(f'unsupported result type: {type(value).__name__}')


def kernel(points,a,b,eps,exclude=None):
    points=np.atleast_2d(points);a=np.atleast_2d(a);b=np.atleast_2d(b)
    if eps<0 or not all(np.all(np.isfinite(v)) for v in [points,a,b,np.asarray(eps)]):raise ValueError('invalid_input')
    dl=b-a;length=np.linalg.norm(dl,axis=1)
    if np.any(length<=1e-12):raise ValueError('degenerate_segment')
    e=dl/length[:,None];r=points[:,None,:]-a[None,:,:]
    s=np.einsum('pej,ej->pe',r,e);t=s-length[None,:]
    perp=r-s[:,:,None]*e[None,:,:];k2=np.sum(perp*perp,axis=2)+eps*eps
    aa=np.sqrt(s*s+k2);bb=np.sqrt(t*t+k2)
    excluded=np.zeros_like(s,dtype=bool) if exclude is None else exclude
    singular=(k2==0)&(s>=0)&(t<=0)&~excluded
    if np.any(singular):raise ValueError('singular_ideal_segment')
    integral=np.zeros_like(s);same=s*t>0
    denominator=aa*bb*(s*bb+t*aa)
    np.divide(length[None,:]*(s+t),denominator,out=integral,where=same)
    mixed=~same&(k2>0)
    value=np.divide(s,aa,out=np.zeros_like(s),where=aa>0)-np.divide(t,bb,out=np.zeros_like(t),where=bb>0)
    np.divide(value,k2,out=integral,where=mixed)
    result=np.cross(e[None,:,:],perp)*integral[:,:,None]/(4*np.pi)
    result[excluded]=0.
    return result


def reference_segment(p,a,b,eps):
    dl=b-a
    def fun(t):
        r=p-a-t*dl
        return np.cross(dl,r)/(np.dot(r,r)+eps*eps)**1.5/(4*np.pi)
    projection=float(np.dot(p-a,dl)/np.dot(dl,dl))
    result,error=quad_vec(fun,0.,1.,epsabs=1e-12,epsrel=1e-12,points=[projection] if 0<projection<1 else None)
    return result


def context_eval(p,a,b,eps,incident=False):
    if not all(np.all(np.isfinite(v)) for v in [p,a,b,eps]) or eps<0:raise ValueError('invalid_input')
    dl=b-a;length=np.linalg.norm(dl)
    if length<=1e-12:raise ValueError('degenerate_segment')
    projection=np.dot(p-a,dl)/(length*length)
    closest=a+np.clip(projection,0.,1.)*dl
    on_segment=np.linalg.norm(p-closest)<=1e-12
    if on_segment and not incident:raise ValueError('nonincident_contact')
    if incident and not on_segment:raise ValueError('incidence_geometry_mismatch')
    return kernel(p,a,b,eps,exclude=np.array([[incident]]))[0,0]


def scalar_tests():
    a=np.zeros(3);b=np.array([1.,0,0]);rows=[]
    for s in [-1.,0.,.5,1.,2.]:
        for d in [0.,.0002,.002,.02,.2,2.]:
            p=np.array([s,d,0.]);actual=kernel(p,a,b,.02)[0,0];ref=reference_segment(p,a,b,.02)
            absolute=float(np.linalg.norm(actual-ref));relative=absolute/max(np.linalg.norm(ref),1e-30)
            rows.append(dict(s=s,d=d,absolute=absolute,relative=relative,eligible=absolute<=1e-12 or relative<=1e-10))
    p=np.array([.4,.03,.02]);v=kernel(p,a,b,.02)[0,0]
    c=np.cos(.7);s=np.sin(.7);r=np.array([[c,-s,0],[s,c,0],[0,0,1.]])
    shift=np.array([.3,-.4,.5]);tests=[]
    def check(name,value,reference):
        diff=float(np.linalg.norm(value-reference));rel=diff/max(np.linalg.norm(reference),1e-30)
        tests.append(dict(name=name,absolute=diff,relative=rel,eligible=diff<=1e-12 or rel<=1e-11))
    check('orientation',kernel(p,b,a,.02)[0,0],-v)
    check('circulation_sign',-kernel(p,a,b,.02)[0,0],-v)
    check('rotation_translation',kernel(p@r.T+shift,a@r.T+shift,b@r.T+shift,.02)[0,0],r@v)
    for scale in [.1,10.]:check(f'scale_{scale}',kernel(scale*p,scale*a,scale*b,scale*.02)[0,0],v/scale)
    mid=.37*b
    check('subdivision',kernel(p,a,mid,.02)[0,0]+kernel(p,mid,b,.02)[0,0],v)
    pfar=np.array([2.,.0001,0.]);ideal=kernel(pfar,a,b,0.)[0,0];smooth=kernel(pfar,a,b,.01)[0,0]
    line_factor=.0001**2/np.sqrt(.0001**4+.01**4)
    extension=dict(finite_segment_distance=float(np.linalg.norm(pfar-b)),infinite_line_distance=.0001,
        integral_kernel_relative_to_ideal=float(np.linalg.norm(smooth-ideal)/np.linalg.norm(ideal)),
        line_correction_relative_to_ideal=float(abs(line_factor-1)))
    refused=[]
    for name,args in [('zero_length',(p,a,a,.02,False)),('negative_core',(p,a,b,-.02,False)),
        ('nonfinite',(np.array([np.nan,0,0]),a,b,.02,False)),('nonincident_contact',(np.array([.5,0,0]),a,b,.02,False)),
        ('wrong_incident_id',(np.array([.5,.1,0]),a,b,.02,True))]:
        try:context_eval(*args)
        except ValueError:refused.append(dict(name=name,detected=True))
        else:refused.append(dict(name=name,detected=False))
    own=context_eval(np.array([.5,0,0]),a,b,.02,True)
    # Two opposite copies of exactly one shared edge.
    near=np.array([.5,.005,0]);same=kernel(near,a,b,.01)[0,0]-kernel(near,a,b,.01)[0,0]
    mismatch=kernel(near,a,b,.01)[0,0]-kernel(near,a,b,.02)[0,0]
    shared=dict(same_core_cancellation=float(np.linalg.norm(same)),different_core_residual=float(np.linalg.norm(mismatch)))
    return dict(quadrature=rows,invariances=tests,extension=extension,refused=refused,own_straight_velocity=own.tolist(),shared=shared)


def ring_tests():
    rows=[]
    for eps in [.1,.05]:
        def fun(theta):
            a=2*np.sin(theta/2)**2
            return a/(2*a+eps*eps)**1.5/(4*np.pi)
        reference=2*quad(fun,0.,np.pi,epsabs=1e-13,epsrel=1e-13,points=[eps,4*eps])[0]
        axis_reference=1/(2*(1+.3**2+eps**2)**1.5)
        cases=[]
        for n in [64,128,256,512]:
            theta=2*np.pi*np.arange(n)/n
            a=np.column_stack([np.cos(theta),np.sin(theta),np.zeros(n)]);b=np.roll(a,-1,axis=0)
            self_velocity=kernel(np.array([1.,0,0]),a,b,eps).sum(axis=1)[0]
            axis_velocity=kernel(np.array([0.,0,.3]),a,b,eps).sum(axis=1)[0]
            cases.append(dict(n=n,self_velocity=self_velocity.tolist(),self_relative_error=float(abs(self_velocity[2]/reference-1)),
                axis_relative_error=float(abs(axis_velocity[2]/axis_reference-1))))
        order=float(np.log2(cases[-2]['self_relative_error']/cases[-1]['self_relative_error']))
        rows.append(dict(eps=eps,reference_self=reference,axis_reference=axis_reference,cases=cases,order=order,
            eligible=cases[-1]['self_relative_error']<=.01 and 1.7<=order<=2.3,whole_ring_zero_error=abs(reference)))
    return rows


def steady_case(n,ratio,span=10000.):
    rho=1.2;speed=8.;alpha=np.pi/180;wake=20*span
    x=np.r_[(np.arange(n)+.25)/n,1.,1.+wake]
    points=np.array([[xx,yy,0.] for xx in x for yy in [-span/2,span/2]])
    faces=[[2*i,2*i+1,2*(i+1)+1,2*(i+1)] for i in range(n+1)]
    edges=sorted({tuple(sorted((f[i],f[(i+1)%4]))) for f in faces for i in range(4)})
    index={edge:i for i,edge in enumerate(edges)};source=np.zeros((len(edges),n));body=source.copy()
    for i,f in enumerate(faces):
        column=min(i,n-1)
        for k in range(4):
            pair=(f[k],f[(k+1)%4]);edge=tuple(sorted(pair));sign=1 if pair==edge else -1
            source[index[edge],column]+=sign
            if i<n:body[index[edge],column]+=sign
    a=np.array([points[i] for i,j in edges]);b=np.array([points[j] for i,j in edges])
    eps=ratio/n
    collocation=np.column_stack([(np.arange(n)+.75)/n,np.zeros(n),np.zeros(n)])
    aic=kernel(collocation,a,b,eps)[:,:,2]@source
    rhs=np.full(n,-speed*np.sin(alpha));circulation=np.linalg.solve(aic,rhs)
    strengths=source@circulation;body_strength=body@circulation
    te=index[(2*n,2*n+1)]
    body_indices=np.array([i for i in range(len(edges)) if np.any(body[i]) and i!=te])
    centers=.5*(a[body_indices]+b[body_indices]);dl=b[body_indices]-a[body_indices]
    excluded=body_indices[:,None]==np.arange(len(edges))[None,:]
    induced=np.einsum('pej,e->pj',kernel(centers,a,b,eps,excluded),strengths)
    flow=speed*np.array([np.cos(alpha),0,np.sin(alpha)])
    forces=rho*body_strength[body_indices,None]*np.cross(flow+induced,dl)
    total=forces.sum(axis=0);moment=np.cross(centers,forces).sum(axis=0)
    lift_axis=np.array([-np.sin(alpha),0,np.cos(alpha)]);flow_axis=flow/speed
    lift=float(total@lift_axis/span);drag=float(total@flow_axis/span);my=float(moment[1]/span)
    gref=np.pi*speed*np.sin(alpha);lref=rho*speed*gref;mref=-rho*speed*np.cos(alpha)*gref/4
    extra=rho*body_strength[te]*np.cross(flow,b[te]-a[te])
    onlynormal=np.array([0.,0.,total[2]])
    residual=float(np.linalg.norm(aic@circulation-rhs)/np.linalg.norm(rhs))
    return dict(geometry_id='steady-quarter-threequarter-flat-v1',n=n,span=span,wake_length=wake,core_ratio=ratio,eps_m=eps,lift_per_span=lift,drag_per_span=drag,moment_y_per_span=my,
        circulation_total=float(circulation[-1]),reference_lift=lref,reference_moment=mref,
        lift_relative=float(abs(lift/lref-1)),moment_relative=float(abs(my/mref-1)),drag_over_reference_lift=float(abs(drag/lref)),
        nonpenetration_relative=residual,shared_te_strength=float(strengths[te]),
        normal_only_drag_over_reference_lift=float(abs(onlynormal@flow_axis/span/lref)),
        uncancelled_bound_te_lift_error=float(abs(((total+extra)@lift_axis/span-lref)/lref)))


def main():
    scalar=scalar_tests();rings=ring_tests();cases=[]
    for ratio in [0.,.001,.01]:
        for n in [8,16,32,64]:
            if time.monotonic()-START>120 or memory.memory_bytes()>1024**3:raise TimeoutError('experiment budget')
            row=steady_case(n,ratio);cases.append(row)
            print(json.dumps(dict(stage='steady',**row)),flush=True)
    alpha=np.pi/180;speed=8.
    integral=quad(lambda theta:speed*np.sin(alpha)*(1+np.cos(theta)),0,np.pi,epsabs=1e-12,epsrel=1e-12)[0]
    moment=quad(lambda theta:.5*(1-np.cos(theta))*speed*np.sin(alpha)*(1+np.cos(theta)),0,np.pi,epsabs=1e-12,epsrel=1e-12)[0]
    checks=dict(segment_quadrature=all(r['eligible'] for r in scalar['quadrature']),
        invariances=all(r['eligible'] for r in scalar['invariances']),invalid_geometry=all(r['detected'] for r in scalar['refused']),
        straight_self_zero=np.linalg.norm(scalar['own_straight_velocity'])==0,
        extension_far_field=scalar['extension']['integral_kernel_relative_to_ideal']<1e-3,
        negative_extension_correction=scalar['extension']['line_correction_relative_to_ideal']>.9,
        ring_self=all(r['eligible'] for r in rings),negative_whole_ring_zero=all(r['whole_ring_zero_error']>1e-3 for r in rings),
        shared_kernel=scalar['shared']['same_core_cancellation']<=1e-12 and scalar['shared']['different_core_residual']>1e-4,
        analytic_integrals=abs(integral-np.pi*speed*np.sin(alpha))<=1e-12 and abs(moment-integral/4)<=1e-12)
    for ratio in [0.,.001,.01]:
        fine=next(c for c in cases if c['core_ratio']==ratio and c['n']==64)
        checks[f'steady_{ratio}']=fine['lift_relative']<=1e-3 and fine['moment_relative']<=1e-3 and fine['drag_over_reference_lift']<=1e-4 and fine['nonpenetration_relative']<=1e-12
    checks['negative_normal_only']=all(c['normal_only_drag_over_reference_lift']>1e-3 for c in cases)
    checks['negative_bound_te']=all(c['uncancelled_bound_te_lift_error']>.5 for c in cases)
    width_cases=[steady_case(32,0.,span) for span in [1000.,2000.,10000.,20000.]]
    checks['width_wake_limit']=all(width_cases[-1][key]<width_cases[0][key] and width_cases[-1][key]<=1e-3 for key in ['lift_relative','moment_relative'])
    source=Path(__file__)
    result=dict(protocol='docs/verification/vortex-core-and-steady-load-experiment.md',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        numpy=np.__version__,scipy=scipy.__version__,memory_helper_sha256=hashlib.sha256(Path(memory.__file__).read_bytes()).hexdigest(),
        scalar=scalar,rings=rings,steady_cases=cases,width_wake_cases=width_cases,checks=checks,all_criteria_pass=bool(all(checks.values())),
        seconds=time.monotonic()-START,peak_process_working_bytes=memory.memory_bytes(),
        limits='same regularized-kernel mathematical references and high-AR steady 2D limit; no physical viscous-core, free-wake, camber/unsteady or ANCF coupled validation')
    (ROOT/'docs/verification/vortex-core-and-steady-load-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],rings=rings),default=json_scalar),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':main()
