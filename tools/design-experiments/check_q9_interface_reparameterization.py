"""Two-element chart/interface experiment, not UVLAS production code.
Protocol: docs/verification/q9-interface-reparameterization-experiment.md.
"""
import hashlib
import itertools
import json
import time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
import scipy
import check_q9_independent_load as common
import check_q9_curved_load as curved
import check_q9_plate_convergence as plate

NDOF=135
KAPPA=.001


def chart_transform(alpha,beta,reverse_x=True):
    if not np.isfinite(alpha) or not np.isfinite(beta) or alpha==0:
        raise ValueError('noninvertible_or_nonfinite_thickness_chart')
    sign=-1 if alpha<0 and reverse_x else 1
    h=np.array([[1,beta,.5*beta*beta],[0,alpha,alpha*beta],[0,0,alpha*alpha]])
    block=np.kron(h,np.eye(3));t=np.zeros((81,81));permutation=np.zeros_like(t)
    for i,(xi,eta) in enumerate(common.nodes):
        j=int(np.flatnonzero(np.max(np.abs(common.nodes-[sign*xi,eta]),axis=1)<1e-14)[0])
        t[9*i:9*i+9,9*j:9*j+9]=block
        permutation[9*i:9*i+9,9*j:9*j+9]=np.eye(9)
    return t,permutation,sign


def canonical_layers(h):
    return [dict(zlo=-h/2,zhi=-h/6,young=1.,rho=1.,id='lower'),
            dict(zlo=-h/6,zhi=h/2,young=2.,rho=3.,id='upper')]


def local_layers(h,alpha,beta,wrong_order=False):
    result=[]
    for layer in canonical_layers(h):
        lo,hi=sorted([(layer['zlo']-beta)/alpha,(layer['zhi']-beta)/alpha])
        result.append(dict(layer,zlo=lo,zhi=hi))
    result.sort(key=lambda item:item['zlo'])
    if wrong_order:
        for actual,wrong in zip(result,canonical_layers(h)):
            actual.update(young=wrong['young'],rho=wrong['rho'],id=wrong['id'])
    return result


def canonical_geometry(e):
    q=np.zeros((9,3,3))
    for i,(xi,eta) in enumerate(common.nodes):
        q[i,0]=[-.25+.5*e+.25*xi,.5*eta,0.]
        q[i,1]=[0,0,1]
    return q.reshape(81)


def gather(e):
    result=np.zeros((81,NDOF))
    for i,(xi,eta) in enumerate(common.nodes.astype(int)):
        node=(eta+1)*5+2*e+xi+1
        result[9*i:9*i+9,9*node:9*node+9]=np.eye(9)
    return result


def global_geometry():
    q=np.zeros((15,3,3))
    for j,y in enumerate(np.linspace(-.5,.5,3)):
        for i,x in enumerate(np.linspace(-.5,.5,5)):
            q[j*5+i,0]=[x,y,0]
            q[j*5+i,1]=[0,0,1]
    return q.reshape(NDOF)


def scalar_shape(xi,eta,z):
    n,nx,ny=common.shape(xi,eta)
    s=np.column_stack((n,z*n,.5*z*z*n)).reshape(27)
    gx=np.column_stack((nx,z*nx,.5*z*z*nx)).reshape(27)
    gy=np.column_stack((ny,z*ny,.5*z*z*ny)).reshape(27)
    gz=np.column_stack((np.zeros(9),n,z*n)).reshape(27)
    return s,gx,gy,gz


def shape_tables(points,zcenter,zhalf):
    shapes=[];gradients=[]
    for xi,eta,t in points:
        s,gx,gy,gz=scalar_shape(xi,eta,zcenter+zhalf*t)
        shapes.append(s);gradients.append(np.column_stack((gx,gy,zhalf*gz)))
    return np.array(shapes),np.array(gradients)


def geometry_certificate(q,zlo,zhi):
    # Exact polynomials for the stored off-centre chart, not a float recentering.
    center=(F.from_float(zlo)+F.from_float(zhi))/2
    half=(F.from_float(zhi)-F.from_float(zlo))/2
    zp={(0,0,0):center,(0,0,1):half}
    thickness=[{(0,0,0):F(1)},zp,curved.times(curved.multiply(zp,zp),F(1,2))]
    lag={-1:[F(0),-F(1,2),F(1,2)],0:[F(1),F(0),-F(1)],1:[F(0),F(1,2),F(1,2)]}
    polys=[{}, {}, {}];q=q.reshape(9,3,3)
    for i,(xi,eta) in enumerate(common.nodes.astype(int)):
        plane={(a,b,0):v*w for a,v in enumerate(lag[xi]) for b,w in enumerate(lag[eta]) if v*w}
        for kind in range(3):
            shape=curved.multiply(plane,thickness[kind])
            for axis in range(3):
                polys[axis]=curved.add(polys[axis],curved.times(shape,F.from_float(float(q[i,kind,axis]))))
    jac=[[curved.derivative(p,axis) for axis in range(3)] for p in polys];det={}
    for p in itertools.permutations(range(3)):
        inversions=sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
        term=curved.multiply(curved.multiply(jac[0][p[0]],jac[1][p[1]]),jac[2][p[2]])
        det=curved.add(det,curved.times(term,F((-1)**inversions)))
    bounds=curved.bernstein_bounds(det);bounds['volume']=float(curved.integrate(det))
    return bounds


def evaluate_layer(q,layer,order=7):
    points,weights=leggauss(order)
    indexes=list(itertools.product(range(order),repeat=3))
    coordinates=np.array([[points[i],points[j],points[k]] for i,j,k in indexes])
    w=np.array([weights[i]*weights[j]*weights[k] for i,j,k in indexes])
    center=(layer['zlo']+layer['zhi'])/2;half=(layer['zhi']-layer['zlo'])/2
    s,g=shape_tables(coordinates,center,half)
    jac=np.einsum('ia,pib->pab',q.reshape(27,3),g);det=np.linalg.det(jac)
    if np.any(det<=0):raise ValueError('nonpositive_reference_jacobian')
    grad=np.linalg.solve(jac.transpose(0,2,1),g.transpose(0,2,1)).transpose(0,2,1)
    b=np.zeros((len(s),6,81))
    b[:,0,0::3]=grad[:,:,0];b[:,1,1::3]=grad[:,:,1];b[:,2,2::3]=grad[:,:,2]
    b[:,3,1::3]=grad[:,:,2];b[:,3,2::3]=grad[:,:,1]
    b[:,4,0::3]=grad[:,:,2];b[:,4,2::3]=grad[:,:,0]
    b[:,5,0::3]=grad[:,:,1];b[:,5,1::3]=grad[:,:,0]
    measure=w*det
    mr=(np.sqrt(layer['rho']*measure)[:,None,None]*np.einsum('pi,ab->paib',s,np.eye(3)).reshape(-1,3,81)).reshape(-1,81)
    kr=(np.sqrt(layer['young']*measure)[:,None,None]*np.einsum('ab,pbc->pac',common.chol,b)).reshape(-1,81)
    return dict(s=s,b=b,x=s@q.reshape(27,3),mass_rows=mr,stiffness_rows=kr,
                layer_id=layer['id'],sampled_det_min=float(det.min()))


def face_force(q,axis,sign,zlo,zhi,young,order=9):
    p,w=leggauss(order);coordinates=[];weights=[]
    center=(zlo+zhi)/2;half=(zhi-zlo)/2
    for i,j in itertools.product(range(order),repeat=2):
        if axis==0:coordinates.append([sign,p[i],p[j]])
        elif axis==1:coordinates.append([p[i],sign,p[j]])
        else:coordinates.append([p[i],p[j],sign])
        weights.append(w[i]*w[j])
    s,g=shape_tables(coordinates,center,half)
    jac=np.einsum('ia,pib->pab',q.reshape(27,3),g);x=s@q.reshape(27,3)
    area=sign*np.cross(jac[:,:,(axis+1)%3],jac[:,:,(axis+2)%3])
    traction=np.zeros_like(area)
    traction[:,0]=-young*KAPPA*x[:,2]*area[:,0]*np.array(weights)
    return (s.T@traction).reshape(81)


def local_system(q,layers,internal_sign):
    cert=geometry_certificate(q,min(l['zlo'] for l in layers),max(l['zhi'] for l in layers))
    if not cert['exact_positive']:raise ValueError('reference_domain_not_positive_certified')
    records=[evaluate_layer(q,l) for l in layers]
    mr=np.vstack([r['mass_rows'] for r in records]);kr=np.vstack([r['stiffness_rows'] for r in records])
    faces={(axis,sign):np.zeros(81) for axis in range(3) for sign in [-1,1]}
    for axis in [0,1]:
        for sign in [-1,1]:
            faces[axis,sign]=sum((face_force(q,axis,sign,l['zlo'],l['zhi'],l['young']) for l in layers),np.zeros(81))
    for sign,layer in [(-1,layers[0]),(1,layers[-1])]:
        faces[2,sign]=face_force(q,2,sign,min(l['zlo'] for l in layers),max(l['zhi'] for l in layers),layer['young'])
    internal=faces[0,internal_sign]
    external=sum((f for key,f in faces.items() if key!=(0,internal_sign)),np.zeros(81))
    return dict(q0=q,k=kr.T@kr,m=mr.T@mr,kr=kr,mr=mr,fields=records,
                external=external,internal=internal,all_force=external+internal,certificate=cert,layers=layers)


def natural_scale(h,nodes):
    return np.tile(np.repeat([1.,1./h,2./h**2],3),nodes)


def trace_metrics(left,right,pl,pr,alpha,beta,h):
    sign=-1 if alpha<0 else 1;w=natural_scale(h,15)
    geometry_error=0.;operator_error=0.
    for eta,zeta in itertools.product([-1.,0.,1.],repeat=2):
        z=h*zeta/2
        sl=np.kron(scalar_shape(1,eta,z)[0].reshape(1,27),np.eye(3))
        sr=np.kron(scalar_shape(-sign,eta,(z-beta)/alpha)[0].reshape(1,27),np.eye(3))
        geometry_error=max(geometry_error,float(np.max(np.abs(sl@left['q0']-sr@right['q0']))))
        operator_error=max(operator_error,float(np.max(np.abs((sl@pl-sr@pr)*w[None,:]))))
    return dict(geometry_absolute=geometry_error,scaled_operator_absolute=operator_error)


def physical_references(h):
    layers=canonical_layers(h)
    mass=sum(l['rho']*(l['zhi']-l['zlo']) for l in layers)
    center=sum(l['rho']*(l['zhi']**2-l['zlo']**2)/2 for l in layers)/mass
    energy=.5*KAPPA**2*sum(l['young']*(l['zhi']**3-l['zlo']**3)/3 for l in layers)
    return dict(mass=mass,centroid=np.array([0.,0.,center]),energy=energy)


def global_rigid(q0,center):
    result=np.zeros((NDOF,6))
    for i,(r,d,c) in enumerate(q0.reshape(15,3,3)):
        result[9*i:9*i+3,:3]=np.eye(3)
        for j,omega in enumerate(np.eye(3)):
            for kind,value in enumerate([r-center,d,c]):result[9*i+3*kind:9*i+3*kind+3,3+j]=np.cross(omega,value)
    return result


def exact_raw():
    result=np.zeros((15,3,3))
    for i,(r,_,_) in enumerate(global_geometry().reshape(15,3,3)):
        x,y,_=r
        result[i,0,2]=.5*KAPPA*(x*x-common.nu*y*y)
        result[i,1]=[-KAPPA*x,common.nu*KAPPA*y,0]
        result[i,2,2]=common.nu*KAPPA
    return result.reshape(NDOF)


def solve_global(k,m,f,systems,maps,h,deadline):
    refs=physical_references(h);q0=global_geometry()
    rigid=global_rigid(q0,refs['centroid']);c=rigid.T@m;scale=1/np.sqrt(np.diag(k))
    cs=c*scale[None,:];cs/=np.linalg.norm(cs,axis=1)[:,None]
    _,singular,vh=np.linalg.svd(cs,full_matrices=True)
    constraint_rank=int(np.count_nonzero(singular>1e-10*singular[0]));z=vh[6:].T
    reduced=z.T@(scale[:,None]*k*scale[None,:])@z
    rs=1/np.sqrt(np.diag(reduced));a=rs[:,None]*reduced*rs[None,:];b=rs*(z.T@(scale*f))
    y,solver=curved.high_precision_solve(a,b,deadline);q=(scale[:,None]*z*rs[None,:])@y
    gauge=np.linalg.solve(c@rigid,c@exact_raw());exact=exact_raw()-rigid@gauge
    uerror=0.;uref=0.;eperror=0.;epref=0.;point_values=[]
    for system,p in zip(systems,maps):
        local=p@q
        # Independent reference sampling uses layer order retained in this chart.
        for layer in system['layers']:
            observation=evaluate_layer(system['q0'],layer,9)
            expected_u,expected_eps=curved.exact_fields(observation['x'],False,np.eye(3),np.zeros(3))
            expected_u-=gauge[:3]+np.cross(gauge[3:],observation['x']-refs['centroid'])
            u=observation['s']@local.reshape(27,3);eps=np.einsum('pab,b->pa',observation['b'],local)
            uerror=max(uerror,float(np.max(np.abs(u-expected_u))));uref=max(uref,float(np.max(np.abs(expected_u))))
            eperror=max(eperror,float(np.max(np.abs(eps-expected_eps))));epref=max(epref,float(np.max(np.abs(expected_eps))))
    # A common physical sampling order is used for comparing reparameterizations.
    for e in range(2):
        canonical_local=gather(e)@q
        for eta,zeta in itertools.product(np.linspace(-1,1,5),repeat=2):
            for xi in [-1.,0.,1.]:
                s=scalar_shape(xi,eta,h*zeta/2)[0]
                point_values.extend(s@canonical_local.reshape(27,3))
    energy=float(.5*q@k@q)
    mass=float(rigid[:,0]@m@rigid[:,0])
    centroid=rigid[:,:3].T@m@q0/mass
    fullres=float(np.linalg.norm(scale*(k@q-f))/np.linalg.norm(scale*f))
    constraint_error=float(np.linalg.norm(cs@(q/scale))/np.linalg.norm(q/scale))
    result=dict(solver=solver,constraint_rank=constraint_rank,constraint_relative=constraint_error,
        full_projected_residual=fullres,displacement_relative=uerror/uref,strain_relative=eperror/epref,
        energy=energy,energy_relative=abs(energy/refs['energy']-1),mass=mass,mass_relative=abs(mass/refs['mass']-1),
        centroid=centroid.tolist(),centroid_scaled_error=float(np.max(np.abs(centroid-refs['centroid']))/h),
        dof_errors_by_unit={name:float(np.max(np.abs((q-exact).reshape(15,3,3)[:,i,:])))
            for i,name in enumerate(['position','first_thickness_derivative','second_thickness_derivative'])})
    passed=solver['eligible'] and constraint_rank==6 and fullres<=1e-8 and constraint_error<=1e-10
    passed=passed and result['displacement_relative']<=1e-6 and result['strain_relative']<=1e-6 and result['energy_relative']<=1e-6
    passed=passed and result['mass_relative']<=1e-10 and result['centroid_scaled_error']<=1e-10 and np.all(np.isfinite(q))
    result['eligible']=bool(passed)
    return result,q,np.array(point_values)


def run_case(h,label,alpha,beta,left,right_base,baseline):
    start=time.monotonic();t,perm,sign=chart_transform(alpha,beta)
    right=local_system(t@canonical_geometry(1),local_layers(h,alpha,beta),-sign)
    pl=gather(0);pr=t@gather(1);systems=[left,right];maps=[pl,pr]
    k=sum((p.T@s['k']@p for s,p in zip(systems,maps)),np.zeros((NDOF,NDOF)))
    m=sum((p.T@s['m']@p for s,p in zip(systems,maps)),np.zeros_like(k))
    f=sum((p.T@s['external'] for s,p in zip(systems,maps)),np.zeros(NDOF))
    canonical_k=gather(0).T@left['k']@gather(0)+gather(1).T@right_base['k']@gather(1)
    canonical_m=gather(0).T@left['m']@gather(0)+gather(1).T@right_base['m']@gather(1)
    canonical_f=gather(0).T@left['external']+gather(1).T@right_base['external']
    comparisons=dict(local_stiffness=common.matrix_error(right_base['k'],t.T@right['k']@t),
        local_mass=common.matrix_error(right_base['m'],t.T@right['m']@t),
        local_force=common.force_error(right_base['k'],right_base['all_force'],t.T@right['all_force']),
        global_stiffness=common.matrix_error(canonical_k,k),global_mass=common.matrix_error(canonical_m,m),
        global_force=common.force_error(canonical_k,canonical_f,f))
    mr=np.vstack([s['mr']@p for s,p in zip(systems,maps)]);kr=np.vstack([s['kr']@p for s,p in zip(systems,maps)])
    ranks=dict(mass=common.rank(mr/np.sqrt(np.diag(m))[None,:]),stiffness=common.rank(kr/np.sqrt(np.diag(k))[None,:]))
    hblock=np.array([[1,beta,.5*beta*beta],[0,alpha,alpha*beta],[0,0,alpha*alpha]])
    interface_t=np.kron(np.eye(3),np.kron(hblock,np.eye(3)))
    we=natural_scale(h,3);wf=natural_scale(h/abs(alpha),3)
    condition=np.column_stack((-interface_t,np.eye(27)))*np.concatenate([we,wf])[None,:]/wf[:,None]
    interface_ranks=common.rank(condition)
    trace=trace_metrics(left,right,pl,pr,alpha,beta,h)
    workscale=natural_scale(h,9)
    virtual=np.random.default_rng(1034).standard_normal(81)*workscale*1e-5
    work0=float(virtual@right_base['all_force']);work1=float((t@virtual)@right['all_force'])
    workden=float(np.linalg.norm(virtual/workscale)*np.linalg.norm(workscale*right_base['all_force']))
    work_error=abs(work0-work1)/workden
    global_scale=natural_scale(h,15)
    il=pl.T@left['internal'];ir=pr.T@right['internal']
    internal_error=float(np.linalg.norm(global_scale*(il+ir))/max(np.linalg.norm(global_scale*il),np.linalg.norm(global_scale*ir)))
    solve,q,points=solve_global(k,m,f,systems,maps,h,start+120)
    response=None
    if baseline is not None:
        response=dict(displacement=float(np.max(np.abs(points-baseline['points']))/np.max(np.abs(baseline['points']))),
                      energy=float(abs(solve['energy']/baseline['energy']-1)))
    record=dict(h=h,chart=label,alpha=alpha,beta=beta,local_layers=right['layers'],geometry_certificate=right['certificate'],
        comparisons=comparisons,ranks=ranks,interface_constraint_ranks=interface_ranks,trace=trace,
        virtual_work_relative=work_error,internal_covector_cancellation_relative=internal_error,solve=solve,
        response_to_identity=response,seconds=time.monotonic()-start,peak_process_working_bytes=plate.memory_bytes())
    passed=all(v<=1e-9 for kind in ['local_stiffness','local_mass','global_stiffness','global_mass'] for v in comparisons[kind].values())
    passed=passed and comparisons['local_force']<=1e-9 and comparisons['global_force']<=1e-9
    passed=passed and all(v==135 for v in ranks['mass'].values()) and all(v==129 for v in ranks['stiffness'].values()) and all(v==27 for v in interface_ranks.values())
    passed=passed and trace['geometry_absolute']<=1e-12 and trace['scaled_operator_absolute']<=1e-12 and work_error<=1e-10 and internal_error<=1e-10 and solve['eligible']
    passed=passed and (response is None or all(v<=1e-6 for v in response.values()))
    passed=passed and record['seconds']<=120 and record['peak_process_working_bytes']<=4*1024**3
    record['eligible']=bool(passed)
    return record,dict(points=points,energy=solve['energy'])


def negative_controls():
    h=.01;left=local_system(canonical_geometry(0),local_layers(h,1.,0.),1)
    right_base=local_system(canonical_geometry(1),local_layers(h,1.,0.),-1);checks=[]
    try:chart_transform(0.,0.);rejected=False
    except ValueError:rejected=True
    checks.append(dict(name='zero_alpha',expected='reject_noninvertible',detected=rejected))
    t,_,_=chart_transform(-1.,0.,reverse_x=False)
    try:local_system(t@canonical_geometry(1),local_layers(h,-1.,0.),-1);rejected=False
    except ValueError as exc:rejected=str(exc) in ('reference_domain_not_positive_certified','nonpositive_reference_jacobian')
    cert=geometry_certificate(t@canonical_geometry(1),-h/2,h/2)
    checks.append(dict(name='unpaired_reversal',expected='reject_orientation',detected=rejected,det_upper=cert['upper']))
    alpha=-1.5;beta=h/5;t,perm,sign=chart_transform(alpha,beta)
    right=local_system(t@canonical_geometry(1),local_layers(h,alpha,beta),-sign)
    naive=trace_metrics(left,right,gather(0),perm@gather(1),alpha,beta,h)
    checks.append(dict(name='raw_coefficient_sharing',expected='detect_trace_mismatch',detected=naive['scaled_operator_absolute']>1e-4,trace=naive))
    w=natural_scale(h,9);v=np.random.default_rng(1034).standard_normal(81)*w*1e-5
    true_work=float((t@v)@right['all_force']);bad_work=float(v@(t@right['all_force']))
    mismatch=abs(true_work-bad_work)/(np.linalg.norm(v/w)*np.linalg.norm(w*right_base['all_force']))
    checks.append(dict(name='wrong_force_transform',expected='detect_virtual_work_mismatch',detected=bool(mismatch>1e-4),mismatch=float(mismatch)))
    t,_,sign=chart_transform(-1.,0.)
    wrong=local_system(t@canonical_geometry(1),local_layers(h,-1.,0.,wrong_order=True),-sign)
    ke=common.matrix_error(right_base['k'],t.T@wrong['k']@t);me=common.matrix_error(right_base['m'],t.T@wrong['m']@t)
    checks.append(dict(name='reassigned_layer_ids',expected='detect_physics_change',detected=max(ke['maximum'],me['maximum'])>1e-3,stiffness=ke,mass=me))
    alpha=-1.5;beta=h/5;t,_,sign=chart_transform(alpha,beta)
    shifted_q=(t@canonical_geometry(1)).reshape(9,3,3);shifted_q[:,0,2]+=.01*h
    shifted=dict(q0=shifted_q.reshape(81))
    mismatch=trace_metrics(left,shifted,gather(0),t@gather(1),alpha,beta,h)
    checks.append(dict(name='physical_interface_gap',expected='reject_geometry_mismatch',detected=mismatch['geometry_absolute']>1e-8,trace=mismatch))
    return checks


def main():
    cases=[];stopped=None
    for h in [.1,.01]:
        left=local_system(canonical_geometry(0),local_layers(h,1.,0.),1)
        right=local_system(canonical_geometry(1),local_layers(h,1.,0.),-1)
        baseline=None
        for name,alpha,beta in [('IDENTITY',1.,0.),('SCALE',2.,0.),('SHIFT',1.,h/4),('REVERSE',-1.,0.),('COMBINED',-1.5,h/5)]:
            try:record,solution=run_case(h,name,alpha,beta,left,right,baseline)
            except Exception as exc:record=dict(h=h,chart=name,eligible=False,status='exception',error=str(exc));solution=None
            cases.append(record)
            print(json.dumps({key:record.get(key) for key in ['h','chart','eligible','seconds','error']}),flush=True)
            if name=='IDENTITY':baseline=solution
            if record.get('status')=='exception' or record.get('seconds',0)>120 or record.get('peak_process_working_bytes',0)>4*1024**3:
                stopped=record;break
        if stopped:break
    negatives=negative_controls() if stopped is None else []
    passed=len(cases)==10 and len(negatives)==6 and all(c['eligible'] for c in cases) and all(c['detected'] for c in negatives)
    script=Path(__file__).resolve()
    dependencies=[Path(module.__file__).resolve() for module in [common,curved,plate]]
    report=dict(protocol='docs/verification/q9-interface-reparameterization-experiment.md',numpy=np.__version__,scipy=scipy.__version__,
        cases=cases,negative_controls=negatives,all_criteria_pass=bool(passed),stopped=stopped,
        source_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),
        dependencies=[dict(path=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in dependencies],
        limitation='Conforming two-element flat isotropic two-layer manufactured equilibrium; no physical thickness change, anisotropic laminate, nonlinear or dynamic connection evidence')
    output=script.parents[2]/'build/analysis/q9-interface-reparameterization-result.json'
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(all_criteria_pass=passed,normal_cases=len(cases),negative_controls=negatives)),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
