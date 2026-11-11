"""Curved/distorted single-element load experiment, not product code.
Protocol: docs/verification/q9-curved-load-experiment.md, sections 1-8.
"""
import hashlib
import itertools
import json
import math
import time
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import lu_factor, lu_solve
import check_q9_independent_load as common
import check_q9_plate_convergence as plate

ZERO=(0,0,0)
CURVATURE=.001
STRAIN=np.array([[.001,.0002,-.0001],[.0002,-.0004,.00015],[-.0001,.00015,.0003]])


def add(a,b):
    result=dict(a)
    for key,value in b.items():
        result[key]=result.get(key,F(0))+value
        if result[key]==0:del result[key]
    return result


def multiply(a,b):
    result={}
    for i,x in a.items():
        for j,y in b.items():
            key=tuple(v+w for v,w in zip(i,j))
            result[key]=result.get(key,F(0))+x*y
    return {k:v for k,v in result.items() if v}


def times(a,v):return {k:x*v for k,x in a.items() if x*v}


def derivative(p,axis):
    result={}
    for key,value in p.items():
        if key[axis]:
            new=list(key);new[axis]-=1
            result[tuple(new)]=value*key[axis]
    return result


def integrate(p):
    result=F(0)
    for key,value in p.items():
        if any(i%2 for i in key):continue
        for i in key:value*=F(2,i+1)
        result+=value
    return result


def geometry_polynomials(q,h):
    lag={-1:[F(0),-F(1,2),F(1,2)],0:[F(1),F(0),-F(1)],1:[F(0),F(1,2),F(1,2)]}
    result=[{}, {}, {}];half=F.from_float(h)/2
    q=q.reshape(9,3,3)
    for node,(xi,eta) in enumerate(common.nodes):
        for kind in range(3):
            shape={(i,j,kind):x*y*half**kind/(2 if kind==2 else 1)
                for i,x in enumerate(lag[int(xi)]) for j,y in enumerate(lag[int(eta)]) if x*y}
            for axis in range(3):
                result[axis]=add(result[axis],times(shape,F.from_float(float(q[node,kind,axis]))))
    return result


def bernstein_bounds(p):
    degree=tuple(max((k[i] for k in p),default=0) for i in range(3))
    shifted={}
    # xi=2t-1, exact power expansion on [0,1]^3.
    for key,value in p.items():
        for j in itertools.product(*(range(k+1) for k in key)):
            factor=math.prod(math.comb(k,l)*2**l*(-1)**(k-l) for k,l in zip(key,j))
            shifted[j]=shifted.get(j,F(0))+value*factor
    coefficients=[]
    for i in itertools.product(*(range(d+1) for d in degree)):
        value=F(0)
        for j,a in shifted.items():
            if a and all(l<=k for l,k in zip(j,i)):
                factor=math.prod(F(math.comb(k,l),math.comb(d,l)) for k,l,d in zip(i,j,degree))
                value+=a*factor
        coefficients.append(value)
    lower=min(coefficients);upper=max(coefficients)
    return dict(degree=list(degree),coefficient_count=len(coefficients),
        exact_positive=bool(lower>0),lower=math.nextafter(float(lower),-math.inf),
        upper=math.nextafter(float(upper),math.inf),
        lower_exact_fraction=str(lower),upper_exact_fraction=str(upper))


def certificate(q,h,rotation,translation):
    polys=geometry_polynomials(q,h)
    jac=[[derivative(p,a) for a in range(3)] for p in polys]
    det={}
    for perm in itertools.permutations(range(3)):
        inversions=sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        term=multiply(multiply(jac[0][perm[0]],jac[1][perm[1]]),jac[2][perm[2]])
        det=add(det,times(term,F((-1)**inversions)))
    bounds=bernstein_bounds(det);volume=integrate(det)
    center=np.array([float(integrate(multiply(p,det))/volume) for p in polys])
    # Original physical Z after pulling the transformed reference back.
    physical_z={}
    for a in range(3):
        shifted=add(polys[a],{ZERO:-F.from_float(float(translation[a]))})
        physical_z=add(physical_z,times(shifted,F.from_float(float(rotation[a,2]))))
    z2=float(integrate(multiply(multiply(physical_z,physical_z),det)))
    bounds.update(volume=float(volume),centroid=center.tolist(),physical_z_second_moment=z2)
    return bounds


def rotation_matrix():
    axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis)
    x,y,z=axis;cross=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
    return np.eye(3)*np.cos(.7)+(1-np.cos(.7))*np.outer(axis,axis)+np.sin(.7)*cross


ROTATION=rotation_matrix();SHIFT=np.array([.3,-.4,.5])


def geometry(name):
    result=np.zeros((9,3,3))
    for i,(xi,eta) in enumerate(common.nodes):
        s=xi/2;t=eta/2
        r=np.array([s,t,0.]);d=np.array([0.,0.,1.]);c=np.zeros(3)
        if name=='SKEW':r[0]+=.4*t
        elif name=='TRAP':r[0]*=1+.6*t
        elif name in ('SADDLE','NORMAL-C'):
            r[2]=.4*s*t;d=np.array([-.4*t,-.4*s,1.])
            if name=='NORMAL-C':d/=np.linalg.norm(d);c=.2*d
        result[i]=[r,d,c]
    return result.reshape(81)


def transform(q,rotation,translation=None):
    result=q.reshape(9,3,3)@rotation.T
    if translation is not None:result[:,0,:]+=translation
    return result.reshape(81)


def shape_tables(points,h):
    scalar=[];gradient=[]
    for xi,eta,zeta in points:
        n,nx,ny=common.shape(xi,eta);z=h*zeta/2
        scalar.append(np.column_stack((n,z*n,.5*z*z*n)).reshape(27))
        gx=np.column_stack((nx,z*nx,.5*z*z*nx)).reshape(27)
        gy=np.column_stack((ny,z*ny,.5*z*z*ny)).reshape(27)
        gz=np.column_stack((np.zeros(9),n,z*n)).reshape(27)*h/2
        gradient.append(np.column_stack((gx,gy,gz)))
    return np.array(scalar),np.array(gradient)


@lru_cache(None)
def volume_table(order,h):
    p,w=leggauss(order)
    indices=list(itertools.product(range(order),repeat=3))
    points=np.array([[p[i],p[j],p[k]] for i,j,k in indices])
    weights=np.array([w[i]*w[j]*w[k] for i,j,k in indices])
    s,g=shape_tables(points,h)
    return s,g,weights


def evaluate(q,h,order):
    s,g,weights=volume_table(order,h)
    jac=np.einsum('ia,pib->pab',q.reshape(27,3),g)
    det=np.linalg.det(jac)
    if np.any(det<=0):raise ValueError('Nonpositive sampled reference Jacobian')
    grad=np.linalg.solve(jac.transpose(0,2,1),g.transpose(0,2,1)).transpose(0,2,1)
    b=np.zeros((len(s),6,81))
    b[:,0,0::3]=grad[:,:,0];b[:,1,1::3]=grad[:,:,1];b[:,2,2::3]=grad[:,:,2]
    b[:,3,1::3]=grad[:,:,2];b[:,3,2::3]=grad[:,:,1]
    b[:,4,0::3]=grad[:,:,2];b[:,4,2::3]=grad[:,:,0]
    b[:,5,0::3]=grad[:,:,1];b[:,5,1::3]=grad[:,:,0]
    measure=weights*det
    mass_rows=(np.sqrt(measure)[:,None,None]*np.einsum('pi,ab->paib',s,np.eye(3)).reshape(-1,3,81)).reshape(-1,81)
    stiffness_rows=(np.sqrt(measure)[:,None,None]*np.einsum('ab,pbc->pac',common.chol,b)).reshape(-1,81)
    return dict(s=s,b=b,jac=jac,x=s@q.reshape(27,3),mass=mass_rows.T@mass_rows,
        stiffness=stiffness_rows.T@stiffness_rows,mass_rows=mass_rows,stiffness_rows=stiffness_rows,
        sampled_det_min=float(det.min()),sampled_jac_condition_max=float(np.linalg.cond(jac).max()),
        volume=float(measure.sum()))


def field(x,affine):
    if affine:
        u=x@STRAIN.T;gradient=np.broadcast_to(STRAIN,(len(x),3,3)).copy()
    else:
        xx,yy,zz=x.T;k=CURVATURE;nu=common.nu
        u=np.column_stack((-k*xx*zz,nu*k*yy*zz,.5*k*(xx*xx-nu*yy*yy+nu*zz*zz)))
        gradient=np.zeros((len(x),3,3))
        gradient[:,0,0]=-k*zz;gradient[:,0,2]=-k*xx
        gradient[:,1,1]=nu*k*zz;gradient[:,1,2]=nu*k*yy
        gradient[:,2,0]=k*xx;gradient[:,2,1]=-nu*k*yy;gradient[:,2,2]=nu*k*zz
    return u,gradient


def exact_fields(x,affine,rotation,translation):
    original=(x-translation)@rotation
    u,gradient=field(original,affine)
    eps=.5*(gradient+gradient.transpose(0,2,1))
    eps=np.einsum('ab,pbc,dc->pad',rotation,eps,rotation)
    eng=np.column_stack((eps[:,0,0],eps[:,1,1],eps[:,2,2],2*eps[:,1,2],2*eps[:,0,2],2*eps[:,0,1]))
    return u@rotation.T,eng


def raw_exact_q(q,affine,rotation,translation):
    nodes=q.reshape(9,3,3);result=np.zeros_like(nodes)
    hessian=np.zeros((3,3,3))
    if not affine:
        k=CURVATURE;nu=common.nu
        hessian[0,0,2]=hessian[0,2,0]=-k
        hessian[1,1,2]=hessian[1,2,1]=nu*k
        hessian[2]=np.diag([k,-nu*k,nu*k])
    for i,(r,d,c) in enumerate(nodes):
        original=(r-translation)@rotation
        u,grad=field(original[None,:],affine)
        db=d@rotation;cb=c@rotation
        result[i,0]=rotation@u[0]
        result[i,1]=rotation@grad[0]@db
        result[i,2]=rotation@(grad[0]@cb+np.einsum('aij,i,j->a',hessian,db,db))
    return result.reshape(81)


def force(q,h,order,affine,rotation,translation):
    p,w=leggauss(order);points=[];weights=[];axes=[];signs=[]
    for axis in range(3):
        others=[i for i in range(3) if i!=axis]
        for sign in [-1,1]:
            for i,j in itertools.product(range(order),repeat=2):
                point=np.zeros(3);point[axis]=sign;point[others]=[p[i],p[j]]
                points.append(point);weights.append(w[i]*w[j]);axes.append(axis);signs.append(sign)
    s,g=shape_tables(points,h)
    jac=np.einsum('ia,pib->pab',q.reshape(27,3),g)
    x=s@q.reshape(27,3);original=(x-translation)@rotation
    if affine:
        stress=common.lam*np.trace(STRAIN)*np.eye(3)+2*common.mu*STRAIN
        stress=np.broadcast_to(stress,(len(x),3,3))
    else:
        stress=np.zeros((len(x),3,3));stress[:,0,0]=-common.young*CURVATURE*original[:,2]
    stress=np.einsum('ab,pbc,dc->pad',rotation,stress,rotation)
    area=np.array([sign*np.cross(j[:,(axis+1)%3],j[:,(axis+2)%3]) for j,axis,sign in zip(jac,axes,signs)])
    traction=np.einsum('pab,pb->pa',stress,area)*np.array(weights)[:,None]
    return (s.T@traction).reshape(81)


def rigid_modes(q,center):
    result=np.zeros((81,6))
    for i,(r,d,c) in enumerate(q.reshape(9,3,3)):
        result[9*i:9*i+3,:3]=np.eye(3)
        for j,omega in enumerate(np.eye(3)):
            for kind,value in enumerate([r-center,d,c]):
                result[9*i+3*kind:9*i+3*kind+3,3+j]=np.cross(omega,value)
    return result


def high_precision_solve(a,b,deadline):
    factor=lu_factor(a)
    with localcontext() as ctx:
        ctx.prec=80
        ad=[[Decimal.from_float(float(v)) for v in row] for row in a]
        bd=[Decimal.from_float(float(v)) for v in b]
        y=[Decimal.from_float(float(v)) for v in lu_solve(factor,b)];history=[]
        for step in range(6):
            if time.monotonic()>deadline:raise TimeoutError('120 second case budget')
            residual=np.array([float(bi-sum((v*z for v,z in zip(row,y)),Decimal(0))) for row,bi in zip(ad,bd)])
            relative=float(np.linalg.norm(residual)/np.linalg.norm(b))
            history.append(dict(corrections=step,relative_residual=relative))
            if relative<=1e-10:break
            if step<5:y=[v+Decimal.from_float(float(d)) for v,d in zip(y,lu_solve(factor,residual))]
        digest=hashlib.sha256('\n'.join(str(v) for v in y).encode('ascii')).hexdigest()
        projected=np.array([float(v) for v in y])
    return projected,dict(history=history,eligible=bool(relative<=1e-10),retained_scaled_solution_sha256=digest,
        projected_residual=float(np.linalg.norm(b-a@projected)/np.linalg.norm(b)))


def solve_case(name,h,rotation,translation):
    start=time.monotonic();q0=transform(geometry(name),rotation,translation);affine=name=='NORMAL-C'
    cert=certificate(q0,h,rotation,translation)
    record=dict(shape=name,h=h,transformed=bool(np.linalg.norm(translation)>0),geometry_certificate=cert)
    if not cert['exact_positive']:return dict(**record,eligible=False,status='geometry_not_certified'),None
    low=evaluate(q0,h,7);high=evaluate(q0,h,9)
    f=force(q0,h,9,affine,rotation,translation);f11=force(q0,h,11,affine,rotation,translation)
    k=low['stiffness'];m=low['mass'];scale=1/np.sqrt(np.diag(k))
    errors=dict(stiffness=common.matrix_error(k,high['stiffness']),mass=common.matrix_error(m,high['mass']),
        force=common.force_error(k,f,f11))
    record.update(quadrature_difference=errors,sampled_det_min=high['sampled_det_min'],
        sampled_jac_condition_max=high['sampled_jac_condition_max'],volume_relative=abs(high['volume']/cert['volume']-1))
    if any(v>1e-8 for kind in ['stiffness','mass'] for v in errors[kind].values()) or errors['force']>1e-8:
        return dict(**record,eligible=False,status='quadrature_failed'),None
    rigid=rigid_modes(q0,np.array(cert['centroid']));constraint=rigid.T@m
    constraint_scaled=constraint*scale[None,:];constraint_scaled/=np.linalg.norm(constraint_scaled,axis=1)[:,None]
    _,singular,vh=np.linalg.svd(constraint_scaled,full_matrices=True)
    rank=int(np.count_nonzero(singular>1e-10*singular[0]))
    z=vh[6:].T
    reduced=z.T@(scale[:,None]*k*scale[None,:])@z
    reduced_scale=1/np.sqrt(np.diag(reduced))
    a=reduced_scale[:,None]*reduced*reduced_scale[None,:]
    rhs=reduced_scale*(z.T@(scale*f))
    y,solve=high_precision_solve(a,rhs,start+120)
    lift=scale[:,None]*z*reduced_scale[None,:];q=lift@y
    exact_raw=raw_exact_q(q0,affine,rotation,translation)
    gauge=np.linalg.solve(constraint@rigid,constraint@exact_raw)
    exact=exact_raw-rigid@gauge
    uref,eref=exact_fields(high['x'],affine,rotation,translation)
    raw_u=high['s']@exact_raw.reshape(27,3)
    raw_eps=np.einsum('pab,b->pa',high['b'],exact_raw)
    repro_u=float(np.max(np.abs(raw_u-uref))/np.max(np.abs(uref)))
    repro_eps=float(np.max(np.abs(raw_eps-eref))/np.max(np.abs(eref)))
    center=np.array(cert['centroid'])
    adjusted_u=uref-gauge[:3]-np.cross(gauge[3:],high['x']-center)
    actual_u=high['s']@q.reshape(27,3);actual_eps=np.einsum('pab,b->pa',high['b'],q)
    displacement=float(np.max(np.abs(actual_u-adjusted_u))/np.max(np.abs(adjusted_u)))
    strain=float(np.max(np.abs(actual_eps-eref))/np.max(np.abs(eref)))
    if affine:
        eng=np.array([STRAIN[0,0],STRAIN[1,1],STRAIN[2,2],2*STRAIN[1,2],2*STRAIN[0,2],2*STRAIN[0,1]])
        reference_energy=float(.5*eng@common.elastic@eng*cert['volume'])
    else:reference_energy=.5*common.young*CURVATURE**2*cert['physical_z_second_moment']
    energy=float(.5*q@k@q);energy_error=abs(energy/reference_energy-1)
    residual=float(np.linalg.norm(scale*(k@q-f))/np.linalg.norm(scale*f))
    rigid_work=float(np.max(np.abs(rigid.T@f)/(np.linalg.norm(rigid/scale[:,None],axis=0)*np.linalg.norm(scale*f))))
    constraint_error=float(np.linalg.norm(constraint_scaled@(q/scale))/np.linalg.norm(q/scale))
    ranks=dict(mass=common.rank(low['mass_rows']/np.sqrt(np.diag(m))[None,:]),
        stiffness=common.rank(low['stiffness_rows']*scale[None,:]))
    finite_q=transform(q0,ROTATION,SHIFT)
    _,g,_=volume_table(9,h)
    finite_jac=np.einsum('ia,pib->pab',finite_q.reshape(27,3),g)
    fdeform=np.linalg.solve(high['jac'].transpose(0,2,1),finite_jac.transpose(0,2,1)).transpose(0,2,1)
    green=.5*(np.einsum('pji,pjk->pik',fdeform,fdeform)-np.eye(3))
    rigid_f=float(np.max(np.abs(fdeform-ROTATION)));rigid_e=float(np.max(np.abs(green)))
    groups={}
    for kind,label in enumerate(['position','first_thickness_derivative','second_thickness_derivative']):
        delta=(q-exact).reshape(9,3,3)[:,kind,:]
        groups[label]=dict(max_abs_error=float(np.max(np.abs(delta))),reference_max=float(np.max(np.abs(exact.reshape(9,3,3)[:,kind,:]))))
    record.update(status='evaluated',solve=solve,constraint_rank=rank,constraint_relative=constraint_error,
        ranks=ranks,rigid_load_work=rigid_work,full_projected_residual=residual,
        reference_reproduction=dict(displacement=repro_u,strain=repro_eps),
        displacement_relative=displacement,strain_relative=strain,energy_relative=energy_error,energy=energy,
        unit_separated_dof_errors=groups,finite_rigid_kinematics=dict(f_error=rigid_f,green_error=rigid_e),
        seconds=time.monotonic()-start,peak_process_working_bytes=plate.memory_bytes())
    passed=rank==6 and all(v==81 for v in ranks['mass'].values()) and all(v==75 for v in ranks['stiffness'].values()) and solve['eligible']
    passed=passed and rigid_work<=1e-10 and residual<=1e-8 and constraint_error<=1e-10 and record['volume_relative']<=1e-10
    passed=passed and repro_u<=1e-10 and repro_eps<=1e-10 and displacement<=1e-6 and strain<=1e-6 and energy_error<=1e-6
    passed=passed and rigid_f<=1e-12 and rigid_e<=1e-12 and record['seconds']<=120 and record['peak_process_working_bytes']<=4*1024**3 and np.all(np.isfinite(q))
    record['eligible']=bool(passed)
    return record,dict(k=k,m=m,f=f,q=q,energy=energy,points=high['s']@q.reshape(27,3))


def main():
    cases=[];comparisons=[];stop=None
    identity=np.eye(3);zero=np.zeros(3);t=np.kron(np.eye(27),ROTATION)
    for name in ['RECT','SKEW','TRAP','SADDLE','NORMAL-C']:
        for h in [.1,.01]:
            pair=[]
            for rotation,translation in [(identity,zero),(ROTATION,SHIFT)]:
                try:record,data=solve_case(name,h,rotation,translation)
                except Exception as exc:record=dict(shape=name,h=h,eligible=False,status='exception',error=str(exc));data=None
                cases.append(record);pair.append(data)
                print(json.dumps({k:record.get(k) for k in ['shape','h','transformed','status','eligible','seconds','error']}),flush=True)
                if record['status']=='exception' or record.get('seconds',0)>120 or record.get('peak_process_working_bytes',0)>4*1024**3:stop=record;break
            if len(pair)==2 and all(v is not None for v in pair):
                base,moved=pair
                errors=dict(stiffness=common.matrix_error(base['k'],t.T@moved['k']@t),mass=common.matrix_error(base['m'],t.T@moved['m']@t),
                    force=common.force_error(base['k'],base['f'],t.T@moved['f']),
                    displacement=float(np.max(np.abs(moved['points']@ROTATION-base['points']))/np.max(np.abs(base['points']))),
                    energy=float(abs(moved['energy']/base['energy']-1)))
                passed=all(v<=1e-9 for kind in ['stiffness','mass'] for v in errors[kind].values()) and errors['force']<=1e-9 and errors['displacement']<=1e-6 and errors['energy']<=1e-8
                comparisons.append(dict(shape=name,h=h,errors=errors,eligible=bool(passed)))
            if stop:break
        if stop:break
    passed=len(cases)==20 and len(comparisons)==10 and all(c['eligible'] for c in cases+comparisons)
    report=dict(protocol='docs/verification/q9-curved-load-experiment.md',cases=cases,coordinate_comparisons=comparisons,
        rotation=ROTATION.tolist(),translation=SHIFT.tolist(),all_criteria_pass=passed,stopped=stop,
        limitation='Single element manufactured linear fields; no generic curved bending convergence or nonlinear tangent evidence')
    path=Path(__file__).resolve().parents[2]/'build/analysis/q9-curved-load-result.json'
    path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(all_criteria_pass=passed,cases=len(cases),coordinate_comparisons=len(comparisons))),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
