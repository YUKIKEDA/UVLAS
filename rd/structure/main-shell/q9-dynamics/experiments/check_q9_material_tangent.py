"""Constitutive/nonlinear tangent experiment, not a production material.
Protocol: rd/structure/main-shell/q9-dynamics/records/q9-material-tangent-experiment.md, sections 1-8.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import hashlib
import itertools
import json
import time
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
import check_q9_independent_load as common
import check_q9_curved_load as curved
import check_q9_interface_reparameterization as interface
import check_q9_plate_convergence as plate

PAIRS=[(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]
H=.02
E0S=np.array([[.0005,-.0002,.0001,.0003,-.0001,.0002],
             [-.0003,.0004,-.0001,-.0002,.0002,-.00015],
             [.0002,-.0001,.0003,.0001,-.0003,.00025]])
ANGLES=np.deg2rad([0.,35.,-20.])


def rz(t):
    c=np.cos(t);s=np.sin(t)
    return np.array([[c,-s,0],[s,c,0],[0,0,1.]])


def ry(t):
    c=np.cos(t);s=np.sin(t)
    return np.array([[c,0,s],[0,1.,0],[-s,0,c]])


def rx(t):
    c=np.cos(t);s=np.sin(t)
    return np.array([[1.,0,0],[0,c,-s],[0,s,c]])


def engineering(e):
    return np.stack([e[...,i,j]*(1 if i==j else 2) for i,j in PAIRS],axis=-1)


def stress_tensor(s):
    result=np.zeros(s.shape[:-1]+(3,3),dtype=s.dtype)
    for a,(i,j) in enumerate(PAIRS):result[...,i,j]=result[...,j,i]=s[...,a]
    return result


def strain_tensor(e):
    factors=np.array([1,1,1,.5,.5,.5])
    return stress_tensor(e*factors)


def orthotropic(young=(130e9,10e9,8e9),shear=(3.5e9,4.5e9,5e9),nu=(.28,.25,.32)):
    if not np.all(np.isfinite(list(young)+list(shear)+list(nu))) or min(young+shear)<=0:
        raise ValueError('nonpositive_or_nonfinite_modulus')
    compliance=np.diag(1/np.array(young+shear))
    compliance[0,1]=compliance[1,0]=-nu[0]/young[0]
    compliance[0,2]=compliance[2,0]=-nu[1]/young[0]
    compliance[1,2]=compliance[2,1]=-nu[2]/young[1]
    scale=1/np.sqrt(np.diag(compliance))
    try:np.linalg.cholesky(scale[:,None]*compliance*scale[None,:])
    except np.linalg.LinAlgError:raise ValueError('nonpositive_compliance')
    return np.linalg.solve(compliance,np.eye(6))


def isotropic(young,nu):
    mu=young/(2*(1+nu));lam=young*nu/((1+nu)*(1-2*nu))
    c=np.diag([2*mu]*3+[mu]*3);c[:3,:3]+=lam
    return c


def fourth(c):
    result=np.zeros((3,3,3,3))
    for a,(i,j) in enumerate(PAIRS):
        for b,(k,l) in enumerate(PAIRS):
            for ii,jj in {(i,j),(j,i)}:
                for kk,ll in {(k,l),(l,k)}:result[ii,jj,kk,ll]=c[a,b]
    return result


def validate(f,a,e0):
    if not all(np.all(np.isfinite(v)) for v in [f,a,e0]):raise ValueError('nonfinite_input')
    if np.max(np.abs(a.T@a-np.eye(3)))>1e-12 or np.linalg.det(a)<=0:raise ValueError('invalid_material_basis')
    if np.max(np.abs(e0-e0.T))>1e-12:raise ValueError('nonsymmetric_eigenstrain')
    if np.linalg.det(f)<=0:raise ValueError('invalid_current_jacobian')
    try:np.linalg.cholesky(np.eye(3)+2*e0)
    except np.linalg.LinAlgError:raise ValueError('invalid_natural_metric')


def point_engineering(f,a,e0,c):
    e=.5*(f.T@f-np.eye(3));em=a.T@e@a-e0
    ev=engineering(em);sv=c@ev;sm=stress_tensor(sv);s=a@sm@a.T;p=f@s
    return dict(energy=.5*ev@sv,sm=sm,s=s,p=p,cauchy=p@f.T/np.linalg.det(f),elastic=em)


def point_tensor(f,a,e0,c4):
    em=a.T@(.5*(f.T@f-np.eye(3)))@a-e0
    sm=np.einsum('ijkl,kl->ij',c4,em);s=a@sm@a.T
    return .5*np.einsum('ij,ijkl,kl',em,c4,em),f@s


def relative(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),np.finfo(float).tiny))


def point_tests(c):
    a=rz(.37)@ry(-.23);e0=strain_tensor(np.array([.001,-.0003,.0002,.0004,-.0002,.0003]))
    f=curved.ROTATION@np.array([[1.002,.001,0],[0,.999,.0005],[0,0,1.001]])
    validate(f,a,e0)
    actual=point_engineering(f,a,e0,c);w,p=point_tensor(f,a,e0,fourth(c))
    delta=np.array([[.0002,.0003,-.0001],[.0003,-.0004,.00015],[-.0001,.00015,.0001]])
    work0=float(actual['sm'].reshape(-1)@delta.reshape(-1));work1=float(engineering(delta)@np.array([actual['sm'][i,j] for i,j in PAIRS]))
    q=rz(.43)@rx(-.31);rotated=point_engineering(q@f,a,e0,c)
    metric=np.eye(3)+2*a@e0@a.T;f0=np.linalg.cholesky(metric).T
    natural=point_engineering(f0,a,e0,c)
    metric_error=float(np.max(np.abs(f0.T@f0-metric)))
    result=dict(energy_tensor_relative=float(abs(actual['energy']/w-1)),force_tensor_relative=relative(actual['p'],p),
        work_relative=abs(work0-work1)/max(abs(work0),abs(work1)),
        rotated_energy_relative=float(abs(rotated['energy']/actual['energy']-1)),
        rotated_pk1_relative=relative(rotated['p'],q@actual['p']),
        rotated_cauchy_relative=relative(rotated['cauchy'],q@actual['cauchy']@q.T),
        rotated_pk2_relative=relative(rotated['s'],actual['s']),natural_stress_max=float(np.max(np.abs(natural['s']))),
        natural_energy=float(natural['energy']),natural_metric_error=metric_error)
    keys=['energy_tensor_relative','force_tensor_relative','work_relative','rotated_energy_relative','rotated_pk1_relative','rotated_cauchy_relative','rotated_pk2_relative']
    result['eligible']=all(result[k]<=1e-10 for k in keys) and result['natural_stress_max']<=1e-3 and result['natural_energy']<=1e-12
    wrong_e=.5*(f.T@f-np.eye(3));wrong_em=a.T@wrong_e@a-e0
    wrongv=np.array([wrong_em[i,j] for i,j in PAIRS]);wrong_energy=float(.5*wrongv@c@wrongv)
    bad_axis=point_engineering(q@f,q@a,e0,c)
    negatives=[dict(name='missing_engineering_shear_factor',mismatch=abs(wrong_energy/actual['energy']-1)),
        dict(name='updating_reference_material_axes',mismatch=float(abs(bad_axis['energy']/rotated['energy']-1)))]
    return result,negatives


def axial_tests():
    young=70e9;area=1e-4;length=1.;e0=-.001
    stretch=np.array([1.,0.]);elastic=.5*(stretch@stretch-1)-e0
    force=young*area*elastic*stretch;kmat=young*area/length*np.outer(stretch,stretch)
    tangent=kmat+young*area/length*elastic*np.eye(2)
    energy=.5*young*area*length*elastic**2
    natural=np.sqrt(1+2*e0);residual=young*area*natural*(.5*(natural**2-1)-e0)
    errors=dict(stress=abs(young*elastic-70e6),force=abs(force[0]-7000),energy=abs(energy-3.5),
        axial_tangent=abs(tangent[0,0]-7007000),transverse_tangent=abs(tangent[1,1]-7000),natural_residual=abs(residual))
    limits=dict(stress=1e-4,force=1e-7,energy=1e-10,axial_tangent=1e-4,transverse_tangent=1e-4,natural_residual=1e-4)
    result=dict(errors=errors,eligible=all(errors[k]<=limits[k] for k in errors),transverse_tangent=float(tangent[1,1]))
    return result,dict(name='omitted_geometric_tangent',missing_transverse_stiffness=float(tangent[1,1]-kmat[1,1]))


def section_tests():
    c=isotropic(70e9,.3);e=1e-4;h=.02;n=np.zeros(6);moment=np.zeros(6);energy=0.
    for lo,hi,sign in [(-h/2,0,1.),(0,h/2,-1.)]:
        e0=np.diag([e,-.3*e,-.3*e])*sign
        state=point_engineering(np.eye(3),np.eye(3),e0,c)
        sv=np.array([state['s'][i,j] for i,j in PAIRS])
        for t,weight in zip(*leggauss(3)):
            z=(lo+hi)/2+(hi-lo)*t/2;dv=weight*(hi-lo)/2
            n+=sv*dv;moment+=z*sv*dv;energy+=float(state['energy'])*dv
    result=dict(membrane=n.tolist(),moment=moment.tolist(),energy=energy,
        eligible=bool(np.max(np.abs(n))<=1e-7 and abs(moment[0]/700-1)<=1e-10 and np.max(np.abs(moment[1:]))<=1e-7 and abs(energy/7-1)<=1e-10))
    return result,dict(name='averaged_layer_eigenstrain',missing_energy=float(energy),missing_moment=float(moment[0]))


def natural_model_comparison():
    c4=fourth(isotropic(70e9,.3));f=np.diag([1.003,.999,1.001]);e=.5*(f.T@f-np.eye(3));records=[]
    for e0 in [.001,.05]:
        lam=np.sqrt(1+2*e0);add=e-e0*np.eye(3);fel=f/lam;ee=.5*(fel.T@fel-np.eye(3))
        wa=.5*np.einsum('ij,ijkl,kl',add,c4,add)
        wm=lam**3*.5*np.einsum('ij,ijkl,kl',ee,c4,ee)
        ratio=float(wm/wa);expected=float(1/lam)
        records.append(dict(e0=e0,ratio=ratio,expected=expected,eligible=abs(ratio/expected-1)<=1e-10))
    return records


def context(q0,order,c,alpha=1.,beta=0.,world=np.eye(3)):
    sign=-1 if alpha<0 else 1;points,weights=leggauss(order)
    intervals=[(-H/2,-H/6),(-H/6,H/6),(H/6,H/2)]
    layers=[]
    for idx,(lo,hi) in enumerate(intervals):
        zlo,zhi=sorted([(lo-beta)/alpha,(hi-beta)/alpha]);layers.append((zlo,zhi,idx))
    layers.sort();records=[]
    for zlo,zhi,idx in layers:
        triples=list(itertools.product(range(order),repeat=3))
        coords=np.array([[points[i],points[j],points[k]] for i,j,k in triples])
        weights3=np.array([weights[i]*weights[j]*weights[k] for i,j,k in triples])
        s,g=interface.shape_tables(coords,(zlo+zhi)/2,(zhi-zlo)/2)
        jac=np.einsum('ia,pib->pab',q0.reshape(27,3),g);det=np.linalg.det(jac)
        if np.any(det<=0):raise ValueError('invalid_reference_geometry')
        grad=np.linalg.solve(jac.transpose(0,2,1),g.transpose(0,2,1)).transpose(0,2,1)
        axes=np.array([world@rz(ANGLES[idx]+.07*sign*xi-.04*eta)@ry(.10*sign*xi) for xi,eta,_ in coords])
        e0=np.broadcast_to(strain_tensor(E0S[idx]),(len(coords),3,3)).copy()
        records.append((weights3*det,grad,axes,e0))
    return dict(w=np.concatenate([r[0] for r in records]),g=np.concatenate([r[1] for r in records]),
        axes=np.concatenate([r[2] for r in records]),e0=np.concatenate([r[3] for r in records]),c=c,c4=fourth(c))


def tensor_evaluation(q,ctx):
    f=np.einsum('ia,pib->pab',q.reshape(27,3),ctx['g'])
    e=.5*(np.einsum('pji,pjk->pik',f,f)-np.eye(3))
    em=np.einsum('pji,pjk,pkl->pil',ctx['axes'],e,ctx['axes'])-ctx['e0']
    sm=np.einsum('ijkl,pkl->pij',ctx['c4'],em)
    s=np.einsum('pij,pjk,plk->pil',ctx['axes'],sm,ctx['axes'])
    p=np.einsum('pij,pjk->pik',f,s)
    energy=.5*np.einsum('p,pij,ijkl,pkl->',ctx['w'],em,ctx['c4'],em)
    force=np.einsum('p,pij,paj->ia',ctx['w'],ctx['g'],p).reshape(81)
    return energy,force


def analytic_evaluation(q,ctx):
    g=ctx['g'];axes=ctx['axes'];w=ctx['w'];c=ctx['c']
    f=np.einsum('ia,pib->pab',q.reshape(27,3),g)
    determinant=np.linalg.det(f)
    if np.any(determinant<=0):raise ValueError('invalid_current_geometry')
    e=.5*(np.einsum('pji,pjk->pik',f,f)-np.eye(3))
    em=np.einsum('pji,pjk,pkl->pil',axes,e,axes)-ctx['e0']
    ev=engineering(em);sv=ev@c.T;sm=stress_tensor(sv)
    s=np.einsum('pij,pjk,plk->pil',axes,sm,axes)
    fm=np.einsum('pij,pjk->pik',f,axes);gm=np.einsum('pij,pjk->pik',g,axes)
    b=np.zeros((len(w),6,81))
    for a in range(3):
        for row,(i,j) in enumerate(PAIRS):
            entry=fm[:,a,i,None]*gm[:,:,j]
            if i!=j:entry=entry+fm[:,a,j,None]*gm[:,:,i]
            b[:,row,a::3]=entry
    force=np.einsum('p,pai,pa->i',w,b,sv)
    kmat=np.einsum('p,pai,ab,pbj->ij',w,b,c,b,optimize=True)
    kgscalar=np.einsum('p,pia,pab,pjb->ij',w,g,s,g,optimize=True)
    kgeo=np.kron(kgscalar,np.eye(3))
    return dict(energy=float(.5*np.einsum('p,pa,pa->',w,ev,sv)),force=force,k=kmat+kgeo,kmat=kmat,kgeo=kgeo,
        sampled_det_f_min=float(determinant.min()),max_elastic_green_norm=float(np.linalg.norm(em,axis=(1,2)).max()))


def weighted_force_error(a,b,scale):return relative(scale*a,scale*b)


def matrix_error(a,b,scale):
    ref=scale[:,None]*a*scale[None,:];diff=scale[:,None]*(a-b)*scale[None,:]
    return dict(relative=float(np.linalg.norm(diff)/np.linalg.norm(ref)),maximum=float(np.max(np.abs(diff))/np.max(np.abs(ref))))


def response_difference(a,b,scale):
    return dict(energy=float(abs(a['energy']/b['energy']-1)),force=weighted_force_error(a['force'],b['force'],scale),
        tangent=matrix_error(a['k'],b['k'],scale))


def comparison_pass(record,tolerance):
    return record['energy']<=tolerance and record['force']<=tolerance and all(v<=tolerance for v in record['tangent'].values())


def run_element(name,c):
    start=time.monotonic();q0=curved.geometry(name)
    raw=curved.raw_exact_q(q0,False,np.eye(3),np.zeros(3))*(.02/curved.CURVATURE)
    pre=q0.reshape(27,3)@np.diag([1.002,.999,1.001])+raw.reshape(27,3)
    q=(pre@curved.ROTATION.T).reshape(81)
    ctx=context(q0,7,c);highctx=context(q0,9,c)
    state=analytic_evaluation(q,ctx);high=analytic_evaluation(q,highctx)
    scale=interface.natural_scale(H,9)
    quadrature=response_difference(state,high,scale)
    wt,ft=tensor_evaluation(q,ctx)
    independent=dict(energy=float(abs(state['energy']/wt-1)),force=weighted_force_error(state['force'],ft,scale))
    directions=[];rng=np.random.default_rng(20261029)
    for _ in range(6):
        v=rng.standard_normal(81)*scale
        df=np.einsum('ia,pib->pab',v.reshape(27,3),ctx['g'])
        v/=np.max(np.linalg.norm(df,axis=(1,2)))
        complex_energy,complex_force=tensor_evaluation(q.astype(complex)+1j*1e-30*v,ctx)
        gradient_error=abs(float(np.imag(complex_energy)/1e-30)-float(v@state['force']))/(np.linalg.norm(scale*state['force'])*np.linalg.norm(v/scale))
        action=state['k']@v
        complex_error=weighted_force_error(np.imag(complex_force)/1e-30,action,scale)
        differences=[]
        for step in [1e-2,1e-3,1e-4,1e-5,1e-6]:
            plus=tensor_evaluation(q+step*v,ctx)[1];minus=tensor_evaluation(q-step*v,ctx)[1]
            error=weighted_force_error((plus-minus)/(2*step),action,scale)
            differences.append(dict(step=step,error=error))
        orders=[float(np.log10(differences[i]['error']/differences[i+1]['error'])) for i in range(2)]
        order_ok=differences[0]['error']<=1e-8 or all(1.8<=v<=2.2 for v in orders)
        passed=gradient_error<=1e-9 and complex_error<=1e-9 and min(d['error'] for d in differences)<=1e-6 and order_ok
        directions.append(dict(energy_direction_error=float(gradient_error),complex_tangent_error=complex_error,
            central_difference=differences,initial_orders=orders,eligible=bool(passed)))
    r=rz(.43)@rx(-.31);shift=np.array([.3,-.4,.5]);t=np.kron(np.eye(27),r)
    current=analytic_evaluation(curved.transform(q,r,shift),ctx)
    current_pulled=dict(energy=current['energy'],force=t.T@current['force'],k=t.T@current['k']@t)
    current_error=response_difference(current_pulled,state,scale)
    transformed_ctx=context(curved.transform(q0,r,shift),7,c,world=r)
    both=analytic_evaluation(curved.transform(q,r,shift),transformed_ctx)
    both_pulled=dict(energy=both['energy'],force=t.T@both['force'],k=t.T@both['k']@t)
    both_error=response_difference(both_pulled,state,scale)
    omega=np.array([[0,-3,2],[3,0,-1],[-2,1,0.]])/np.sqrt(14)
    v=(q.reshape(27,3)@omega.T).reshape(81);acc=(v.reshape(27,3)@omega.T).reshape(81)
    first=float(v@state['force']);second=float(v@state['k']@v);extra=float(acc@state['force'])
    first_error=abs(first)/(np.linalg.norm(v/scale)*np.linalg.norm(scale*state['force']))
    second_error=abs(second+extra)/(abs(second)+abs(extra))
    alpha=-1.5;beta=H/5;chart,_,_=interface.chart_transform(alpha,beta)
    otherctx=context(chart@q0,7,c,alpha=alpha,beta=beta)
    chart_state=analytic_evaluation(chart@q,otherctx)
    chart_pulled=dict(energy=chart_state['energy'],force=chart.T@chart_state['force'],k=chart.T@chart_state['k']@chart)
    chart_error=response_difference(chart_pulled,state,scale)
    cert0=curved.certificate(q0,H,np.eye(3),np.zeros(3));cert1=curved.certificate(q,H,np.eye(3),np.zeros(3))
    record=dict(shape=name,quadrature=quadrature,tensor_path=independent,directions=directions,
        current_rotation=current_error,whole_coordinate_rotation=both_error,thickness_chart=chart_error,
        rigid_rotation_work=dict(first_relative=first_error,second_relative=second_error),
        energy=state['energy'],sampled_det_f_min=state['sampled_det_f_min'],max_elastic_green_norm=state['max_elastic_green_norm'],
        reference_certificate=cert0,current_certificate=cert1,seconds=time.monotonic()-start,peak_process_working_bytes=plate.memory_bytes())
    passed=comparison_pass(quadrature,1e-7) and all(v<=1e-10 for v in independent.values()) and all(d['eligible'] for d in directions)
    passed=passed and all(comparison_pass(v,1e-9) for v in [current_error,both_error,chart_error]) and first_error<=1e-9 and second_error<=1e-9
    passed=passed and cert0['exact_positive'] and cert1['exact_positive'] and record['seconds']<=120 and record['peak_process_working_bytes']<=4*1024**3
    record['eligible']=bool(passed)
    return record


def invalid_inputs(c):
    f=np.eye(3);a=np.eye(3);e0=np.zeros((3,3));records=[]
    checks=[('negative_modulus',lambda:orthotropic(young=(-1.,10e9,8e9))),
            ('nonpositive_compliance',lambda:orthotropic(nu=(10.,.25,.32))),
            ('nonorthogonal_material_basis',lambda:validate(f,np.diag([1.,1.,2.]),e0)),
            ('invalid_current_jacobian',lambda:validate(np.diag([-1.,1.,1.]),a,e0)),
            ('invalid_natural_metric',lambda:validate(f,a,np.diag([-.6,0,0])))]
    for name,operation in checks:
        try:operation();detected=False
        except ValueError:detected=True
        records.append(dict(name=name,detected=detected))
    return records


def main():
    c=orthotropic();point,material_negatives=point_tests(c);axial,geo_negative=axial_tests();section,section_negative=section_tests()
    natural=natural_model_comparison();cases=[];stopped=None
    for name in ['RECT','NORMAL-C']:
        try:record=run_element(name,c)
        except Exception as exc:record=dict(shape=name,eligible=False,error=str(exc))
        cases.append(record);print(json.dumps({k:record.get(k) for k in ['shape','eligible','seconds','error']}),flush=True)
        if 'error' in record or record.get('seconds',0)>120 or record.get('peak_process_working_bytes',0)>4*1024**3:stopped=record;break
    negatives=invalid_inputs(c)
    for n in material_negatives:n['detected']=bool(n['mismatch']>1e-4)
    geo_negative['detected']=bool(abs(geo_negative['missing_transverse_stiffness']-7000)<=1e-4)
    section_negative['detected']=bool(abs(section_negative['missing_energy']-7)<=1e-10 and abs(section_negative['missing_moment']-700)<=1e-7)
    negatives+=material_negatives+[geo_negative,section_negative]
    passed=point['eligible'] and axial['eligible'] and section['eligible'] and all(r['eligible'] for r in natural) and len(cases)==2 and all(r['eligible'] for r in cases) and all(r['detected'] for r in negatives)
    source=Path(__file__).resolve();modules=[common,curved,interface,plate]
    result=dict(protocol='rd/structure/main-shell/q9-dynamics/records/q9-material-tangent-experiment.md',numpy=np.__version__,point=point,
        axial=axial,section=section,natural_model_comparison=natural,cases=cases,negative_controls=negatives,
        all_criteria_pass=bool(passed),stopped=stopped,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependencies=[dict(path=Path(m.__file__).name,sha256=hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest()) for m in modules],
        limitation='Synthetic material/element evaluation and derivatives, not physical material validation, equilibrium, eigenvalues, transient or restart')
    path=rd_paths.ROOT/'build/analysis/q9-material-tangent-result.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(all_criteria_pass=passed,point=point,negative_controls=negatives)),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
