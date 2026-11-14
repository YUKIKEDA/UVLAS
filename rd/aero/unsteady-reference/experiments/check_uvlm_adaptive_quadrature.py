"""Predeclared off-surface adaptive W-P1 design experiment; no product API."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib
import json
import math
import sys
import time
import numpy as np
import scipy
import check_uvlm_wp1_three_dimensional as base

ROOT=rd_paths.ROOT
START=time.monotonic()
RECORDS=[]

def segment_distance(p,a,b):
    edge=b-a
    t=np.clip(np.dot(p-a,edge)/np.dot(edge,edge),0.,1.)
    return float(np.linalg.norm(p-a-t*edge))

def triangle_distance(p,x):
    e,n,_=base.triangle(x)
    r=p-x[0]
    uv=np.linalg.solve(e.T@e,e.T@r)
    if min(uv)>=0 and sum(uv)<=1:
        return float(abs(r@n))
    return min(segment_distance(p,x[i],x[(i+1)%3]) for i in range(3))

def compensated(values):
    return np.array([[math.fsum(v[j,k] for v in values) for k in range(3)] for j in range(2)])

class Limit(Exception):
    pass

def evaluate(patch,point,eps=0.,atol=1e-10,rtol=1e-8,max_depth=20,max_leaves=100000,
             max_evaluations=20000000,self_policy='off_surface_only'):
    leaves={};count=0;deepest=0;last=None;minimum=math.inf;audits=0
    def failure(status):
        return dict(status=status,diagnostic_value=None if last is None else last.tolist(),
            evaluations=count,leaves=len(leaves),max_depth=deepest,audits=audits)
    try:
        patch=np.asarray(patch,dtype=float);point=np.asarray(point,dtype=float)
        try: base.validate(patch,np.array([1.,0.]),point[None,:],eps,self_policy)
        except ValueError as exc:
            msg=str(exc)
            status=('unsupported_self_policy' if 'self' in msg else
                    'contact' if 'contact' in msg else 'degenerate_geometry' if 'degenerate' in msg else 'invalid_input')
            return failure(status)
        length=max(np.linalg.norm(patch[(i+1)%4]-patch[i]) for i in range(4))
        origin=patch[0].copy();x=(patch-origin)/length;p=(point-origin)/length;core=eps/length
        mu=np.array([[1.,0.],[1.,0.],[0.,1.],[0.,1.]])
        def integrate(leaf,order):
            nonlocal count
            kind,xx,coeff,depth,ident=leaf
            required=order if kind=='e' else order**2
            if count+required>max_evaluations: raise Limit('evaluation_limit')
            count+=required
            if kind=='e':
                q,w=base.gauss(order);edge=xx[1]-xx[0]
                y=xx[0]+q[:,None]*edge;r=p-y
                denom=(np.sum(r*r,axis=1)+core**2)**1.5
                strength=(1-q[:,None])*coeff[0]+q[:,None]*coeff[1]
                return np.einsum('q,qb,qc->bc',w,strength,np.cross(edge,r)/denom[:,None])/(4*np.pi)
            y,w,_,_=base.quadrature(xx,order);r=p-y
            kernel=np.cross(coeff[:,None,:],r[None,:,:])/(np.sum(r*r,axis=1)+core**2)[None,:,None]**1.5
            return np.einsum('q,bqc->bc',w,kernel)/(4*np.pi)
        def make(kind,xx,coeff,depth,ident):
            nonlocal deepest,minimum
            deepest=max(deepest,depth)
            leaf=(kind,xx,coeff,depth,ident)
            distance=segment_distance(p,*xx) if kind=='e' else triangle_distance(p,xx)
            minimum=min(minimum,distance)
            h=np.linalg.norm(xx[1]-xx[0]) if kind=='e' else max(np.linalg.norm(xx[i]-xx[(i+1)%3]) for i in range(3))
            ratio=h/np.sqrt(distance**2+core**2)
            low=integrate(leaf,8);high=integrate(leaf,16)
            return dict(leaf=leaf,value=high,error=16*np.linalg.norm(high-low,axis=1),ratio=float(ratio))
        for i in range(4):
            ident=f'e{i}'
            leaves[ident]=make('e',x[[i,(i+1)%4]],mu[[i,(i+1)%4]],0,ident)
        for i,ids in enumerate(base.TRIS):
            xx=x[ids];e,n,_=base.triangle(xx)
            gradient=(e@np.linalg.solve(e.T@e,mu[ids][1:]-mu[ids][0])).T
            vort=np.cross(n,gradient)
            ident=f't{i}';leaves[ident]=make('t',xx,vort,0,ident)
        if len(leaves)>max_leaves: raise Limit('leaf_limit')
        while True:
            if time.monotonic()-START>180: raise Limit('experiment_time_limit')
            ordered=[leaves[k] for k in sorted(leaves)]
            total=compensated([v['value'] for v in ordered]);last=total/length
            target=atol+rtol*np.linalg.norm(total,axis=1)
            error=np.array([math.fsum(v['error'][j] for v in ordered) for j in range(2)])
            magnitude=np.array([math.fsum(np.linalg.norm(v['value'][j]) for v in ordered) for j in range(2)])
            floor=8*np.finfo(float).eps*magnitude
            if np.any(floor>target): raise Limit('roundoff_limited')
            unsafe=[v for v in ordered if v['ratio']>2]
            audit_error=None
            if not unsafe and np.all(error<=target):
                audits+=1
                for v in ordered:
                    audited=integrate(v['leaf'],24)
                    v['audit_error']=np.linalg.norm(audited-v['value'],axis=1)
                audit_error=np.array([math.fsum(v['audit_error'][j] for v in ordered) for j in range(2)])
                if np.all(audit_error<=target/4):
                    ids=[v['leaf'][4] for v in ordered]
                    tree=hashlib.sha256(json.dumps(ids,separators=(',',':')).encode()).hexdigest()
                    return dict(status='estimated_converged',value=last.tolist(),estimated_error=(error/length).tolist(),
                        audit_error=(audit_error/length).tolist(),roundoff_floor=(floor/length).tolist(),
                        target=(target/length).tolist(),evaluations=count,leaves=len(leaves),max_depth=deepest,
                        max_ratio=max(v['ratio'] for v in ordered),minimum_distance=minimum*length,
                        tree_sha256=tree,audits=audits,parent_length=float(length))
            if unsafe:
                chosen=min(unsafe,key=lambda v:(-v['ratio'],v['leaf'][4]))
            else:
                field='audit_error' if audit_error is not None else 'error'
                chosen=min(ordered,key=lambda v:(-float(np.max(v[field]/target)),v['leaf'][4]))
            kind,xx,coeff,depth,ident=chosen['leaf']
            if depth>=max_depth: raise Limit('depth_limit')
            growth=1 if kind=='e' else 3
            if len(leaves)+growth>max_leaves: raise Limit('leaf_limit')
            children=[]
            if kind=='e':
                midpoint=xx.mean(axis=0);mm=coeff.mean(axis=0)
                children=[(np.array([xx[0],midpoint]),np.array([coeff[0],mm])),
                          (np.array([midpoint,xx[1]]),np.array([mm,coeff[1]]))]
            else:
                a,b,c=xx;ab=(a+b)/2;bc=(b+c)/2;ca=(c+a)/2
                children=[(np.array(y),coeff) for y in [[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]]]
            # Build every child before publishing the replacement.
            prepared={f'{ident}.{j}':make(kind,yy,cc,depth+1,f'{ident}.{j}') for j,(yy,cc) in enumerate(children)}
            del leaves[ident];leaves.update(prepared)
    except Limit as exc:
        return failure(str(exc))

def consume(result,coefficients):
    if result['status']!='estimated_converged': raise ValueError('ineligible_velocity')
    return np.asarray(coefficients)@np.asarray(result['value'])

def reference_check(result,reference,length):
    if result['status']!='estimated_converged':return dict(eligible=False,error=None,estimator_eligible=False)
    actual=np.linalg.norm(np.array(result['value'])-reference,axis=1)
    scale=np.maximum(np.linalg.norm(reference,axis=1),1/length)
    bound=np.array(result['estimated_error'])+result['roundoff_floor']+1e-12*scale
    return dict(eligible=bool(np.all(actual<=1e-10/length+1e-8*np.linalg.norm(reference,axis=1))),
        error=actual.tolist(),relative_error=(actual/np.maximum(np.linalg.norm(reference,axis=1),1e-300)).tolist(),
        estimator_eligible=bool(np.all(actual<=bound)))

def main():
    previous_path=ROOT/'rd/aero/kernel/records/uvlm-near-field-screening-result.json'
    previous=json.loads(previous_path.read_text());patch=np.array(previous['patch'])
    for row in previous['cases']:
        p=np.array([1-row['distance'],0.,0.])
        reference=np.array([[0.,0.,b['analytic']] for b in row['basis']])
        result=evaluate(patch,p)
        check=reference_check(result,reference,previous['parent_length'])
        RECORDS.append(dict(group='planar',distance=row['distance'],point=p.tolist(),result=result,comparison=check))
        print(json.dumps(dict(group='planar',distance=row['distance'],status=result['status'],leaves=result['leaves'])),flush=True)
    independent=[]
    for name,xx in base.SHAPES.items():
        length=max(np.linalg.norm(xx[(i+1)%4]-xx[i]) for i in range(4))
        p=np.array([.65,.1,.2])
        for eps in [0.,.08]:
            refs=[];diffs=[]
            for mu in np.eye(2):
                low=base.area_reference(xx,mu,p,eps,order=96)[0]
                high=base.area_reference(xx,mu,p,eps,order=128)[0]
                refs.append(high);diffs.append(float(np.linalg.norm(high-low)))
            reference=np.array(refs);result=evaluate(xx,p,eps)
            check=reference_check(result,reference,length)
            check['reference_converged']=bool(np.all(np.array(diffs)<=(1e-10/length+1e-8*np.linalg.norm(reference,axis=1))/20))
            independent.append(dict(shape=name,eps=eps,reference=reference.tolist(),reference_change=diffs,result=result,comparison=check))
    constant=[]
    for name,xx in base.SHAPES.items():
        length=max(np.linalg.norm(xx[(i+1)%4]-xx[i]) for i in range(4))
        tri=xx[base.TRIS[0]];_,normal,_=base.triangle(tri)
        for kind,anchor in [('face',tri.mean(axis=0)),('edge',(tri[0]+tri[1])/2),('vertex',tri[0])]:
            p=anchor+1e-3*length*normal
            for eps in [0.,.02*length]:
                result=evaluate(xx,p,eps)
                reference=base.core.kernel([p],xx,np.roll(xx,-1,axis=0),eps)[0].sum(axis=0)
                if result['status']=='estimated_converged':
                    error=float(np.linalg.norm(consume(result,[1.,1.])-reference))
                    threshold=2e-10/length+1e-8*np.linalg.norm(reference)
                    eligible=error<=threshold
                else:error=None;eligible=False
                constant.append(dict(shape=name,kind=kind,eps=eps,point=p.tolist(),reference=reference.tolist(),result=result,error=error,eligible=bool(eligible)))
        print(json.dumps(dict(group='constant',shape=name,completed=len(constant))),flush=True)
    p=np.array([.99,0.,0.]);baseline=evaluate(patch,p);again=evaluate(patch,p)
    rz=np.array([[np.cos(.4),-np.sin(.4),0],[np.sin(.4),np.cos(.4),0],[0,0,1]])
    ry=np.array([[np.cos(-.3),0,np.sin(-.3)],[0,1,0],[-np.sin(-.3),0,np.cos(-.3)]])
    rotation=rz@ry;translation=np.array([.3,-.2,.5]);transforms=[]
    for name,xx,pp,ref in [('rotation',patch@rotation.T+translation,p@rotation.T+translation,np.array(baseline['value'])@rotation.T),
        ('small',patch*1e-3,p*1e-3,np.array(baseline['value'])/1e-3),
        ('large',patch*1e3,p*1e3,np.array(baseline['value'])/1e3)]:
        r=evaluate(xx,pp)
        error=None if r['status']!='estimated_converged' else float(np.linalg.norm(np.array(r['value'])-ref)/max(np.linalg.norm(ref),1e-300))
        transforms.append(dict(kind=name,result=r,relative_error=error,eligible=error is not None and error<=1e-8))
    failures=[];out=np.array([12.,34.,56.]);out_hash=hashlib.sha256(out.tobytes()).hexdigest()
    patch_hash=hashlib.sha256(patch.tobytes()).hexdigest()
    for name,kwargs,expected in [('depth',dict(max_depth=0),'depth_limit'),('leaf',dict(max_leaves=6),'leaf_limit'),
        ('evaluation',dict(max_evaluations=1),'evaluation_limit'),
        ('roundoff',dict(atol=0.,rtol=1e-18),'roundoff_limited')]:
        pp=np.array([0.,0.,0.]) if name=='roundoff' else np.array([.999,0.,0.])
        r=evaluate(patch,pp,**kwargs)
        try:out[:]=consume(r,[1.,1.]);blocked=False
        except ValueError:blocked=True
        failures.append(dict(name=name,result=r,expected=expected,eligible=r['status']==expected and blocked))
    for name,xx,pp,eps,kwargs in [
        ('contact_zero',patch,np.array([1.2,0.,0.]),0.,{}),
        ('contact_core',patch,np.array([1.2,0.,0.]),.08,{}),
        ('degenerate',np.zeros((4,3)),p,0.,{}),
        ('nonfinite',patch,np.array([np.nan,0.,0.]),0.,{}),
        ('self',patch,p,0.,dict(self_policy='principal_value')),
        ('negative_core',patch,p,-.1,{})]:
        r=evaluate(xx,pp,eps,**kwargs)
        try:out[:]=consume(r,[1.,1.]);blocked=False
        except ValueError:blocked=True
        failures.append(dict(name=name,result=r,eligible=r['status']!='estimated_converged' and blocked))
    test_tri=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    distances=[dict(point=pp,expected=ref,actual=triangle_distance(np.array(pp),test_tri)) for pp,ref in
        [([.2,.2,.3],.3),([.5,-.2,.3],np.sqrt(.13)),([-.1,-.2,.3],np.sqrt(.14))]]
    complex_weights=np.array([1+.2j,-.4+.7j]);cancel_weights=np.array([1.,-1.])
    linearity=all(np.array_equal(consume(baseline,c),c@np.array(baseline['value'])) for c in [complex_weights,cancel_weights])
    checks=dict(planar_accuracy=all(r['comparison']['eligible'] for r in RECORDS),
        planar_estimator=all(r['comparison']['estimator_eligible'] for r in RECORDS),
        independent_reference_convergence=all(r['comparison']['reference_converged'] for r in independent),
        independent_accuracy=all(r['comparison']['eligible'] for r in independent),
        independent_estimator=all(r['comparison']['estimator_eligible'] for r in independent),
        constant_near_field=all(r['eligible'] for r in constant),
        covariance=all(r['eligible'] for r in transforms),deterministic=baseline==again,
        linearity=linearity,failures=all(r['eligible'] for r in failures),
        buffers_unchanged=hashlib.sha256(out.tobytes()).hexdigest()==out_hash and hashlib.sha256(patch.tobytes()).hexdigest()==patch_hash,
        distances=all(abs(r['actual']-r['expected'])<=1e-13 for r in distances))
    source=Path(__file__).resolve()
    loaded=sorted({Path(m.__file__).resolve() for m in list(sys.modules.values()) if getattr(m,'__file__',None)
        and rd_paths.is_experiment(m.__file__) and Path(m.__file__).resolve()!=source})
    result=dict(protocol='rd/aero/unsteady-reference/records/uvlm-adaptive-quadrature-experiment.md',
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in loaded},
        baseline_result_sha256=hashlib.sha256(previous_path.read_bytes()).hexdigest(),
        planar=RECORDS,independent=independent,constant=constant,transformations=transforms,
        failures=failures,distances=distances,deterministic_baseline=baseline,
        complex_velocity=[[float(v.real),float(v.imag)] for v in consume(baseline,complex_weights)],
        cancellation_velocity=consume(baseline,cancel_weights).tolist(),checks=checks,all_criteria_pass=all(checks.values()),
        seconds=time.monotonic()-START,peak_process_working_bytes=base.core.memory.memory_bytes(),numpy=np.__version__,scipy=scipy.__version__,
        limits='limited off-surface evaluator only; no assembled loads, cache lifecycle, tangent or self induction')
    (ROOT/'rd/aero/unsteady-reference/records/uvlm-adaptive-quadrature-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=base.core.json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],peak_bytes=result['peak_process_working_bytes'])),flush=True)
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':
    try:main()
    except Exception as exc:
        (ROOT/'build/analysis/uvlm-adaptive-quadrature-stopped.json').write_text(json.dumps(dict(error=repr(exc),planar=RECORDS,seconds=time.monotonic()-START),indent=2,default=base.core.json_scalar)+'\n',encoding='utf-8')
        raise
