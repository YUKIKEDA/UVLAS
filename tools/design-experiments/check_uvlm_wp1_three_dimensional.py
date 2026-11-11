"""Predeclared finite 3D W-P1 sheet, core consistency and strip-state trial."""
import hashlib
import json
import sys
import time
from functools import lru_cache
from pathlib import Path
import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
import check_vortex_core_load as core
import check_aerogrid_geometry as geometry
import check_uvlm_trailing_edge as planar

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()
RECORDS=[]
TRIS=[[0,2,1],[0,3,2]]
POINTS=np.array([[.65,.1,.2],[1.2,0.,.45],[1.6,-.9,.25],[1.25,.9,-.4]])
SHAPES=dict(flat=np.array([[1.,-.6,0.],[1.,.6,0.],[1.4,.6,0.],[1.4,-.6,0.]]),
    skew=np.array([[1.,-.6,0.],[1.2,.6,0.],[1.65,.5,0.],[1.4,-.6,0.]]),
    warped=np.array([[1.,-.6,0.],[1.,.6,.03],[1.4,.6,.18],[1.4,-.6,-.08]]))


def budget():
    if time.monotonic()-START>180 or core.memory.memory_bytes()>512*1024**2:
        raise RuntimeError('experiment_budget_exceeded')


@lru_cache(None)
def gauss(n):
    q,w=leggauss(n)
    return (q+1)/2,w/2


def triangle(x):
    e=np.column_stack([x[1]-x[0],x[2]-x[0]])
    cross=np.cross(e[:,0],e[:,1]);jac=np.linalg.norm(cross)
    if jac<=1e-12:raise ValueError('degenerate_facet')
    return e,cross/jac,jac


def validate(quad,values,points,eps,self_policy):
    if self_policy!='off_surface_only':raise ValueError('unsupported_self_evaluation')
    if quad.shape!=(4,3) or np.asarray(values).shape!=(2,):raise ValueError('invalid_patch_shape')
    if eps<0 or not all(np.all(np.isfinite(v)) for v in [quad,values,points,np.asarray(eps)]):raise ValueError('nonfinite_or_negative_core')
    for ids in TRIS:
        x=quad[ids];e,n,_=triangle(x)
        for p in points:
            r=p-x[0];coordinates=np.linalg.solve(e.T@e,e.T@r)
            if abs(r@n)<=1e-11 and min(coordinates)>=-1e-11 and sum(coordinates)<=1+1e-11:
                raise ValueError('nonincident_surface_contact')


def quadrature(x,n):
    _,normal,jac=triangle(x);q,w=gauss(n)
    u,v=np.meshgrid(q,q,indexing='ij');wu,wv=np.meshgrid(w,w,indexing='ij')
    bary=np.stack([(1-u)*(1-v),u,(1-u)*v],axis=-1).reshape(-1,3)
    return bary@x,(wu*wv*(1-u)*jac).reshape(-1),bary,normal


def line(points,a,b,ma,mb,eps,order=48):
    q,w=gauss(order);samples=a[None,:]+q[:,None]*(b-a)[None,:]
    r=points[:,None,:]-samples[None,:,:];den=np.sum(r*r,axis=2)+eps*eps
    values=ma*(1-q)+mb*q
    return np.einsum('q,pqj->pj',w,np.cross(b-a,r)*values[None,:,None]/den[:,:,None]**1.5)/(4*np.pi)


def triangle_primary(points,x,values,eps,line_order=48,surface_order=24,boundary=True,sheet=True):
    e,n,_=triangle(x)
    result=np.zeros((len(points),3),dtype=np.result_type(values,float))
    if boundary:
        for i,j in [(0,2),(2,1),(1,0)]:result+=line(points,x[i],x[j],values[i],values[j],eps,line_order)
    if sheet:
        gradient=e@np.linalg.solve(e.T@e,np.array([values[1]-values[0],values[2]-values[0]]))
        vorticity=np.cross(n,gradient)
        if np.linalg.norm(vorticity)>0:
            y,w,_,_=quadrature(x,surface_order);r=points[:,None,:]-y[None,:,:]
            den=np.sum(r*r,axis=2)+eps*eps
            result+=np.einsum('q,pqj->pj',w,np.cross(vorticity,r)/den[:,:,None]**1.5)/(4*np.pi)
    return result


def primary(quad,values,points=POINTS,eps=0.,self_policy='off_surface_only',sheet=True,line_order=48,surface_order=24):
    quad=np.asarray(quad);points=np.atleast_2d(points);values=np.asarray(values)
    validate(quad,values,points,eps,self_policy)
    mu=np.array([values[0],values[0],values[1],values[1]])
    result=np.zeros((len(points),3),dtype=np.result_type(values,float))
    for ids in TRIS:result+=triangle_primary(points,quad[ids],mu[ids],eps,line_order,surface_order,sheet=sheet)
    return result


def area_reference(quad,values,points=POINTS,eps=0.,order=64,soft_gradient=False):
    quad=np.asarray(quad);points=np.atleast_2d(points);values=np.asarray(values)
    validate(quad,values,points,eps,'off_surface_only')
    mu=np.array([values[0],values[0],values[1],values[1]])
    result=np.zeros((len(points),3),dtype=np.result_type(values,float))
    for ids in TRIS:
        y,w,bary,n=quadrature(quad[ids],order);r=points[:,None,:]-y[None,:,:]
        r2=np.sum(r*r,axis=2);den=r2+eps*eps;rn=np.einsum('pqj,j->pq',r,n)
        coefficient=r2+eps*eps if soft_gradient else r2-2*eps*eps
        kernel=(coefficient[:,:,None]*n[None,None,:]-3*r*rn[:,:,None])/den[:,:,None]**2.5
        result+=np.einsum('q,pqj->pj',w*(bary@mu[ids]),kernel)/(4*np.pi)
    return result


def error(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-15))


def point_eligible(a,b):
    d=np.linalg.norm(a-b,axis=1);ref=np.linalg.norm(b,axis=1)
    return bool(np.all((d<=1e-10)|(d<=1e-8*ref)))


def subdivide(triangles,values):
    x=triangles;a,b,c=x[:,0],x[:,1],x[:,2];ab=(a+b)/2;bc=(b+c)/2;ca=(c+a)/2
    mu=values;ma,mb,mc=mu[:,0],mu[:,1],mu[:,2];mab=(ma+mb)/2;mbc=(mb+mc)/2;mca=(mc+ma)/2
    xx=np.stack([np.stack([a,ab,ca],1),np.stack([ab,b,bc],1),np.stack([ca,bc,c],1),np.stack([ab,bc,ca],1)],axis=1)
    mm=np.stack([np.stack([ma,mab,mca],1),np.stack([mab,mb,mbc],1),np.stack([mca,mbc,mc],1),np.stack([mab,mbc,mca],1)],axis=1)
    return xx.reshape(-1,3,3),mm.reshape(-1,3)


def constant_subtriangles(tris,mu,eps):
    starts=tris[:,[0,2,1],:].reshape(-1,3);ends=tris[:,[2,1,0],:].reshape(-1,3)
    strength=np.repeat(mu.mean(axis=1),3)
    return np.einsum('pej,e->pj',core.kernel(POINTS,starts,ends,eps),strength)


def fixture_triangles(vertices,nodes,coeff):
    cells=[];nrows,nspan=coeff.shape
    for r in range(nrows-1):
        for j in range(nspan):
            ids=[r*(nspan+1)+j,r*(nspan+1)+j+1,(r+1)*(nspan+1)+j+1,(r+1)*(nspan+1)+j]
            mu=np.array([coeff[r,j],coeff[r,j],coeff[r+1,j],coeff[r+1,j]])
            for t in TRIS:cells.append((np.array(ids)[t],mu[t]))
    return cells


def edge_assembly(cells,epsilons):
    edges={}
    for (ids,mu),eps in zip(cells,epsilons):
        for i,j in [(0,2),(2,1),(1,0)]:
            a,b=int(ids[i]),int(ids[j]);key=tuple(sorted((a,b)));sign=1 if (a,b)==key else -1
            value=np.array([mu[i],mu[j]]) if sign==1 else np.array([mu[j],mu[i]])
            if key not in edges:edges[key]=dict(mu=np.zeros(2),eps=eps,count=0)
            elif edges[key]['eps']!=eps:raise ValueError('incompatible_edge_kernel')
            edges[key]['mu']+=sign*value;edges[key]['count']+=1
    return edges


def grid_velocity(vertices,coeff,points,eps,reference=False):
    nrows,nspan=coeff.shape;result=np.zeros((len(points),3))
    for r in range(nrows-1):
        for j in range(nspan):
            ids=[r*(nspan+1)+j,r*(nspan+1)+j+1,(r+1)*(nspan+1)+j+1,(r+1)*(nspan+1)+j]
            values=[coeff[r,j],coeff[r+1,j]]
            if reference:result+=area_reference(vertices[ids],values,points,eps,80)
            else:result+=primary(vertices[ids],values,points,eps,line_order=96,surface_order=48)
    return result


def topology():
    vertices=np.array([[x,y,.04*r+.03*r*eta] for r,x in enumerate([1.,1.3,1.7]) for y,eta in zip([-.6,0.,.6],[-1.,0.,1.])])
    same=np.repeat(np.array([.8,.3,.1])[:,None],2,axis=1);coeff=np.array([[.8,-.2],[.3,.4],[.1,.05]])
    results=[]
    for name,a in [('matching',same),('different_strips',coeff)]:
        cells=fixture_triangles(vertices,None,a);edges=edge_assembly(cells,[.08]*len(cells))
        separate=np.zeros_like(POINTS);merged=np.zeros_like(POINTS)
        for ids,mu in cells:
            separate+=triangle_primary(POINTS,vertices[ids],mu,.08)
            merged+=triangle_primary(POINTS,vertices[ids],mu,.08,boundary=False)
        for (i,j),value in edges.items():merged+=line(POINTS,vertices[i],vertices[j],*value['mu'],.08)
        internal=max(float(np.max(np.abs(v['mu']))) for v in edges.values() if v['count']==2)
        results.append(dict(name=name,merged_difference=error(merged,separate),maximum_internal_coefficient=internal))
    averaged=np.repeat(coeff.mean(axis=1)[:,None],2,axis=1)
    original=grid_velocity(vertices,coeff,POINTS,.08);wrong=grid_velocity(vertices,averaged,POINTS,.08)
    refused=False
    cells=fixture_triangles(vertices,None,coeff)
    try:edge_assembly(cells,[.04]+[.08]*(len(cells)-1))
    except ValueError:refused=True
    old=coeff.copy();digest=hashlib.sha256(old.tobytes()).hexdigest();te=np.array([.45,.35])
    def trial(a):
        if a.shape!=(3,2):raise ValueError('wrong_mu_shape')
        return np.vstack([te,a[:-1]])
    first=trial(old);second=trial(old);wrong_shape=False
    try:trial(np.zeros((3,3)))
    except ValueError:wrong_shape=True
    return dict(vertices=vertices.tolist(),coefficients=coeff.tolist(),cases=results,
        shared_core_mismatch_refused=refused,wrong_vertex_mu_shape_refused=wrong_shape,
        averaged_strips_difference=error(wrong,original),trial_repeat_equal=bool(np.array_equal(first,second)),
        accepted_unchanged=hashlib.sha256(old.tobytes()).hexdigest()==digest,
        trial_mu=first.tolist(),transport_correct=bool(np.array_equal(first[0],te) and np.array_equal(first[1:],old[:-1])))


def manufactured():
    grid=geometry.make_grid(2,2,2.,wake=False);points=grid['collocation'];normals=grid['normals']
    abb=geometry.aic(grid,0.)
    by_edge=np.stack([line(points,grid['points'][i],grid['points'][j],1.,1.,0.,96) for i,j in grid['edges']],axis=1)
    ref=np.einsum('pej,pj->pe',by_edge,normals)@grid['source']
    vertices=np.array([[x,y,.03*r+.02*r*eta] for r,x in enumerate([1.,1.2,1.45]) for y,eta in zip([-1.,0.,1.],[-1.,0.,1.])])
    w=np.empty((4,6));wr=np.empty_like(w)
    for dof in range(6):
        a=np.zeros((3,2));a.flat[dof]=1.
        w[:,dof]=np.einsum('pj,pj->p',grid_velocity(vertices,a,points,0.),normals)
        wr[:,dof]=np.einsum('pj,pj->p',grid_velocity(vertices,a,points,0.,True),normals)
    t=np.zeros((2,4));t[0,2]=1;t[1,3]=1
    target=np.array([.3,.5,.4,.2]);old=np.array([.2,-.1,.15,-.05])
    rhs=ref@target+wr[:,:2]@(t@target)+wr[:,2:]@old
    matrix=abb+w[:,:2]@t;b=rhs-w[:,2:]@old
    solved=np.linalg.solve(matrix,b)
    return dict(bound_points=grid['points'].tolist(),wake_vertices=vertices.tolist(),target=target.tolist(),solved=solved.tolist(),
        bound_influence_difference=error(abb,ref),wake_influence_difference=error(w,wr),
        solution_relative=error(solved,target),residual_relative=float(np.linalg.norm(matrix@solved-b)/np.linalg.norm(b)),
        front_operator_shape=list(w[:,:2].shape),mu_shape=[3,2],geometry_shape=[3,3,3],
        limitation='manufactured normal-velocity RHS, not a physical finite-wing load benchmark')


def invalid_cases():
    quad=SHAPES['flat'];rows=[]
    for eps in [0.,.08]:
        for label,p in [('interior',[1.2,0.,0.]),('edge',[1.,0.,0.]),('vertex',quad[0])]:
            refused=False
            try:primary(quad,[1.,0.],np.array([p]),eps)
            except ValueError:refused=True
            rows.append(dict(name=label,eps=eps,refused=refused))
    for label,patch,points,eps,policy in [
        ('degenerate',np.zeros((4,3)),POINTS,0.,'off_surface_only'),
        ('negative_core',quad,POINTS,-.1,'off_surface_only'),
        ('nonfinite',quad,np.full((1,3),np.nan),0.,'off_surface_only'),
        ('self_policy',quad,POINTS,0.,'self')]:
        refused=False
        try:primary(patch,[1.,0.],points,eps,policy)
        except ValueError:refused=True
        rows.append(dict(name=label,refused=refused))
    return rows


def main():
    subdivision=[];constant=[];negative=[]
    for name,quad in SHAPES.items():
        for eps in [0.,.08]:
            basis=[]
            for values in [[1.,0.],[0.,1.],[1.,1.]]:
                budget();value=primary(quad,values,eps=eps);reference=area_reference(quad,values,eps=eps)
                RECORDS.append(dict(shape=name,eps=eps,coefficients=values,primary=value.tolist(),reference=reference.tolist(),
                    relative_difference=error(value,reference),eligible=point_eligible(value,reference)))
                basis.append(value)
                if values!=[1.,1.]:
                    edge=primary(quad,values,eps=eps,sheet=False)
                    soft=area_reference(quad,values,eps=eps,soft_gradient=True)
                    negative.append(dict(shape=name,eps=eps,coefficients=values,
                        missing_sheet_difference=error(edge,value),soft_gradient_difference=error(soft,value)))
                    mu=np.array([values[0],values[0],values[1],values[1]])
                    tris=np.array([quad[t] for t in TRIS]);muv=np.array([mu[t] for t in TRIS]);errors=[]
                    for level in range(6):
                        approx=constant_subtriangles(tris,muv,eps);errors.append(error(approx,value))
                        if level<5:tris,muv=subdivide(tris,muv)
                    orders=[float(np.log2(errors[n-1]/errors[n])) for n in [4,5]]
                    subdivision.append(dict(shape=name,eps=eps,coefficients=values,errors=errors,orders=orders))
            ring=core.kernel(POINTS,quad,np.roll(quad,-1,axis=0),eps).sum(axis=1)
            constant.append(dict(shape=name,eps=eps,constant_ring_difference=error(basis[2],ring),
                L_plus_R_difference=error(basis[0]+basis[1],ring),eligible=point_eligible(basis[2],ring) and point_eligible(basis[0]+basis[1],ring)))
            print(json.dumps(dict(stage='patch',shape=name,eps=eps,maximum_reference_difference=max(r['relative_difference'] for r in RECORDS if r['shape']==name and r['eps']==eps))),flush=True)
    a=.4;b=-.3;rz=np.array([[np.cos(a),-np.sin(a),0.],[np.sin(a),np.cos(a),0.],[0.,0.,1.]])
    ry=np.array([[np.cos(b),0.,np.sin(b)],[0.,1.,0.],[-np.sin(b),0.,np.cos(b)]])
    r=rz@ry;shift=np.array([.3,-.2,.5]);quad=SHAPES['warped'];values=[.7,-.2]
    base=primary(quad,values,eps=.08)
    rotated=primary(quad@r.T+shift,values,POINTS@r.T+shift,.08)
    scaled=primary(2*quad,values,2*POINTS,.16)
    complex_values=np.array([.7+.3j,-.2+.15j]);complex_v=primary(quad,complex_values,eps=.08)
    separate=primary(quad,complex_values.real,eps=.08)+1j*primary(quad,complex_values.imag,eps=.08)
    transformations=dict(rotation=error(rotated,base@r.T),scale=error(scaled,base/2),complex_linearity=error(complex_v,separate))
    collapse=[];p=np.array([[.2,0.,0.],[.75,0.,0.]])
    for values in [[1.,0.],[0.,1.]]:
        velocity=primary(SHAPES['flat'],values,p,0.)
        left,right=planar.weights(1.,.4,p[:,0],.6);ref=values[0]*left+values[1]*right
        collapse.append(dict(coefficients=values,z_difference=error(velocity[:,2],ref),maximum_inplane=float(np.max(np.abs(velocity[:,:2])))))
    topo=topology();manufact=manufactured();invalid=invalid_cases()
    checks=dict(independent_surface_integral=all(r['eligible'] for r in RECORDS),
        constant_ring_and_basis=all(r['eligible'] for r in constant),
        subdivision_accuracy=all(r['errors'][-1]<=.001 for r in subdivision),
        subdivision_orders=all(min(r['orders'])>=1.7 for r in subdivision),
        geometric_covariance=max(transformations.values())<=1e-11,
        planar_reduction=all(r['z_difference']<=1e-9 and r['maximum_inplane']<=1e-10 for r in collapse),
        negative_missing_sheet=max(r['missing_sheet_difference'] for r in negative)>.01,
        negative_soft_gradient=max(r['soft_gradient_difference'] for r in negative if r['eps']>0)>.01,
        edge_aggregation=all(r['merged_difference']<=1e-11 for r in topo['cases']),
        matching_trace_cancellation=topo['cases'][0]['maximum_internal_coefficient']<=1e-12,
        negative_strip_average=topo['averaged_strips_difference']>.001,
        shared_core_refusal=topo['shared_core_mismatch_refused'],mu_shape_refusal=topo['wrong_vertex_mu_shape_refused'],
        trial_transport=topo['trial_repeat_equal'] and topo['accepted_unchanged'] and topo['transport_correct'],
        manufactured_solution=manufact['solution_relative']<=1e-8,
        manufactured_residual=manufact['residual_relative']<=1e-11,
        invalid_inputs=all(r['refused'] for r in invalid))
    source=Path(__file__).resolve()
    loaded=sorted({Path(module.__file__).resolve() for module in list(sys.modules.values()) if getattr(module,'__file__',None)
        and Path(module.__file__).resolve().parent==source.parent and Path(module.__file__).resolve()!=source})
    result=dict(protocol='docs/verification/uvlm-wp1-three-dimensional-experiment.md',representation_id='W-P1-stream/P0-span-facet-v1',
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in loaded},
        shapes={name:x.tolist() for name,x in SHAPES.items()},observation_points=POINTS.tolist(),
        patch_cases=RECORDS,constant_cases=constant,subdivision_cases=subdivision,negative_controls=negative,
        transformations=transformations,planar_reduction=collapse,topology=topo,manufactured=manufact,invalid_cases=invalid,
        checks=checks,all_criteria_pass=bool(all(checks.values())),numpy=np.__version__,scipy=scipy.__version__,
        seconds=time.monotonic()-START,peak_process_working_bytes=core.memory.memory_bytes(),
        limits='off-surface induction and manufactured coupling only; no physical finite-wing loads, self/free wake or ANCF coupling')
    (ROOT/'docs/verification/uvlm-wp1-three-dimensional-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=core.json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,manufactured=manufact,seconds=result['seconds']),default=core.json_scalar),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':
    try:main()
    except Exception as error:
        (ROOT/'build/analysis/uvlm-wp1-three-dimensional-stopped.json').write_text(json.dumps(dict(error=repr(error),partial_patch_cases=RECORDS,seconds=time.monotonic()-START),indent=2,allow_nan=False,default=core.json_scalar)+'\n',encoding='utf-8')
        raise
