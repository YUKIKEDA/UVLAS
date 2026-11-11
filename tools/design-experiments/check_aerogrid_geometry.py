"""Unified candidate geometry; docs/verification/aerogrid-camber-addedmass-experiment.md."""
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad
import check_vortex_core_load as core
import check_q9_interface_reparameterization as interface
import check_q9_curved_load as curved
import check_q9_independent_load as common
import check_q9_plate_convergence as memory

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()
GEOMETRY='dual-quarter-v1'


def make_grid(m,n,span,camber=0.,wake=True,mapper=None,normal_mapper=None):
    if mapper is not None and normal_mapper is None:raise ValueError('mapped_surface_requires_its_normal')
    xv=np.r_[(np.arange(m)+.25)/m,1.];etas=np.linspace(-1.,1.,n+1)
    material=np.array([[x,y] for x in xv for y in etas])
    def analytic(x,e):return np.array([x,span*e/2,4*camber*x*(1-x)])
    points=np.array([(mapper(x,e) if mapper else analytic(x,e)) for x,e in material])
    if wake:points=np.concatenate([points,points[-(n+1):]+np.array([20*span,0.,0.])])
    faces=[];owners=[]
    for i in range(m+int(wake)):
        for j in range(n):
            a=i*(n+1)+j;b=a+1;c=a+n+2;d=a+n+1
            faces.append([a,b,c,d]);owners.append(min(i,m-1)*n+j)
    edges=sorted({tuple(sorted((f[k],f[(k+1)%4]))) for f in faces for k in range(4)})
    ei={edge:i for i,edge in enumerate(edges)};source=np.zeros((len(edges),m*n));body=source.copy()
    for row,f in enumerate(faces):
        for k in range(4):
            pair=(f[k],f[(k+1)%4]);edge=tuple(sorted(pair));sign=1 if pair==edge else -1
            source[ei[edge],owners[row]]+=sign
            if row<m*n:body[ei[edge],owners[row]]+=sign
    te=set(range(m*(n+1),(m+1)*(n+1)))
    body_ids=np.array([i for i,e in enumerate(edges) if np.any(body[i]) and not set(e).issubset(te)])
    triangles=[];tri_owner=[]
    for p,(a,b,c,d) in enumerate(faces[:m*n]):
        triangles.extend([[a,c,b],[a,d,c]]);tri_owner.extend([p,p])
    coll_material=np.array([[(i+.75)/m,(etas[j]+etas[j+1])/2] for i in range(m) for j in range(n)])
    collocation=np.array([(mapper(x,e) if mapper else analytic(x,e)) for x,e in coll_material])
    normals=np.array([(normal_mapper(x,e) if normal_mapper else np.array([-4*camber*(1-2*x),0.,1.])) for x,e in coll_material])
    normals/=np.linalg.norm(normals,axis=1)[:,None]
    return dict(spec=GEOMETRY,m=m,n=n,material=material,coll_material=coll_material,points=points,
        faces=np.array(faces[:m*n]),edges=edges,source=source,body=body,body_ids=body_ids,
        triangles=np.array(triangles),tri_owner=np.array(tri_owner),collocation=collocation,normals=normals)


def aic(grid,eps):
    a=np.array([grid['points'][i] for i,j in grid['edges']]);b=np.array([grid['points'][j] for i,j in grid['edges']])
    influences=core.kernel(grid['collocation'],a,b,eps)
    return np.einsum('pej,pj->pe',influences,grid['normals'])@grid['source']


def circulating_forces(grid,gamma,flow,eps):
    a=np.array([grid['points'][i] for i,j in grid['edges']]);b=np.array([grid['points'][j] for i,j in grid['edges']])
    ids=grid['body_ids'];centers=.5*(a[ids]+b[ids]);dl=b[ids]-a[ids]
    excluded=ids[:,None]==np.arange(len(grid['edges']))[None,:]
    induced=np.einsum('pej,e->pj',core.kernel(centers,a,b,eps,excluded),grid['source']@gamma)
    strengths=grid['body']@gamma
    forces=1.2*strengths[ids,None]*np.cross(flow+induced,dl)
    return centers,forces


def noncirculating(grid,dot):
    p=grid['points'][grid['triangles']]
    area=.5*np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);centers=p.mean(axis=1)
    return centers,1.2*dot[grid['tri_owner'],None]*area,area


def q9_operator(x,eta):
    shape=interface.scalar_shape(2*x-1,eta,.003)[0]
    return np.kron(shape.reshape(1,27),np.eye(3))


def mapping_checks():
    q=np.zeros((9,3,3))
    for i,(xi,eta) in enumerate(common.nodes):
        x=(xi+1)/2
        q[i,0]=[x,eta/2,.06*x*(1-x)+.04*eta*eta+.02*x*eta]
        q[i,1]=[.04*xi,.02*eta,1+.03*xi*eta]
        q[i,2]=[.02,.01*xi,-.01*eta]
    q=q.reshape(81)
    mapper=lambda x,e:q9_operator(x,e)@q
    def normal_mapper(x,e):
        _,gx,gy,_=interface.scalar_shape(2*x-1,e,.003)
        tx=2*np.kron(gx.reshape(1,27),np.eye(3))@q
        ty=np.kron(gy.reshape(1,27),np.eye(3))@q
        return np.cross(tx,ty)
    grid=make_grid(3,2,1.,wake=False,mapper=mapper,normal_mapper=normal_mapper)
    normal_error=0.
    for (x,e),normal in zip(grid['coll_material'],grid['normals']):
        delta=1e-6
        tx=(mapper(x+delta,e)-mapper(x-delta,e))/(2*delta)
        ty=(mapper(x,e+delta)-mapper(x,e-delta))/(2*delta)
        value=np.cross(tx,ty);value/=np.linalg.norm(value)
        normal_error=max(normal_error,float(np.linalg.norm(value-normal)))
    try:make_grid(3,2,1.,wake=False,mapper=mapper)
    except ValueError:missing_normal_refused=True
    else:missing_normal_refused=False
    h=np.array([q9_operator(x,e) for x,e in grid['material']]);tri=grid['triangles']
    tri_h=h[tri].mean(axis=1);ids=grid['body_ids'];edge_nodes=np.array([grid['edges'][i] for i in ids])
    edge_h=h[edge_nodes].mean(axis=1)
    rng=np.random.default_rng(20261030);edge_force=rng.standard_normal((len(ids),3))
    dot=.2+.03*np.arange(6)
    def generalized(state):
        points=np.einsum('pai,i->pa',h,state)
        p=points[tri];area=.5*np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])
        force=1.2*dot[grid['tri_owner'],None]*area
        return np.einsum('pai,pa->i',tri_h,force)+np.einsum('pai,pa->i',edge_h,edge_force),force
    direct,force=generalized(q)
    nodal=np.zeros_like(grid['points'])
    for nodes,f in zip(edge_nodes,edge_force):
        for node in nodes:nodal[node]+=.5*f
    for nodes,f in zip(tri,force):
        for node in nodes:nodal[node]+=f/3
    scattered=np.einsum('pai,pa->i',h,nodal)
    scale=interface.natural_scale(.02,9)
    qerror=float(np.linalg.norm(scale*(direct-scattered))/np.linalg.norm(scale*direct))
    velocity=1e-3*scale*rng.standard_normal(81)
    power_direct=float(velocity@direct)
    tri_velocity=np.einsum('pai,i->pa',tri_h,velocity);edge_velocity=np.einsum('pai,i->pa',edge_h,velocity)
    power_load=float(np.sum(tri_velocity*force)+np.sum(edge_velocity*edge_force))
    center=np.array([.2,-.1,.3]);rigid=curved.rigid_modes(q,center)
    tri_points=grid['points'][tri].mean(axis=1);edge_points=grid['points'][edge_nodes].mean(axis=1)
    resultant=force.sum(axis=0)+edge_force.sum(axis=0)
    moment=np.cross(tri_points-center,force).sum(axis=0)+np.cross(edge_points-center,edge_force).sum(axis=0)
    wrench=np.r_[resultant,moment];wrench_error=float(np.linalg.norm(rigid.T@direct-wrench)/max(1.,np.linalg.norm(wrench)))
    direction=1e-4*scale*rng.standard_normal(81)
    dp=np.einsum('pai,i->pa',h,direction)[tri];p=grid['points'][tri]
    da=.5*(np.cross(dp[:,1]-dp[:,0],p[:,2]-p[:,0])+np.cross(p[:,1]-p[:,0],dp[:,2]-dp[:,0]))
    analytic=np.einsum('pai,pa->i',tri_h,1.2*dot[grid['tri_owner'],None]*da)
    cs=generalized(q.astype(complex)+1j*1e-25*direction)[0].imag/1e-25
    derivative=float(np.linalg.norm(scale*(cs-analytic))/np.linalg.norm(scale*analytic))
    averaged_material=grid['material'][tri].mean(axis=1)
    wrong_h=np.array([q9_operator(x,e) for x,e in averaged_material])
    wrong_positions=np.einsum('pai,i->pa',wrong_h,q)
    mismatch=float(np.max(np.linalg.norm(wrong_positions-tri_points,axis=1)))
    wrong_power=float(abs(np.sum(np.einsum('pai,i->pa',wrong_h-tri_h,velocity)*force)))
    angle=.5;r=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
    shift=np.array([.3,-.4,.5]);rotated=curved.transform(q,r,shift)
    geometric=float(np.max(np.abs(np.einsum('pai,i->pa',h,rotated)-(grid['points']@r.T+shift))))
    return dict(geometry_id=GEOMETRY,material_vortex_coordinates=grid['material'].tolist(),
        material_collocation_coordinates=grid['coll_material'].tolist(),triangle_indices=tri.tolist(),
        dual_load_error=qerror,power_error=float(abs(power_direct-power_load)/max(1.,abs(power_load))),
        wrench_error=wrench_error,area_tangent_relative=derivative,geometry_covariance_error=geometric,
        normal_finite_difference_error=normal_error,missing_normal_refused=missing_normal_refused,
        wrong_smooth_centroid_position_error=mismatch,wrong_smooth_centroid_power_difference=wrong_power,
        eligible=bool(max(qerror,abs(power_direct-power_load)/max(1.,abs(power_load)),wrench_error,geometric)<=1e-11 and derivative<=1e-10 and normal_error<=1e-8 and missing_normal_refused),
        negative_detected=bool(mismatch>1e-8 or wrong_power>1e-8))


def camber_case(m,n):
    span=10000.;speed=8.;grid=make_grid(n,1,span,m,wake=True);eps=.001/n
    matrix=aic(grid,eps);flow=np.array([speed,0.,0.]);rhs=-grid['normals']@flow
    circulation=np.linalg.solve(matrix,rhs);points,forces=circulating_forces(grid,circulation,flow,eps)
    total=forces.sum(axis=0);moment=np.cross(points-np.array([.25,0,0]),forces).sum(axis=0)
    q=.5*1.2*speed**2
    cl=float(total[2]/(q*span));cm=float(moment[1]/(q*span));cd=float(total[0]/(q*span))
    return dict(geometry_id=GEOMETRY,camber=m,n=n,CL=cl,Cm_quarter=cm,CD=cd,
        CL_reference=4*np.pi*m,Cm_reference=-np.pi*m,CL_relative=float(abs(cl/(4*np.pi*m)-1)),
        Cm_relative=float(abs(cm/(-np.pi*m)-1)),nonpenetration_relative=float(np.linalg.norm(matrix@circulation-rhs)/np.linalg.norm(rhs)))


def added_mass_case(n):
    span=10000.;grid=make_grid(n,1,span,wake=False);matrix=aic(grid,0.)
    rhs=np.ones(n);dot=np.linalg.solve(matrix,rhs)
    points,forces,area=noncirculating(grid,dot)
    total=forces.sum(axis=0)/span;moment=np.cross(points,forces).sum(axis=0)/span
    reference=1.2*np.pi/4
    wrong_force=float(1.2*np.sum(dot/n));wrong_moment=float(-1.2*np.sum(dot*(np.arange(n)+.5)/n**2))
    return dict(geometry_id=GEOMETRY,n=n,Fz_per_span=float(total[2]),My_le_per_span=float(moment[1]),
        reference_added_mass_per_span=reference,force_relative=float(abs(total[2]/(-reference)-1)),
        moment_relative=float(abs(moment[1]/(reference/2)-1)),
        nonpenetration_relative=float(np.linalg.norm(matrix@dot-rhs)/np.linalg.norm(rhs)),
        wrong_physical_panel_force_relative=abs(wrong_force/(-reference)-1),
        wrong_physical_panel_moment_relative=abs(wrong_moment/(reference/2)-1),
        vortex_coordinates=grid['material'].tolist(),collocation_coordinates=grid['coll_material'].tolist())


def main():
    mapping=mapping_checks();camber=[];added=[]
    for m in [.02,.01,.005]:
        for n in [16,32,64,128]:
            if time.monotonic()-START>120 or memory.memory_bytes()>1024**3:raise TimeoutError('experiment budget')
            case=camber_case(m,n);camber.append(case);print(json.dumps(dict(stage='camber',**case)),flush=True)
    for n in [1,2,4,8,16,32,64]:
        case=added_mass_case(n);added.append(case)
        print(json.dumps(dict(stage='added_mass',**{k:v for k,v in case.items() if 'coordinates' not in k})),flush=True)
    integralf=quad(lambda theta:-.6*np.sin(theta)**2,0,np.pi,epsabs=1e-12,epsrel=1e-12)[0]
    integralm=quad(lambda theta:.3*(1-np.cos(theta))*np.sin(theta)**2,0,np.pi,epsabs=1e-12,epsrel=1e-12)[0]
    fine=camber[-1];afine=added[-1]
    checks=dict(mapping=mapping['eligible'],negative_centroid=mapping['negative_detected'],
        camber_lift=fine['CL_relative']<=.02,camber_moment=fine['Cm_relative']<=.02,
        camber_residual=all(c['nonpenetration_relative']<=1e-12 for c in camber),
        added_mass_force=afine['force_relative']<=1e-3,added_mass_moment=afine['moment_relative']<=1e-3,
        added_mass_residual=all(c['nonpenetration_relative']<=1e-12 for c in added),
        analytic_integrals=abs(integralf+1.2*np.pi/4)<=1e-12 and abs(integralm-1.2*np.pi/8)<=1e-12,
        negative_wrong_area=any(max(c['wrong_physical_panel_force_relative'],c['wrong_physical_panel_moment_relative'])>.1 for c in added))
    source=Path(__file__)
    result=dict(protocol='docs/verification/aerogrid-camber-addedmass-experiment.md',geometry_id=GEOMETRY,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),numpy=np.__version__,scipy=scipy.__version__,
        dependency_sha256={name:hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in ['check_vortex_core_load.py','check_q9_interface_reparameterization.py','check_q9_curved_load.py','check_q9_independent_load.py','check_q9_plate_convergence.py']},
        mapping=mapping,camber_cases=camber,added_mass_cases=added,checks=checks,all_criteria_pass=bool(all(checks.values())),
        seconds=time.monotonic()-START,peak_process_working_bytes=memory.memory_bytes(),
        limits='one geometry candidate, small-camber theory and instantaneous noncirculatory component; no wake-time integration or ANCF coupled dynamics')
    (ROOT/'docs/verification/aerogrid-camber-addedmass-result.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=core.json_scalar)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],mapping=mapping),default=core.json_scalar),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':main()
