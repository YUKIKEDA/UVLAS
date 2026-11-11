"""Multi-element manufactured equilibrium and prestress spectrum experiment.
Protocol: docs/verification/q9-multielement-prestress-convergence.md.
Not a product nonlinear/multielement solver.
"""
import hashlib
import itertools
import json
import math
import time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import scipy
import skfem
from scipy.sparse import coo_matrix,diags,bmat,block_diag
from scipy.sparse.linalg import splu,eigsh
from scipy.optimize import linear_sum_assignment
from numpy.polynomial.legendre import leggauss
from skfem import MeshHex,ElementVector,ElementHex2,Basis,BilinearForm,asm
import check_q9_independent_load as common
import check_q9_curved_load as curved
import check_q9_plate_convergence as plate
import check_q9_interface_reparameterization as interface
import check_q9_material_tangent as material
import check_q9_prestress_equilibrium as prestress

L=1.;H=.1;A=.0005;B=.0005;YOUNG=70e9;NU=.3;RHO=2700.
MU=YOUNG/(2*(1+NU));LAME=YOUNG*NU/((1+NU)*(1-2*NU))
E0=np.diag([-.001,.0003,.0003]);C=material.isotropic(YOUNG,NU)
W9=interface.natural_scale(H,1)


def exact_fields(x):
    t=x/L;f=1+3*A/L*t*t;g=3*B/L*t*t
    fp=6*A/L**2*t;gp=6*B/L**2*t
    ff=np.zeros(np.shape(x)+(3,3));ff[...,0,0]=f;ff[...,1,1]=ff[...,2,2]=1;ff[...,2,0]=g
    ee=.5*(np.einsum('...ji,...jk->...ik',ff,ff)-np.eye(3))
    elastic=ee-E0;trace=np.trace(elastic,axis1=-2,axis2=-1)
    ss=LAME*trace[...,None,None]*np.eye(3)+2*MU*elastic
    pp=np.einsum('...ij,...jk->...ik',ff,ss)
    s11p=(LAME+2*MU)*(f*fp+g*gp)
    body=np.zeros(np.shape(x)+(3,));body[...,0]=-(fp*ss[...,0,0]+f*s11p)
    body[...,2]=-(gp*(ss[...,0,0]+MU)+g*s11p)
    energy=.5*LAME*trace**2+MU*np.sum(elastic*elastic,axis=(-2,-1))
    displacement=np.stack([A*t**3,np.zeros_like(t),B*t**3],axis=-1)
    return dict(f=ff,e=ee,s=ss,p=pp,body=body,energy=energy,u=displacement)


def context(n,order):
    pts,weights=leggauss(order);ijk=list(itertools.product(range(order),repeat=3))
    coords=np.array([[pts[i],pts[j],pts[k]] for i,j,k in ijk])
    w=np.array([weights[i]*weights[j]*weights[k] for i,j,k in ijk])*H/(8*n*n)
    s,g=interface.shape_tables(coords,0.,H/2)
    g=g*np.array([2*n,2*n,2/H])[None,None,:]
    ctx=dict(w=w,g=g,axes=np.broadcast_to(np.eye(3),(len(s),3,3)).copy(),
        e0=np.broadcast_to(E0,(len(s),3,3)).copy(),c=C,c4=material.fourth(C))
    return ctx,s


def mesh_data(n):
    side=2*n+1;q0=np.zeros((side,side,3,3))
    for j,y in enumerate(np.linspace(-.5,.5,side)):
        for i,x in enumerate(np.linspace(0,1,side)):
            q0[j,i,0]=[x,y,0];q0[j,i,1]=[0,0,1]
    offsets=(common.nodes+1).astype(int);cells=[]
    for j in range(n):
        for i in range(n):
            nodes=(2*j+offsets[:,1])*side+2*i+offsets[:,0]
            cells.append(dict(i=i,j=j,dofs=(9*nodes[:,None]+np.arange(9)).reshape(81)))
    fixed=np.concatenate([np.arange(9*j*side,9*j*side+9) for j in range(side)])
    free=np.setdiff1d(np.arange(9*side*side),fixed)
    return q0.reshape(-1),cells,free,fixed


def face_load(q0,axis,sign,order=5):
    pts,weights=leggauss(order);points=[];w=[]
    for i,j in itertools.product(range(order),repeat=2):
        point=np.zeros(3);point[axis]=sign
        point[[a for a in range(3) if a!=axis]]=[pts[i],pts[j]]
        points.append(point);w.append(weights[i]*weights[j])
    s,g=interface.shape_tables(points,0.,H/2)
    jac=np.einsum('ia,pib->pab',q0.reshape(27,3),g)
    area=sign*np.cross(jac[:,:,(axis+1)%3],jac[:,:,(axis+2)%3])
    x=s@q0.reshape(27,3);p=exact_fields(x[:,0])['p']
    traction=np.einsum('pab,pb->pa',p,area)*np.array(w)[:,None]
    return (s.T@traction).reshape(81)


def assemble_matrices(values,cells,size):
    rows=np.concatenate([np.repeat(c['dofs'],81) for c in cells])
    cols=np.concatenate([np.tile(c['dofs'],81) for c in cells])
    data=np.concatenate([v.reshape(-1) for v in values])
    result=coo_matrix((data,(rows,cols)),shape=(size,size)).tocsr();result.eliminate_zeros()
    return result


def loads(q0,cells,n,ctx,s):
    local=[];global_force=np.zeros(len(q0))
    for cell in cells:
        qe=q0[cell['dofs']];x=s@qe.reshape(27,3)
        force=(s.T@(ctx['w'][:,None]*exact_fields(x[:,0])['body'])).reshape(81)
        boundaries=[]
        if cell['i']==n-1:boundaries.append((0,1))
        if cell['j']==0:boundaries.append((1,-1))
        if cell['j']==n-1:boundaries.append((1,1))
        boundaries.extend([(2,-1),(2,1)])
        for axis,sign in boundaries:force+=face_load(qe,axis,sign)
        local.append(force);np.add.at(global_force,cell['dofs'],force)
    return local,global_force


def evaluate(q,cells,ctx,need_k=True):
    force=np.zeros(len(q));matrices=[];energies=[]
    for cell in cells:
        state=material.analytic_evaluation(q[cell['dofs']],ctx)
        np.add.at(force,cell['dofs'],state['force']);energies.append(state['energy'])
        if need_k:matrices.append(state['k'])
    return force,(assemble_matrices(matrices,cells,len(q)) if need_k else None),math.fsum(energies)


def orientation_bound(q,n):
    polys=curved.geometry_polynomials(q,H);denominators=[F(1,2*n),F(1,2*n),F.from_float(H)/2]
    bounds=[]
    for a in range(3):
        for b in range(3):
            p=curved.times(curved.derivative(polys[a],b),1/denominators[b])
            if a==b:p=curved.add(p,{(0,0,0):F(-1)})
            interval=curved.bernstein_bounds(p)
            bounds.append(max(abs(interval['lower']),abs(interval['upper'])))
    squared=[math.nextafter(v*v,math.inf) for v in bounds]
    radius=math.nextafter(math.sqrt(math.nextafter(math.fsum(squared),math.inf)),math.inf)
    lower=None
    if radius<1:
        exact=F.from_float(H)/F(8*n*n)*(1-F.from_float(radius))**3
        lower=math.nextafter(float(exact),-math.inf)
    return dict(eligible=bool(radius<1),deformation_from_identity_upper=radius,parent_det_lower=lower)


def exact_det_certificate(q):
    return interface.geometry_certificate(q,-H/2,H/2)


def modes(k,m,kcount=6):
    d=1/np.sqrt(m.diagonal());scale=diags(d)
    a=(scale@k@scale).tocsc();b=(scale@m@scale).tocsc()
    values,vectors=eigsh(a,k=kcount,M=b,sigma=0.,which='LM',tol=1e-10,
        v0=np.random.default_rng(20261030).standard_normal(k.shape[0]))
    order=np.argsort(values);values=values[order];vectors=d[:,None]*vectors[:,order]
    residuals=[]
    for value,v in zip(values,vectors.T):
        denominator=np.linalg.norm(k@v)+abs(value)*np.linalg.norm(m@v)
        residuals.append(float(np.linalg.norm(k@v-value*(m@v))/denominator))
    orth=float(np.max(np.abs(vectors.T@m@vectors-np.eye(kcount))))
    return values,vectors,dict(eligible=bool(np.all(values>0) and max(residuals)<=1e-8 and orth<=1e-8),
        eigenpair_residual_max=max(residuals),mass_orthogonality=orth)


def static_errors(q,q0,cells,n):
    ctx,s=context(n,7);uerr=0.;unorm=0.;eerr=0.;enorm=0.;energy=0.;eref=0.
    for cell in cells:
        qe=q[cell['dofs']];ref=q0[cell['dofs']];xyz=s@ref.reshape(27,3)
        exact=exact_fields(xyz[:,0]);u=s@(qe-ref).reshape(27,3)
        f=np.einsum('ia,pib->pab',qe.reshape(27,3),ctx['g'])
        e=.5*(np.einsum('pji,pjk->pik',f,f)-np.eye(3))
        uerr+=float(np.sum(ctx['w']*np.sum((u-exact['u'])**2,axis=1)))
        unorm+=float(np.sum(ctx['w']*np.sum(exact['u']**2,axis=1)))
        eerr+=float(np.sum(ctx['w']*np.sum((e-exact['e'])**2,axis=(1,2))))
        enorm+=float(np.sum(ctx['w']*np.sum(exact['e']**2,axis=(1,2))))
        energy+=material.analytic_evaluation(qe,ctx)['energy'];eref+=float(ctx['w']@exact['energy'])
    return dict(displacement_l2=float(np.sqrt(uerr/unorm)),green_l2=float(np.sqrt(eerr/enorm)),
                energy_relative=float(abs(energy/eref-1)),reference_energy=eref)


def wrench(q,force):
    a=q.reshape(-1,3,3);f=force.reshape(-1,3,3)
    return np.concatenate([f[:,0,:].sum(axis=0),np.cross(a,f).sum(axis=(0,1))])


def run_internal(n):
    start=time.monotonic();q0,cells,free,fixed=mesh_data(n);ctx,s=context(n,5)
    local_load,load=loads(q0,cells,n,ctx,s);w=np.tile(W9,len(q0)//9);d=w[free]
    local_mass=np.kron(RHO*(s.T@(ctx['w'][:,None]*s)),np.eye(3))
    mass=assemble_matrices([local_mass]*len(cells),cells,len(q0));q=q0.copy();history=[];converged=False
    for iteration in range(21):
        if time.monotonic()-start>240:raise TimeoutError('internal case budget')
        fint,k,energy=evaluate(q,cells,ctx);residual=fint-load
        rhs=d*residual[free];kfree=k[free,:][:,free];scaled=diags(d)@kfree@diags(d)
        delta=splu(scaled.tocsc()).solve(-rhs);dq=np.zeros(len(q));dq[free]=d*delta
        rnorm=float(np.linalg.norm(rhs)/np.linalg.norm(d*load[free]))
        correction=max(float(np.max(np.abs(s@dq[c['dofs']].reshape(27,3)))) for c in cells)
        row=dict(iteration=iteration,residual=rnorm,candidate_physical_step=correction);history.append(row)
        if rnorm<=1e-9 and correction<=1e-10:converged=True;break
        if iteration==20:break
        slope=float(rhs@delta)
        if slope>=0:break
        accepted=False
        for reduction in range(13):
            step=2.**(-reduction);candidate=q+step*dq;change=0.;admissible=True
            for cell,fe in zip(cells,local_load):
                qe=candidate[cell['dofs']]
                f=np.einsum('ia,pib->pab',qe.reshape(27,3),ctx['g'])
                if np.any(np.linalg.det(f)<=0):admissible=False;break
                change+=prestress.potential_increment(q[cell['dofs']],qe,ctx,fe)
            if admissible and change<=1e-4*step*slope:
                q=candidate;row.update(step=step,potential_change=change);accepted=True;break
        if not accepted:break
    report=dict(n=n,elements=len(cells),total_coordinates=len(q0),free_coordinates=len(free),newton=history,converged=converged)
    if not converged:return dict(**report,eligible=False,reason='Newton_not_converged'),None
    certificates=[orientation_bound(q[c['dofs']],n) for c in cells]
    geometry=dict(eligible=all(c['eligible'] for c in certificates),
        maximum_radius=max(c['deformation_from_identity_upper'] for c in certificates),
        minimum_parent_det_lower=min((c['parent_det_lower'] for c in certificates if c['parent_det_lower'] is not None),default=None))
    if n==2:
        exact=[exact_det_certificate(q[c['dofs']]) for c in cells]
        geometry['det_polynomial_crosscheck']=all(c['exact_positive'] for c in exact)
        geometry['eligible']=geometry['eligible'] and geometry['det_polynomial_crosscheck']
    fint,k,_=evaluate(q,cells,ctx);residual=fint-load;reaction=np.zeros_like(q);reaction[fixed]=residual[fixed]
    wr=wrench(q,reaction);balance=wrench(q,load+reaction)
    force_scale=7e6
    reactions=dict(force=(wr[:3]).tolist(),moment=wr[3:].tolist(),
        force_reference_error=float(np.linalg.norm(wr[:3]-[-force_scale,0,0])/force_scale),
        balance=float(np.linalg.norm(balance)/force_scale))
    values,vectors,modal=modes(k[free,:][:,free],mass[free,:][:,free])
    full_vectors=np.zeros((len(q),6));full_vectors[free]=vectors
    coefficients=full_vectors.reshape(2*n+1,2*n+1,3,3,6)
    errors=static_errors(q,q0,cells,n)
    right=q.reshape(2*n+1,2*n+1,3,3)[n,-1,0]-q0.reshape(2*n+1,2*n+1,3,3)[n,-1,0]
    report.update(geometry=geometry,static_errors=errors,reaction=reactions,modes=modal,
        frequencies_hz=(np.sqrt(values)/(2*np.pi)).tolist(),right_center_displacement=right.tolist(),
        seconds=time.monotonic()-start,peak_process_working_bytes=plate.memory_bytes())
    report['eligible']=bool(geometry['eligible'] and modal['eligible'] and reactions['force_reference_error']<=1e-6 and reactions['balance']<=1e-8 and report['seconds']<=240 and report['peak_process_working_bytes']<=4*1024**3)
    return report,dict(coefficients=coefficients,q=q,q0=q0,cells=cells,ctx=ctx,free=free,fixed=fixed,load=load,mass=mass,k=k,values=values)


def external_reference(n,verify_component_form=False):
    start=time.monotonic()
    mesh=MeshHex.init_tensor(np.linspace(0,1,n+1),np.linspace(-.5,.5,n+1),np.array([-H/2,H/2]))
    basis=Basis(mesh,ElementHex2(),intorder=9)
    cached=exact_fields(basis.global_coordinates()[0])
    f_coefficient=np.moveaxis(cached['f'],[-2,-1],[0,1])
    s_coefficient=np.moveaxis(cached['s'],[-2,-1],[0,1])
    del cached
    fft=np.einsum('ai...,bi...->ab...',f_coefficient,f_coefficient)
    blocks=[]
    for a in range(3):
        row=[]
        for b in range(3):
            coefficients=[]
            for i,j in itertools.product(range(3),repeat=2):
                value=LAME*f_coefficient[a,i]*f_coefficient[b,j]+MU*f_coefficient[a,j]*f_coefficient[b,i]
                if i==j:value=value+MU*fft[a,b]
                if a==b:value=value+s_coefficient[i,j]
                if np.any(value):coefficients.append((i,j,value))
            @BilinearForm
            def block_tangent(u,v,w):
                if time.monotonic()-start>240:raise TimeoutError('external assembly budget')
                return sum((value*v.grad[i]*u.grad[j] for i,j,value in coefficients))
            row.append(asm(block_tangent,basis).tocsr())
        blocks.append(row)
    @BilinearForm
    def scalar_mass(u,v,w):return RHO*u*v
    k=bmat(blocks,format='csr');ms=asm(scalar_mass,basis).tocsr();m=block_diag([ms,ms,ms],format='csr')
    crosscheck=None
    if verify_component_form:
        vector_basis=Basis(mesh,ElementVector(ElementHex2()),intorder=9)
        @BilinearForm
        def direct_tangent(u,v,w):
            a=np.einsum('ai...,aj...->ij...',f_coefficient,u.grad);b=np.einsum('ai...,aj...->ij...',f_coefficient,v.grad)
            a=.5*(a+np.swapaxes(a,0,1));b=.5*(b+np.swapaxes(b,0,1))
            return LAME*np.einsum('ii...->...',a)*np.einsum('ii...->...',b)+2*MU*np.sum(a*b,axis=(0,1))+np.einsum('ai...,ij...,aj...->...',u.grad,s_coefficient,v.grad)
        @BilinearForm
        def direct_mass(u,v,w):return RHO*np.einsum('i...,i...->...',u,v)
        kd=asm(direct_tangent,vector_basis).toarray();md=asm(direct_mass,vector_basis).toarray()
        permutation=np.zeros((vector_basis.N,3*basis.N))
        for component,ids in enumerate(vector_basis.split_indices()):
            for idx in ids:
                xyz=vector_basis.doflocs[:,idx]
                matches=np.flatnonzero(np.max(np.abs(basis.doflocs-xyz[:,None]),axis=0)<1e-12)
                if len(matches)!=1:raise ValueError('component direct-sum correspondence')
                permutation[idx,component*basis.N+int(matches[0])]=1.
        differences=dict(stiffness=common.matrix_error(permutation.T@kd@permutation,k.toarray()),mass=common.matrix_error(permutation.T@md@permutation,m.toarray()))
        crosscheck=dict(n=n,vector_dofs=int(vector_basis.N),differences=differences,
            eligible=all(v<=1e-12 for x in differences.values() for v in x.values()))
    left=np.flatnonzero(np.abs(basis.doflocs[0])<1e-12)
    fixed=np.concatenate([left+component*basis.N for component in range(3)]);free=np.setdiff1d(np.arange(3*basis.N),fixed)
    values,vectors,checks=modes(k[free,:][:,free],m[free,:][:,free])
    full=np.zeros((3*basis.N,6));full[free]=vectors
    samples=np.zeros((2*n+1,2*n+1,3,3,6))
    for component in range(3):
        for row in range(basis.N):
            x,y,z=basis.doflocs[:,row];ix=int(round(2*n*x));iy=int(round(2*n*(y+.5)));iz=int(round(2*z/H))+1
            samples[iy,ix,iz,component,:]=full[component*basis.N+row]
    coefficients=np.stack([samples[:,:,1],(samples[:,:,2]-samples[:,:,0])/H,
        (samples[:,:,2]+samples[:,:,0]-2*samples[:,:,1])/(H*H/4)],axis=2)
    record=dict(n=n,elements=n*n,total_coordinates=int(3*basis.N),free_coordinates=len(free),modes=checks,basis_component_crosscheck=crosscheck,
        frequencies_hz=(np.sqrt(values)/(2*np.pi)).tolist(),seconds=time.monotonic()-start,peak_process_working_bytes=plate.memory_bytes())
    record['eligible']=bool(checks['eligible'] and (crosscheck is None or crosscheck['eligible']) and record['seconds']<=240 and record['peak_process_working_bytes']<=4*1024**3)
    return record,dict(coefficients=coefficients,values=values)


def common_sampling():
    p,w=leggauss(3);points=[];weights=[]
    for j,i in itertools.product(range(48),repeat=2):
        for ix,iy,iz in itertools.product(range(3),repeat=3):
            points.append([(i+(p[ix]+1)/2)/48,-.5+(j+(p[iy]+1)/2)/48,H*p[iz]/2])
            weights.append(RHO*H*w[ix]*w[iy]*w[iz]/(8*48*48))
    return np.array(points),np.array(weights)


def sampled_modes(coefficients,n,points,weights):
    x,y,z=points.T;i=np.minimum((x*n).astype(int),n-1);j=np.minimum(((y+.5)*n).astype(int),n-1)
    xi=2*(x*n-i)-1;eta=2*((y+.5)*n-j)-1
    lx=[.5*xi*(xi-1),1-xi*xi,.5*xi*(xi+1)];ly=[.5*eta*(eta-1),1-eta*eta,.5*eta*(eta+1)]
    zz=[np.ones(len(z)),z,.5*z*z];result=np.zeros((len(z),3,6))
    for ox,oy in (common.nodes+1).astype(int):
        for kind in range(3):
            factor=lx[ox]*ly[oy]*zz[kind]
            result+=factor[:,None,None]*coefficients[2*j+oy,2*i+ox,kind]
    return (result*np.sqrt(weights)[:,None,None]).reshape(-1,6)


def match_modes(reference,candidate,ref_values,values):
    mac=np.abs(reference.T@candidate)**2
    rows,cols=linear_sum_assignment(-mac);mapping=np.empty(6,dtype=int);mapping[rows]=cols
    freq=np.sqrt(values[mapping])/(2*np.pi);ref_freq=np.sqrt(ref_values)/(2*np.pi)
    errors=np.abs(freq/ref_freq-1);groups=[];i=0
    while i<4:
        end=i+1
        while end<6 and (ref_freq[end]-ref_freq[end-1])/ref_freq[end-1]<=1e-3:end+=1
        a=reference[:,i:end];b=candidate[:,mapping[i:end]]
        a=a@np.linalg.inv(np.linalg.cholesky(a.T@a).T);b=b@np.linalg.inv(np.linalg.cholesky(b.T@b).T)
        singular=np.linalg.svd(a.T@b,compute_uv=False)
        groups.append(dict(first=i,last=end-1,max_sin_angle=float(np.sqrt(max(0.,1-float(singular.min())**2)))))
        i=end
    return dict(mapping=mapping.tolist(),matched_frequencies_hz=freq.tolist(),first4_frequency_relative=errors[:4].tolist(),
        matched_mac=mac[np.arange(6),mapping].tolist(),groups=groups)


def negative_controls(aux):
    q=aux['q'];ctx=aux['ctx'];cells=aux['cells'];free=aux['free'];mass=aux['mass'];d=np.tile(W9,len(q)//9)[free]
    scale=np.linalg.norm(d*aux['load'][free]);zeroctx=dict(ctx,e0=np.zeros_like(ctx['e0']))
    zero_force,_,_=evaluate(q,cells,zeroctx,need_k=False)
    reset=float(np.linalg.norm(d*(zero_force-aux['load'])[free])/scale)
    material_only=[]
    for cell in cells:material_only.append(material.analytic_evaluation(q[cell['dofs']],ctx)['kmat'])
    no_geo=assemble_matrices(material_only,cells,len(q))
    vmat,_,_=modes(no_geo[free,:][:,free],mass[free,:][:,free])
    omitted=float(np.max(np.abs(np.sqrt(vmat[:4]/aux['values'][:4])-1)))
    _,initial_k,_=evaluate(aux['q0'],cells,ctx)
    vi,_,_=modes(initial_k[free,:][:,free],mass[free,:][:,free])
    stale=float(np.max(np.abs(np.sqrt(vi[:4]/aux['values'][:4])-1)))
    return [dict(name='reset_initial_strain',mismatch=reset,detected=reset>.01),
        dict(name='omitted_geometric_tangent',mismatch=omitted,detected=omitted>1e-3),
        dict(name='reference_tangent_instead_of_equilibrium',mismatch=stale,detected=stale>1e-3)]


def main():
    start=time.monotonic();internal=[];external=[];aux=[];ext=[];stopped=None;comparisons=[];reference_comparison=None;negative=[];component_crosscheck=None
    source=Path(__file__).resolve()
    def checkpoint(stage):
        path=source.parents[2]/'build/analysis/q9-multielement-prestress-progress.json';path.parent.mkdir(parents=True,exist_ok=True)
        partial=dict(status='running',stage=stage,internal=internal,external=external,component_crosscheck=component_crosscheck,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
        path.write_text(json.dumps(partial,indent=2)+'\n',encoding='utf8')
    try:
        check,_=external_reference(2,verify_component_form=True);component_crosscheck=check['basis_component_crosscheck']
        checkpoint('component_contraction_crosscheck')
        if not check['eligible']:raise RuntimeError('component contraction mismatch')
        for n in [2,4,8,16]:
            record,data=run_internal(n);internal.append(record);aux.append(data)
            checkpoint('internal_'+str(n))
            print(json.dumps(dict(stage='internal',n=n,eligible=record['eligible'],seconds=record.get('seconds'),errors=record.get('static_errors'))),flush=True)
            if not record['eligible']:stopped='internal_failed';break
        if stopped is None:
            for n in [12,24]:
                record,data=external_reference(n);external.append(record);ext.append(data)
                checkpoint('external_'+str(n))
                print(json.dumps(dict(stage='external',n=n,eligible=record['eligible'],seconds=record['seconds'])),flush=True)
                if not record['eligible']:stopped='reference_failed';break
        if stopped is None:
            points,weights=common_sampling();profiles=sampled_modes(ext[-1]['coefficients'],24,points,weights)
            coarse=sampled_modes(ext[0]['coefficients'],12,points,weights)
            reference_comparison=match_modes(profiles,coarse,ext[-1]['values'],ext[0]['values'])
            for record,data in zip(internal,aux):
                field=sampled_modes(data['coefficients'],record['n'],points,weights)
                comparisons.append(dict(n=record['n'],**match_modes(profiles,field,ext[-1]['values'],data['values'])))
            negative=negative_controls(aux[-1])
    except Exception as exc:stopped=str(exc)
    checks={}
    if stopped is None:
        errors=internal[-1]['static_errors'];previous=internal[-2]['static_errors']
        order=float(np.log2(previous['displacement_l2']/errors['displacement_l2']))
        ref_error=np.array(reference_comparison['first4_frequency_relative']);internal_error=np.array(comparisons[-1]['first4_frequency_relative'])
        last_change=np.abs(np.array(comparisons[-1]['matched_frequencies_hz'][:4])/np.array(comparisons[-2]['matched_frequencies_hz'][:4])-1)
        checks=dict(static_finest=all(errors[k]<=.01 for k in ['displacement_l2','green_l2','energy_relative']),static_last_order=order,
            static_order_pass=bool(errors['displacement_l2']<=1e-8 or order>=1.5),reference_change_pass=bool(np.all(ref_error<=.002)),
            frequency_budget_pass=bool(np.all(ref_error+internal_error<=.004)),frequency_budget_indicators=(ref_error+internal_error).tolist(),
            last_frequency_change=last_change.tolist(),last_frequency_change_pass=bool(np.all(last_change<=.004)),
            finest_subspaces_pass=all(g['max_sin_angle']<=.05 for g in comparisons[-1]['groups']))
    elapsed=time.monotonic()-start
    passed=stopped is None and bool(component_crosscheck and component_crosscheck['eligible']) and len(internal)==4 and len(external)==2 and all(r['eligible'] for r in internal+external)
    passed=passed and all(checks.get(k,False) for k in ['static_finest','static_order_pass','reference_change_pass','frequency_budget_pass','last_frequency_change_pass','finest_subspaces_pass']) and len(negative)==3 and all(r['detected'] for r in negative) and elapsed<=900
    source=Path(__file__).resolve();modules=[common,curved,plate,interface,material,prestress]
    result=dict(protocol='docs/verification/q9-multielement-prestress-convergence.md',numpy=np.__version__,scipy=scipy.__version__,skfem=skfem.__version__,
        physical_inputs=dict(length=L,thickness=H,young=YOUNG,nu=NU,rho=RHO,initial_strain=E0.tolist(),a=A,b=B),
        component_crosscheck=component_crosscheck,internal=internal,external=external,mode_comparisons=comparisons,reference_comparison=reference_comparison,checks=checks,
        negative_controls=negative,seconds=elapsed,stopped=stopped,all_criteria_pass=bool(passed),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),dependencies=[dict(path=Path(m.__file__).name,sha256=hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest()) for m in modules],
        limitation='Flat isotropic manufactured nonuniform extension/shear, fixed thickness-quadratic space; no pure thin-bending, curve/laminate, wing, nonlinear transient or aeroelastic acceptance')
    output=source.parents[2]/'build/analysis/q9-multielement-prestress-result.json';output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(all_criteria_pass=bool(passed),stopped=stopped,seconds=elapsed,checks=checks)),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
