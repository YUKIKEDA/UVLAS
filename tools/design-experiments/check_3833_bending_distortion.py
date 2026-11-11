"""Prescribed-field design experiment, not a product element or equilibrium solver.
Protocol: docs/verification/3833-bending-distortion-experiment.md.
"""
import json
import platform
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

nodes8=np.array([[-1,-1],[1,-1],[1,1],[-1,1],[0,-1],[1,0],[0,1],[-1,0]],float)
nodes9=np.vstack((nodes8,[[0,0]]))
young=1.; nu=.3; kappa=.001; distortion=.3
mu=young/(2*(1+nu)); lam=young*nu/((1+nu)*(1-2*nu))
elastic=np.diag([2*mu]*3+[mu]*3); elastic[:3,:3]+=lam

def lagrange(t,i):
    if i==-1:return .5*t*(t-1),t-.5
    if i==0:return 1-t*t,-2*t
    return .5*t*(t+1),t+.5

def shape(x,y,nnode):
    if nnode==9:
        n=[]; nx=[]; ny=[]
        for a,b in nodes9:
            la,da=lagrange(x,a); lb,db=lagrange(y,b)
            n.append(la*lb); nx.append(da*lb); ny.append(la*db)
        return np.array(n),np.array(nx),np.array(ny)
    n=np.zeros(8); nx=np.zeros(8); ny=np.zeros(8)
    for i,(a,b) in enumerate(nodes8[:4]):
        n[i]=.25*(1+a*x)*(1+b*y)*(a*x+b*y-1)
        nx[i]=a*.25*(1+b*y)*(2*a*x+b*y)
        ny[i]=b*.25*(1+a*x)*(a*x+2*b*y)
    n[4:]=[.5*(1-x*x)*(1-y),.5*(1+x)*(1-y*y),.5*(1-x*x)*(1+y),.5*(1-x)*(1-y*y)]
    nx[4:]=[-x*(1-y),.5*(1-y*y),-x*(1+y),-.5*(1-y*y)]
    ny[4:]=[-.5*(1-x*x),-y*(1+x),.5*(1-x*x),-y*(1-x)]
    return n,nx,ny

def jacobian(q,n,nx,ny,z,h):
    section=q[:,0,:]+z*q[:,1,:]+.5*z*z*q[:,2,:]
    dz=q[:,1,:]+z*q[:,2,:]
    return np.column_stack((nx@section,ny@section,(n@dz)*h/2))

def voigt(e):return np.array([e[0,0],e[1,1],e[2,2],2*e[1,2],2*e[0,2],2*e[0,1]])

def bending(nnode,geometry,h,order):
    nodes=nodes8 if nnode==8 else nodes9
    q0=np.zeros((nnode,3,3)); du=np.zeros_like(q0)
    for i,(xi,eta) in enumerate(nodes):
        x=.5*xi; y=.5*eta
        if geometry=='SKEW':x+=.2*eta
        if geometry=='TRAP':x*=1+distortion*eta
        q0[i,0]=[x,y,0]; q0[i,1]=[0,0,1]
        du[i,0]=[0,0,.5*kappa*(x*x-nu*y*y)]
        du[i,1]=[-kappa*x,nu*kappa*y,0]
        du[i,2]=[0,0,nu*kappa]
    rules=[leggauss(k) for k in order]
    total=shear=independent_shear=0.; maximum_shear=reference_shear_error=0.; minimum_det=np.inf
    for xi,wx in zip(*rules[0]):
        for eta,wy in zip(*rules[1]):
            n,nx,ny=shape(xi,eta,nnode)
            if max(abs(n.sum()-1),abs(nx.sum()),abs(ny.sum()))>1e-12:raise RuntimeError('Partition identity')
            for zeta,wz in zip(*rules[2]):
                z=h*zeta/2
                j0=jacobian(q0,n,nx,ny,z,h)
                det=float(np.linalg.det(j0)); minimum_det=min(minimum_det,det)
                if det<=0:raise RuntimeError('Invalid reference determinant')
                grad=np.linalg.solve(j0.T,jacobian(du,n,nx,ny,z,h).T).T
                strain=voigt(.5*(grad+grad.T))
                weight=wx*wy*wz*det
                total+=.5*float(strain@elastic@strain)*weight
                shear+=.5*mu*float(strain[3:]@strain[3:])*weight
                maximum_shear=max(maximum_shear,float(np.max(np.abs(strain[3:]))))
                if geometry=='TRAP' and nnode==8:
                    d=distortion
                    gx=kappa*d*d/(2*(1+d*eta))*xi*(1-eta*eta)
                    gy=kappa*d*d/2*(eta*(1-xi*xi)-d*xi*xi*(1-eta*eta)/(1+d*eta))
                    independent_shear+=.5*mu*(gx*gx+gy*gy)*wx*wy*wz*((1+d*eta)*h/8)
                    reference_shear_error=max(reference_shear_error,abs(strain[4]-gx),abs(strain[3]-gy),abs(strain[5]))
    exact=young*kappa*kappa*h**3/24
    return dict(nodes=nnode,geometry=geometry,h=h,quadrature=list(order),reference_bending_energy=exact,
                energy=total,energy_ratio=total/exact,shear_energy=shear,shear_to_bending=shear/exact,
                independent_shear_energy=independent_shear,independent_shear_error=float(reference_shear_error),
                maximum_shear=maximum_shear,minimum_det_j0=float(minimum_det))

def affine_curve(nnode):
    nodes=nodes8 if nnode==8 else nodes9
    q0=np.zeros((nnode,3,3))
    for i,(x,y) in enumerate(nodes):
        q0[i,0]=[x,y,.2*x*y]
        d=np.array([-.2*y,-.2*x,1.]);d/=np.linalg.norm(d)
        q0[i,1]=d;q0[i,2]=.2*d
    target=np.diag([np.sqrt(1.002),1.,1.]); q=q0@target.T
    q[:,0,:]+=[.3,-.4,.5]
    f_error=e_error=0.; det_min=np.inf; condition_max=0.
    points,weights=leggauss(5)
    for x in points:
        for y in points:
            n,nx,ny=shape(x,y,nnode)
            for zz in points:
                z=.01*zz
                j0=jacobian(q0,n,nx,ny,z,.02);j=jacobian(q,n,nx,ny,z,.02)
                det_min=min(det_min,float(np.linalg.det(j0)))
                if det_min<=0:raise RuntimeError('Invalid curved reference')
                condition_max=max(condition_max,float(np.linalg.cond(j0)))
                f=np.linalg.solve(j0.T,j.T).T
                e=.5*(f.T@f-np.eye(3))
                f_error=max(f_error,float(np.max(np.abs(f-target))))
                e_error=max(e_error,float(np.max(np.abs(e-np.diag([.001,0,0])))))
    n,nx,ny=shape(.5,.5,nnode);z=.004
    contribution=.5*z*z*(n@q0[:,2,:])
    full=n@(q0[:,0,:]+z*q0[:,1,:]+.5*z*z*q0[:,2,:])
    linear=n@(q0[:,0,:]+z*q0[:,1,:])
    c_error=float(np.max(np.abs((full-linear)-contribution)))
    return dict(nodes=nnode,max_F_error=f_error,max_E_error=e_error,minimum_det_j0=det_min,
                maximum_condition_j0=condition_max,quadratic_thickness_position_error=c_error,
                quadratic_thickness_contribution_norm=float(np.linalg.norm(contribution)))

checks={}
for nnode,nodes in [(8,nodes8),(9,nodes9)]:
    checks[f'nodal_{nnode}']=all(np.max(np.abs(shape(x,y,nnode)[0]-np.eye(nnode)[i]))<=1e-12 for i,(x,y) in enumerate(nodes))
results=[]
for nnode in [8,9]:
    for geometry in ['RECT','SKEW','TRAP']:
        for h in [1.,.1,.01,.001]:
            for order in [(3,3,2),(5,5,5),(7,7,7)]:
                r=bending(nnode,geometry,h,order);results.append(r)
                key=f'{nnode}_{geometry}_{h}_{order}'
                if geometry=='TRAP' and nnode==8:
                    checks[key]=bool(r['independent_shear_error']<=1e-14 and
                        abs(r['shear_energy']/r['independent_shear_energy']-1)<=1e-9 and
                        abs((r['energy']-r['reference_bending_energy'])/r['independent_shear_energy']-1)<=1e-9)
                else:
                    checks[key]=bool(abs(r['energy_ratio']-1)<=1e-10 and r['maximum_shear']<=1e-14)
curved=[affine_curve(n) for n in [8,9]]
for r in curved:
    checks[f'curve_{r["nodes"]}']=bool(r['max_F_error']<=1e-12 and r['max_E_error']<=1e-12 and
        r['quadratic_thickness_position_error']<=1e-12 and r['quadratic_thickness_contribution_norm']>1e-8)

report=dict(protocol='docs/verification/3833-bending-distortion-experiment.md',
            python=platform.python_version(),numpy=np.__version__,platform=platform.platform(),
            inputs=dict(a=1.,b=1.,Y=young,nu=nu,kappa=kappa,distortion=distortion,skew=.2),
            all_protocol_checks_pass=all(checks.values()),checks=checks,bending=results,curved=curved,
            limitation='Prescribed-field reproduction; not an equilibrium or mesh-convergence solution')
out=Path(__file__).resolve().parents[2]/'build'/'analysis'/'3833-bending-distortion-result.json'
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
summary=[dict(nodes=r['nodes'],h=r['h'],ratio=r['energy_ratio'],shear_ratio=r['shear_to_bending'])
         for r in results if r['geometry']=='TRAP' and r['quadrature']==[7,7,7]]
print(json.dumps(dict(all_protocol_checks_pass=report['all_protocol_checks_pass'],count=len(checks),
                     failed=[k for k,v in checks.items() if not v],trap_7_point=summary,curved=curved),indent=2))
if not all(checks.values()):raise SystemExit(1)
