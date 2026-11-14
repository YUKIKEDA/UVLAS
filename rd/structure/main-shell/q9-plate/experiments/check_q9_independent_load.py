"""External linear FEM comparison, not a UVLAS element implementation.
Protocol: rd/structure/main-shell/q9-plate/records/q9-independent-load-experiment.md.
Environment: requirements-reference.txt.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.sparse import csr_matrix
import skfem
from skfem import MeshHex, ElementHex2, ElementVector, Basis, FacetBasis, LinearForm, BilinearForm, asm
from skfem.helpers import dot
from skfem.models.elasticity import linear_elasticity, lame_parameters
from numpy.polynomial.legendre import leggauss

nodes=np.array([[-1,-1],[1,-1],[1,1],[-1,1],[0,-1],[1,0],[0,1],[-1,0],[0,0]],float)
young=1.;nu=.3;rho=1.;kappa=.001
mu=young/(2*(1+nu));lam=young*nu/((1+nu)*(1-2*nu))
elastic=np.diag([2*mu]*3+[mu]*3);elastic[:3,:3]+=lam
chol=np.linalg.cholesky(elastic).T
points,weights=leggauss(5)

def lag(t,i):
    if i==-1:return .5*t*(t-1),t-.5
    if i==0:return 1-t*t,-2*t
    return .5*t*(t+1),t+.5

def shape(x,y):
    n=[];nx=[];ny=[]
    for a,b in nodes:
        la,da=lag(x,a);lb,db=lag(y,b)
        n.append(la*lb);nx.append(da*lb);ny.append(la*db)
    return np.array(n),np.array(nx),np.array(ny)

def kinematics(x,y,z):
    n,nx,ny=shape(x,y)
    scalar=np.column_stack((n,z*n,.5*z*z*n)).reshape(27)
    s=np.kron(scalar.reshape(1,27),np.eye(3))
    dx=np.column_stack((nx,z*nx,.5*z*z*nx)).reshape(27)*2
    dy=np.column_stack((ny,z*ny,.5*z*z*ny)).reshape(27)*2
    dz=np.column_stack((np.zeros(9),n,z*n)).reshape(27)
    b=np.zeros((6,81))
    for i,(a,c,d) in enumerate(zip(dx,dy,dz)):
        b[:,3*i]=[a,0,0,0,d,c]
        b[:,3*i+1]=[0,c,0,d,0,a]
        b[:,3*i+2]=[0,0,d,c,a,0]
    return s,b

def assemble_own(h,order=(5,5,5)):
    mr=[];br=[]
    rules=[leggauss(n) for n in order]
    for x,wx in zip(*rules[0]):
        for y,wy in zip(*rules[1]):
            for zz,wz in zip(*rules[2]):
                s,b=kinematics(x,y,h*zz/2)
                w=wx*wy*wz*h/8
                mr.append(np.sqrt(rho*w)*s);br.append(np.sqrt(w)*chol@b)
    m=np.vstack(mr);b=np.vstack(br)
    return m.T@m,b.T@b,m,b

def own_load(h):
    f=np.zeros(81)
    for face in [-1,1]:
        for y,wy in zip(points,weights):
            for zz,wz in zip(points,weights):
                z=h*zz/2;s,_=kinematics(face,y,z)
                f+=s.T@np.array([-young*kappa*z*face,0,0])*(wy*wz*h/4)
    return f

def exact_displacement(p):
    x,y,z=p
    u=np.array([-kappa*x*z,nu*kappa*y*z,.5*kappa*(x*x-nu*y*y+nu*z*z)])
    return u+np.cross(np.array([nu*kappa/4,kappa/4,0.]),p)

def exact_q():
    q=np.zeros((9,3,3));omega=np.array([nu*kappa/4,kappa/4,0.])
    for i,(xi,eta) in enumerate(nodes):
        x=xi/2;y=eta/2
        q[i,0]=exact_displacement(np.array([x,y,0.]))
        q[i,1]=[-kappa*x,nu*kappa*y,0.];q[i,1]+=np.cross(omega,[0,0,1.])
        q[i,2]=[0,0,nu*kappa]
    return q.reshape(81)

def rank(a):
    s=np.linalg.svd(a,compute_uv=False)
    return {str(t):int(np.count_nonzero(s>t*s[0])) for t in [1e-11,1e-10,1e-9]}

def matrix_error(a,b):
    scale=1/np.sqrt(np.diag(a))
    reference=scale[:,None]*a*scale[None,:]
    difference=scale[:,None]*(a-b)*scale[None,:]
    return dict(relative=float(np.linalg.norm(difference)/np.linalg.norm(reference)),maximum=float(np.max(np.abs(difference))))

def force_error(k,a,b):
    scale=1/np.sqrt(np.diag(k))
    return float(np.linalg.norm(scale*(a-b))/np.linalg.norm(scale*a))

def solve_scaled(k,f,fixed,external=False):
    free=np.setdiff1d(np.arange(len(f)),fixed)
    kff=k[np.ix_(free,free)];scale=1/np.sqrt(np.diag(kff))
    a=scale[:,None]*kff*scale[None,:];rhs=scale*f[free]
    cond=float(np.linalg.cond(a))
    if cond>1e12:return dict(status='hold_condition',condition=cond),None
    y=skfem.solve(csr_matrix(a),rhs) if external else np.linalg.solve(a,rhs)
    u=np.zeros(len(f));u[free]=scale*y
    residual=float(np.linalg.norm(a@y-rhs)/np.linalg.norm(rhs))
    return dict(status='evaluated',condition=cond,residual=residual),u

@BilinearForm
def mass(u,v,w):return rho*dot(u,v)

@LinearForm
def surface(v,w):return -young*kappa*w.x[2]*w.n[0]*v[0]

fixed_q=np.array([8*9,8*9+1,8*9+2,5*9+1,5*9+2,6*9+2])
q_reference=exact_q()
def main():
    cases=[]
    for h in [.1,.01,.001]:
        mesh=MeshHex.init_tensor(np.array([-.5,.5]),np.array([-.5,.5]),np.array([-h/2,h/2]))
        elem=ElementVector(ElementHex2())
        basis=Basis(mesh,elem,intorder=9)
        if basis.N!=81 or basis.X.shape[1]!=125:raise RuntimeError('Unexpected external element or quadrature')
        kp=asm(linear_elasticity(*lame_parameters(young,nu)),basis).toarray()
        mp=asm(mass,basis).toarray()
        faces=mesh.facets_satisfying(lambda x:np.isclose(np.abs(x[0]),.5))
        fb=FacetBasis(mesh,elem,facets=faces,intorder=9)
        fp=asm(surface,fb)
        transform=np.zeros((81,81));p_reference=np.zeros(81)
        for comp,ids in enumerate(basis.split_indices()):
            for row in ids:
                xyz=basis.doflocs[:,row];x,y,z=xyz
                z_index=int(round(2*z/h))
                if z_index not in [-1,0,1] or abs(2*z/h-z_index)>1e-12:
                    raise RuntimeError('External thickness node mismatch')
                z=h*z_index/2
                matches=np.flatnonzero(np.max(np.abs(nodes-[2*x,2*y]),axis=1)<1e-12)
                if len(matches)!=1:raise RuntimeError('External node correspondence ambiguous')
                i=int(matches[0]);transform[row,[9*i+comp,9*i+3+comp,9*i+6+comp]]=[1,z,.5*z*z]
                p_reference[row]=exact_displacement(xyz)[comp]
        if np.max(np.abs(transform@q_reference-p_reference))>1e-14:raise RuntimeError('Analytic field and T disagree')
        mq,kq,m_observe,b_observe=assemble_own(h)
        fq=own_load(h)
        kt=transform.T@kp@transform;mt=transform.T@mp@transform;ft=transform.T@fp
        errors=dict(stiffness=matrix_error(kq,kt),mass=matrix_error(mq,mt),force=force_error(kq,fq,ft))
        fixed_p=[]
        for qcol in fixed_q:
            rows=np.flatnonzero(np.abs(transform[:,qcol])>0)
            exact_rows=[r for r in rows if np.count_nonzero(transform[r])==1]
            if len(exact_rows)!=1:raise RuntimeError('Gauge correspondence ambiguous')
            fixed_p.append(exact_rows[0])
        rigid=np.zeros((81,6))
        for i,(xi,eta) in enumerate(nodes):
            rigid[9*i:9*i+3,:3]=np.eye(3)
            for j,omega in enumerate(np.eye(3)):
                rigid[9*i:9*i+3,3+j]=np.cross(omega,[xi/2,eta/2,0])
                rigid[9*i+3:9*i+6,3+j]=np.cross(omega,[0,0,1.])
        rigid_work=float(np.linalg.norm(rigid.T@fq)/(np.linalg.norm(rigid)*np.linalg.norm(fq)))
        record=dict(h=h,comparison=errors,rigid_load_work=rigid_work,solutions={})
        if h==.1:
            _,_,_,reduced=assemble_own(h,(3,3,2))
            record['ranks']=dict(mass=rank(m_observe),stiffness=rank(b_observe),stiffness_reduced=rank(reduced))
        baseline_ok=all(v<=1e-10 for e in [errors['stiffness'],errors['mass']] for v in e.values()) and errors['force']<=1e-10 and rigid_work<=1e-12
        record['comparison_pass']=baseline_ok
        if not baseline_ok:
            record['solutions']['stopped']='matrix_or_force_comparison_failed'
            cases.append(record);continue
        for name,k,f,fixed,is_external,t in [('own_q',kq,fq,fixed_q,False,transform),('external_p',kp,fp,np.array(fixed_p),True,np.eye(81)),('external_in_q',kt,ft,fixed_q,True,transform)]:
            result,u=solve_scaled(k,f,fixed,is_external)
            if u is not None:
                p=t@u
                q=u if name!='external_p' else np.linalg.solve(transform,u)
                displacement=float(np.max(np.abs(p-p_reference))/np.max(np.abs(p_reference)))
                reference_energy=young*kappa*kappa*h**3/24
                energy=float(.5*u@k@u)
                strain_error=0.
                for x in points:
                    for y in points:
                        for zz in points:
                            z=h*zz/2;_,b=kinematics(x,y,z)
                            strain_error=max(strain_error,float(np.max(np.abs(b@q-np.array([-kappa*z,nu*kappa*z,nu*kappa*z,0,0,0])))))
                result.update(displacement_relative=displacement,energy_relative=abs(energy/reference_energy-1),
                              strain_error=strain_error,normalized_strain_error=strain_error/(kappa*h))
                result['criteria_pass']=bool(result['residual']<=1e-10 and displacement<=1e-6 and result['energy_relative']<=1e-6 and strain_error<=1e-5*kappa*h)
            record['solutions'][name]=result
        cases.append(record)

    report=dict(protocol='rd/structure/main-shell/q9-plate/records/q9-independent-load-experiment.md',python=platform.python_version(),
                numpy=np.__version__,scipy=scipy.__version__,skfem=skfem.__version__,platform=platform.platform(),
                inputs=dict(a=1.,b=1.,Y=young,nu=nu,rho=rho,kappa=kappa),cases=cases,
                limitation='Linear homogeneous single element; no nonlinear, laminate, curved or convergence evidence')
    out=rd_paths.ROOT/'build'/'analysis'/'q9-independent-load-result.json'
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
