"""IL-04: published external basis transformed before integration.
Single element diagnostic only; not a product multi-element basis.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import json
import platform
import sys
from decimal import Decimal, localcontext
from pathlib import Path
import numpy as np
import scipy
import skfem
from skfem import MeshHex,ElementHex2,ElementVector,Basis,FacetBasis,asm
from skfem.models.elasticity import linear_elasticity,lame_parameters
import check_q9_independent_load as common

use_refinement='--refine' in sys.argv

def refined_solve(k,f,fixed,external):
    free=np.setdiff1d(np.arange(len(f)),fixed)
    kff=k[np.ix_(free,free)];scale=1/np.sqrt(np.diag(kff))
    a=scale[:,None]*kff*scale[None,:];rhs=scale*f[free]
    condition=float(np.linalg.cond(a))
    if condition>1e12:return dict(status='hold_condition',condition=condition),None
    def solve(r):
        return skfem.solve(common.csr_matrix(a),r) if external else np.linalg.solve(a,r)
    y=solve(rhs);history=[]
    for step in range(4):
        residual64=rhs-a@y
        with localcontext() as ctx:
            ctx.prec=80
            yd=[Decimal.from_float(float(v)) for v in y]
            residual80=np.array([float(Decimal.from_float(float(b))-
                sum((Decimal.from_float(float(v))*z for v,z in zip(row,yd)),Decimal(0)))
                for row,b in zip(a,rhs)])
        r64=float(np.linalg.norm(residual64)/np.linalg.norm(rhs))
        r80=float(np.linalg.norm(residual80)/np.linalg.norm(rhs))
        history.append(dict(corrections=step,residual64=r64,residual80=r80))
        if r64<=1e-10 and r80<=1e-10:break
        if step<3:y+=solve(residual80)
    q=np.zeros(len(f));q[free]=scale*y
    return dict(status='evaluated',condition=condition,residual=r64,residual_high_precision_sum=r80,
                corrections=len(history)-1,refinement_history=history),q

class BeforeIntegration(ElementHex2):
    def __init__(self,coefficients):
        self.coefficients=coefficients

    def lbasis(self,x,i):
        value=gradient=None
        for k in np.flatnonzero(self.coefficients[:,i]):
            v,d=super().lbasis(x,int(k))
            a=self.coefficients[k,i]
            value=a*v if value is None else value+a*v
            gradient=a*d if gradient is None else gradient+a*d
        return value,gradient

cases=[]
for h in [.1,.01,.001]:
    mesh=MeshHex.init_tensor(np.array([-.5,.5]),np.array([-.5,.5]),np.array([-h/2,h/2]))
    ordinary=Basis(mesh,ElementVector(ElementHex2()),intorder=9)
    xyz=ordinary.doflocs[:,ordinary.element_dofs[::3,0]].T
    coefficient=np.zeros((27,27));semantic=[]
    for j,(x,y,z) in enumerate(xyz):
        plane=np.flatnonzero(np.max(np.abs(common.nodes-[2*x,2*y]),axis=1)<1e-12)
        level=int(round(2*z/h))
        if len(plane)!=1 or level not in [-1,0,1] or abs(2*z/h-level)>1e-12:
            raise RuntimeError('External local reference-node correspondence failed')
        kind=level+1
        semantic.append((int(plane[0]),kind))
        for k,(xx,yy,zz) in enumerate(xyz):
            if max(abs(x-xx),abs(y-yy))<1e-12:
                zk=h*round(2*zz/h)/2
                coefficient[k,j]=[1.,zk,.5*zk*zk][kind]
    external_element=ElementVector(BeforeIntegration(coefficient))
    basis=Basis(mesh,external_element,intorder=9)
    permutation=np.zeros((81,81))
    for j,(plane,kind) in enumerate(semantic):
        for comp in range(3):
            permutation[basis.element_dofs[3*j+comp,0],9*plane+3*kind+comp]=1
    if not (np.all(permutation.sum(0)==1) and np.all(permutation.sum(1)==1)):
        raise RuntimeError('Semantic permutation is not bijective')
    kr=permutation.T@asm(linear_elasticity(*lame_parameters(common.young,common.nu)),basis).toarray()@permutation
    mr=permutation.T@asm(common.mass,basis).toarray()@permutation
    faces=mesh.facets_satisfying(lambda x:np.isclose(np.abs(x[0]),.5))
    fb=FacetBasis(mesh,external_element,facets=faces,intorder=9)
    fr=permutation.T@asm(common.surface,fb)
    mq,kq,_,_=common.assemble_own(h);fq=common.own_load(h)
    errors=dict(stiffness=common.matrix_error(kq,kr),mass=common.matrix_error(mq,mr),force=common.force_error(kq,fq,fr))
    comparison=all(v<=1e-10 for e in [errors['stiffness'],errors['mass']] for v in e.values()) and errors['force']<=1e-10
    record=dict(h=h,comparison=errors,comparison_pass=comparison,solutions={})
    if comparison:
        for name,k,f,external in [('own_q',kq,fq,False),('external_preintegral_q',kr,fr,True)]:
            result,q=(refined_solve(k,f,common.fixed_q,external) if use_refinement
                      else common.solve_scaled(k,f,common.fixed_q,external))
            if q is not None:
                point_error=0.;point_ref=0.
                for x,y in common.nodes:
                    for z in [-h/2,0,h/2]:
                        s,_=common.kinematics(x,y,z)
                        reference=common.exact_displacement(np.array([x/2,y/2,z]))
                        point_error=max(point_error,float(np.max(np.abs(s@q-reference))))
                        point_ref=max(point_ref,float(np.max(np.abs(reference))))
                strain_error=0.
                for x in common.points:
                    for y in common.points:
                        for zz in common.points:
                            z=h*zz/2;_,b=common.kinematics(x,y,z)
                            reference=np.array([-common.kappa*z,common.nu*common.kappa*z,common.nu*common.kappa*z,0,0,0])
                            strain_error=max(strain_error,float(np.max(np.abs(b@q-reference))))
                energy=float(.5*q@k@q);reference_energy=common.young*common.kappa**2*h**3/24
                result.update(displacement_relative=point_error/point_ref,energy_relative=abs(energy/reference_energy-1),
                              strain_error=strain_error,normalized_strain_error=strain_error/(common.kappa*h))
                result['criteria_pass']=bool(result['residual']<=1e-10 and
                    result.get('residual_high_precision_sum',0)<=1e-10 and result['displacement_relative']<=1e-6 and
                    result['energy_relative']<=1e-6 and result['normalized_strain_error']<=1e-5)
            record['solutions'][name]=result
    else:record['solutions']['stopped']='matrix_or_force_comparison_failed'
    cases.append(record)

passed=all(c['comparison_pass'] and all(r.get('criteria_pass',False) for r in c['solutions'].values()) for c in cases)
report=dict(protocol='rd/structure/main-shell/q9-plate/records/q9-independent-load-experiment.md#'+('il-05' if use_refinement else 'il-04'),
            python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,skfem=skfem.__version__,
            all_criteria_pass=passed,cases=cases,
            independence='External ElementHex2 shape/gradient and external constitutive/forms/assembly/solver; only constant basis transform is custom',
            limitation='Single homogeneous linear element, one bending load; not mesh convergence or a product basis')
out=rd_paths.ROOT/'build'/'analysis'/('q9-refined-reference-result.json' if use_refinement else 'q9-preintegral-reference-result.json')
out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
print(json.dumps(report,indent=2))
if not passed:raise SystemExit(1)
