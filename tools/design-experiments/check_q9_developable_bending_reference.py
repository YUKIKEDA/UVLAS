"""Analytic developable bending reference versus independent 3D and Decimal energy."""
from pathlib import Path
from decimal import Decimal,localcontext
import json,itertools,time
import numpy as np
from numpy.polynomial.legendre import leggauss
import q9_curved_energy_scaling as material
core=material.core;ROOT=core.ROOT;F=core.F
def dec(v):return Decimal.from_float(float(v))
def tables(k,x,z):
    J=np.array([[1-k*z,0.,-k*x],[0.,1.,0.],[k*x,0.,1.]])
    Q=np.array([[-k*x*x-z,0.,-x],[0.,0.,0.],[x-2*k*x*z,0.,-k*x*x]])
    grad=Q@np.linalg.inv(J);E=(grad+grad.T)/2;a=np.array([1.,0.,k*x]);closed=-z*np.outer(a,a)/(1+k*k*x*x-k*z)
    return J,E,closed
def main():
    target=F/'q9-developable-bending-reference-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();cases=[]
    for k,h in itertools.product([0.,.4],[.01,.001,.0001]):
        error=0.;shear=0.;middle=0.
        for x,z in itertools.product([0.,.13,.49,.83,1.],[-h/2,0.,h/2]):
            J,E,C=tables(k,x,z);error=max(error,float(np.linalg.norm(E-C)));cov=J.T@E@J;shear=max(shear,float(np.linalg.norm(cov[[0,1],2])));middle=max(middle,float(np.linalg.norm(E)) if z==0 else 0.)
        xp,xw=leggauss(32);zp,zw=leggauss(12);energy=0.
        for i,j in itertools.product(range(32),range(12)):
            x=(xp[i]+1)/2;z=h*zp[j]/2;J,E,_=tables(k,x,z);e=np.einsum('aij,ij->a',core.method.MB,E);energy+=xw[i]/2*zw[j]*h/2*np.linalg.det(J)*(e@material.CM@e)/2
        with localcontext() as ctx:
            ctx.prec=80;nodes,weights=leggauss(36);kd,hd=dec(k),dec(h);integral=Decimal(0);max_terms=0
            for t,w in zip(nodes,weights):
                x=(dec(t)+1)/2;D=1+kd*kd*x*x;xi=kd*hd/(2*D);power=Decimal(1);series=Decimal(0)
                for q in range(100):
                    term=power/Decimal(2*q+3);series+=term;power*=xi*xi
                    if power<Decimal('1e-85'):max_terms=max(max_terms,q+1);break
                integral+=dec(w)/2*D*hd**3/4*series
            reference=(dec(core.shapes.lam)+2*dec(material.MU))/2*integral;relative=float(abs(dec(energy)/reference-1));flat=abs(energy/((core.shapes.lam+2*material.MU)*h**3/24)-1) if k==0 else None
        direction=max(abs(np.array([-k*x,0.,1.])@np.array([-x,0.,-k*x*x])) for x in [0.,.13,.49,.83,1.]);clamp=np.linalg.norm(np.array([-k*0**3/3,0.,0**2/2]));lower=1-abs(k)*h/2;checks=dict(field=error<=1e-11,shear=shear<=1e-11,membrane=middle<=1e-11,direction=direction<=1e-11,clamp=clamp==0,geometry=lower>0,energy=relative<=1e-11,flat=flat is None or flat<=1e-11)
        cases.append(dict(kappa=k,h=h,field_error=error,covariant_shear=shear,midstrain=middle,direction_variation=float(direction),jacobian_lower=lower,energy=energy,decimal_energy=str(reference),energy_relative_difference=relative,series_terms=max_terms,energy_over_h3=energy/h**3,checks={key:bool(value) for key,value in checks.items()}))
    result=dict(source_sha256=core.digest(__file__),definition_sha256=core.digest(ROOT/'docs/q9-developable-bending-reference.md'),protocol_sha256=core.digest(F/'q9-developable-bending-reference-protocol.md'),dependencies={Path(m.__file__).name:core.digest(m.__file__) for m in [material,core,core.shapes]},cases=cases,seconds=time.monotonic()-start,peak_bytes=core.memory.memory_bytes(),all_criteria_pass=bool(len(cases)==6 and all(all(c['checks'].values()) for c in cases)),product_shell_adopted=False,load_equilibrium_reference_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
