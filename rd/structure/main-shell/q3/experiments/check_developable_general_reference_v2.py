"""Kinematics, Decimal energy, analytic body loads and continuum weak work."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from decimal import Decimal,localcontext
import hashlib,itertools,json,time
import numpy as np
import developable_general_reference as ref
import developable_general_decimal_work as work
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def decimal_energy(kind,k,h):
    with localcontext() as ctx:
        ctx.prec=80;k=Decimal(str(k));h=Decimal(str(h));lam,mu=ref.constants(Decimal);nodes,weights=ref.decimal_gauss(32);integral=Decimal(0);terms=0
        for p,w in zip(nodes,weights):
            x=(p+1)/2;D=1+k*k*x*x;xi=k*h/(2*D);power=Decimal(1);series=Decimal(0)
            for q in range(100):
                series+=power/(2*q+3);power*=xi*xi
                if power<Decimal('1e-85'):terms=max(terms,q+1);break
            v=Decimal(1) if kind=='quadratic' else x*x;integral+=w/2*v*v*D*h**3/4*series
        return (lam+2*mu)/2*integral,terms
def main():
    dest=ROOT/'rd/structure/main-shell/q3/records/developable-general-reference-v2-result.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();records=[];stop=None
    try:
        for kind,k,h in itertools.product(['quadratic','quartic'],[0.,.4],[.01,.001,.0001]):
            points=np.array([[x,0.,z] for x,z in itertools.product(np.linspace(0,1,17),np.linspace(-h/2,h/2,5))]);x,y,z=points.T;J,E=ref.direct(kind,k,points);closed=ref.strain(kind,k,x,z);field=float(np.max(abs(E-closed)));cov=J.transpose(0,2,1)@E@J;shear=float(np.max(abs(cov[:,[0,1],2])));thickness=float(np.max(abs(cov[:,2,2])));mid=float(np.max(abs(E[z==0])));_,clamp=ref.positions(kind,k,np.zeros(5),np.zeros(5),np.linspace(-h/2,h/2,5));fixed=float(np.max(abs(clamp)));b=ref.body(kind,k,x,z);Sx=ref.stress(kind,k,x+1e-25j,z).imag/1e-25;Sz=ref.stress(kind,k,x,z+1e-25j).imag/1e-25;Ji,_=ref.inverse(k,x,z);other=-np.einsum('pij,pj->pi',Sx,Ji[:,0,:])-np.einsum('pij,pj->pi',Sz,Ji[:,2,:]);body=float(np.max(abs(b-other))/max(1.,float(np.max(abs(b)))))
            xp,xw=np.polynomial.legendre.leggauss(24);zp,zw=np.polynomial.legendre.leggauss(12);p=np.array([[(xx+1)/2,0.,zz*h/2] for xx,zz in itertools.product(xp,zp)]);w=np.array([a*b*h/4 for a,b in itertools.product(xw,zw)]);Jq,Eq=ref.direct(kind,k,p);energy=float(w@(np.linalg.det(Jq)*ref.energy(Eq))/2);de,terms=decimal_energy(kind,k,h);relative=float(abs(Decimal.from_float(energy)/de-1));lam,mu=ref.constants();flat=abs(energy/((lam+2*mu)*h**3/(24 if kind=='quadratic' else 120))-1) if k==0 else 0.
            low=work.weak(kind,k,h,(24,4,8));high=work.weak(kind,k,h,(32,5,12));delta=float(np.sqrt(np.sum((low['internal']-high['internal'])**2)/np.sum(high['internal']**2)));weak={key:high[key] for key in ['relative','cauchy','allowed_reaction','reference']};checks=dict(field=field<=1e-11,kinematic=max(shear,thickness,mid,fixed)<=1e-11,geometry=1-abs(k)*h/2>0,energy=relative<=1e-10,flat=flat<=1e-10,body=body<=1e-11,weak=high['relative']<=1e-8,cauchy=high['cauchy']<=1e-10,quadrature=delta<=1e-8,reaction=high['allowed_reaction']==0,resources=time.monotonic()-start<=180 and memory.memory_bytes()<=1024**3)
            records.append(dict(kind=kind,kappa=k,h=h,field_error=field,covariant_shear=shear,thickness_strain=thickness,midstrain=mid,clamp=fixed,body_difference=body,energy=energy,decimal_energy=str(de),energy_difference=relative,series_terms=terms,weak=weak,quadrature=delta,checks={key:bool(v) for key,v in checks.items()}));print(kind,k,h,records[-1]['checks'],high['relative'],flush=True)
            if not all(checks.values()):stop='required_gate_failed';break
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(ref.__file__),Path(work.__file__),ROOT/'rd/structure/main-shell/q3/records/developable-general-reference-decimal-work-protocol.md',ROOT/'rd/structure/main-shell/q3/records/developable-general-reference-result.json',ROOT/'rd/structure/main-shell/q3/developable-general-bending-reference-contract.md',ROOT/'rd/structure/main-shell/q3/records/developable-general-reference-experiment.md']
    dest.write_text(json.dumps(dict(records=records,all_checks_pass=len(records)==12 and stop is None,stop_reason=stop,longdouble_mantissa_bits=np.finfo(ref.LD).nmant+1,weak_work_decimal_digits=80,physical_adoption=False,discrete_equilibrium_solved=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
