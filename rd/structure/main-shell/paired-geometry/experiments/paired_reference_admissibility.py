"""Conservative global reference-regularity certificate on a polynomial chart."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import math
import polynomial_reference_surface as geometry

def derivative(coeff,axis):
    out={}
    for key,value in coeff.items():
        if key[axis]:
            target=list(key);target[axis]-=1;out[tuple(target)]=value*key[axis]
    return out
def combine(a,b,factor):
    result=dict(a)
    for k,v in b.items():result[k]=result.get(k,0)+factor*v
    return {k:v for k,v in result.items() if v}
def certify(parameters,h,variable=False):
    if not math.isfinite(h) or h<=0:raise ValueError('nonpositive_reference_thickness')
    if not all(math.isfinite(x) for x in parameters):raise ValueError('nonfinite_reference_geometry')
    offset,width=geometry.chart(parameters);skew=geometry.shear(parameters);f=geometry.coefficients(parameters)
    def dx(c):return {k:v/width for k,v in derivative(c,0).items()}
    def dy(c):return combine(derivative(c,1),dx(c),-skew)
    bounds=[sum(abs(v) for v in c.values()) for c in [dx(dx(f)),dx(dy(f)),dy(dy(f))]]
    B=math.nextafter(math.sqrt(float(bounds[0]**2+2*bounds[1]**2+bounds[2]**2)),math.inf)
    Hmax=h*(1.15 if variable else 1.);theta=math.nextafter(Hmax*B/2,math.inf)
    if theta>=1:raise ValueError('uncertified_reference_normal_tube')
    lower=math.nextafter(h*float(width)/2*(1-theta)**2,0.)
    if lower<=0:raise ValueError('nonpositive_reference_measure_bound')
    return dict(hessian_frobenius_upper=B,max_half_thickness_curvature=theta,global_det_J_lower=lower,
        sufficient_not_necessary=True,domain='s,t in [0,1]; zeta in [-1,1]')
