"""Exact input-geometry integration and the full-to-half dual observation map."""
from bisect import bisect_left,bisect_right
from fractions import Fraction
import math
import numpy as np
import avl_field_overlay_exact as geo
import check_avl_conservative_observation as obs

LAYOUT=['tau0','tau1','kbar','jbar','jcbar']
def public_vector(bands):return np.array([float(b[k]) for k in LAYOUT for b in bands])
def folded_order(m,native):
    if native%m:raise ValueError('shape')
    half=native//m
    return np.array([(half-1-s)*m+j for s in range(half) for j in range(m)])
def weights(rows,points,case,m):
    n=len(rows)
    if rows.shape!=(n,21) or points.shape!=(n,3) or n%(2*m) or not np.all(np.isfinite(rows)) or not np.all(np.isfinite(points)):raise ValueError('geometry_shape_or_finiteness')
    a=rows[:,7:10];b=rows[:,10:13];force=2*np.cross(np.array([1.,0.,0.]),b-a);couples=obs.wrench.convert(a,b,points,force)[:,1]
    domain=geo.wing_domain(case);bounds=obs.BOUNDS;matrix=np.zeros((100,n))
    for i,(aa,bb) in enumerate(zip(a,b)):
        low=geo.D(aa[1]);high=geo.D(bb[1]);length=high-low
        if length<=0:raise ValueError('nonpositive_span')
        front=geo.line(aa,bb);density=geo.D(couples[i])/length
        for j in range(max(0,bisect_right(bounds,low)-1),min(20,bisect_left(bounds,high))):
            lo=max(low,bounds[j]);hi=min(high,bounds[j+1])
            if hi<=lo:continue
            h=bounds[j+1]-bounds[j];center=(bounds[j]+bounds[j+1])/2;te=domain.piece(center)[3]
            t0=hi-lo;t1=((hi-center)**2-(lo-center)**2)/2
            kval=obs.ip((te[0]-front[0],te[1]-front[1]),lo,hi)
            line=obs.ip((front[0]-obs.OX,front[1]),lo,hi)
            for k,value in enumerate([t0/h,6*t1/(h*h),kval/h,line/h,(line-density*t0/2)/h]):matrix[20*k+j,i]=float(value)
    native=n//2;order=folded_order(m,native);reduced=matrix[:,native:].copy();reduced[:,order]+=matrix[:,:native]
    return matrix,reduced,order
def direct(rows,points,g,case,m):
    domain=geo.wing_domain(case);a=rows[:,7:10];b=rows[:,10:13];force=2*g[:,None]*np.cross(np.array([1.,0.,0.]),b-a)
    couples=obs.wrench.convert(a,b,points,force)[:,1];strips=[];densities=[]
    for i in range(0,len(rows),m):
        low=geo.D(a[i,1]);high=geo.D(b[i,1]);fronts=tuple(geo.line(aa,bb) for aa,bb in zip(a[i:i+m],b[i:i+m]))
        strips.append(geo.Strip(low,high,fronts,tuple(float(v) for v in g[i:i+m])))
        densities.append(sum((geo.D(v) for v in couples[i:i+m]),Fraction(0))/(high-low))
    field=geo.Field(domain,tuple(strips));return public_vector(obs.observe(obs.profile(field,tuple(densities))))
def channels(vector,widths):
    root=np.sqrt(widths)
    return dict(tau0=root*vector[:20],tau_p1=np.concatenate([root*vector[:20],np.sqrt(widths/3)*vector[20:40]]),
        k0=root*vector[40:60],j_line=root*vector[60:80],j_point=root*vector[80:100])
def totals(vector,widths):
    return dict(lift=2*math.fsum(vector[:20]*widths)/obs.S,
        moment_line=-2*math.fsum(vector[60:80]*widths)/(obs.S*obs.C),moment_point=-2*math.fsum(vector[80:100]*widths)/(obs.S*obs.C))
def apply(matrix,q):return np.array([math.fsum(float(a)*float(b) for a,b in zip(row,q)) for row in matrix])
def manufactured_geometry():
    rows=np.zeros((2,21));rows[0,7:10]=[2.975,-4,0];rows[0,10:13]=[.375,0,0];rows[1,7:10]=[.375,0,0];rows[1,10:13]=[2.975,4,0]
    # Keep points on their bound segments, with mirrored parameters.
    points=np.array([rows[0,7:10]+(.5-.1/2.6)*(rows[0,10:13]-rows[0,7:10]),rows[1,7:10]+(.5+.1/2.6)*(rows[1,10:13]-rows[1,7:10])])
    full,reduced,order=weights(rows,points,'FW-V01',1);value=apply(reduced,np.array([1.]))
    expected=[];domain=geo.wing_domain('FW-V01')
    for k in LAYOUT:
        for lo,hi in zip(obs.BOUNDS[:-1],obs.BOUNDS[1:]):
            center=float((lo+hi)/2);front=.375+.65*abs(center);te=float(geo.at(domain.piece(geo.D(center))[3],geo.D(center)))
            expected.append(dict(tau0=1.,tau1=0.,kbar=te-front,jbar=front-.375,jcbar=front-.375+.1)[k])
    error=float(np.max(abs(value-expected)));assert error<=1e-12
    return dict(analytic_geometry_error=error,all_criteria_pass=True)
