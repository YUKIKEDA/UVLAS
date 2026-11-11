"""Analytic SAME regularized W-P1 boundary plus constant sheet integral.

Inputs/outputs are local design probes, not public solver APIs. Jacobians are
with respect to point[3], vertices[9], core[1], with coefficients held fixed.
"""
import numpy as np
import uvlm_dual13 as ad
from check_uvlm_adaptive_quadrature import triangle_distance

def segment(p,a,b,eps,ma,mb):
    edge=b-a;L=ad.norm(edge);t=ad.scale(edge,1/L);r=p-a;s=ad.dot(r,t);z=s-L;perp=r-ad.scale(t,s);k2=ad.dot(perp,perp)+eps*eps;ra=ad.sqrt(s*s+k2);rb=ad.sqrt(z*z+k2)
    if ad.val(s*z)>0:I0=L*(s+z)/(ra*rb*(s*rb+z*ra))
    else:I0=(s/ra-z/rb)/k2
    # Exact positive weighted centroid: I1/I0=L*ra/(ra+rb).
    # This avoids cancellation in s*I0+1/ra-1/rb.
    strength=(ma*rb+mb*ra)/(ra+rb)
    return ad.scale(ad.cross(t,perp),I0*strength)

def logarithmic_edge(p,a,b,eps):
    edge=b-a;L=ad.norm(edge);t=ad.scale(edge,1/L);r=p-a;s=ad.dot(r,t);z=L-s;perp=r-ad.scale(t,s);k2=ad.dot(perp,perp)+eps*eps;ra=ad.sqrt(s*s+k2);rb=ad.sqrt(z*z+k2);aa=ad.absolute(s);bb=ad.absolute(z)
    gap=aa+bb-L+k2/(ra+aa)+k2/(rb+bb)
    if ad.val(gap)<=0:raise ValueError('unresolved_edge_logarithm')
    return ad.log1p(2*L/gap)

def sheet_vector(p,x,eps):
    e=x[1]-x[0];f=x[2]-x[0];cr=ad.cross(e,f);jac=ad.norm(cr);n=ad.scale(cr,1/jac);d=ad.dot(p-x[0],n);out=np.array([ad.Dual(0.) for _ in range(3)],object)
    for i,j in [(0,1),(1,2),(2,0)]:
        t=ad.scale(x[j]-x[i],1/ad.norm(x[j]-x[i]));out+=ad.scale(ad.cross(t,n),logarithmic_edge(p,x[i],x[j],eps))
    h2=d*d+eps*eps
    if ad.val(h2)==0:
        # Coplanar OUTSIDE the triangle. The h->0 solid-angle quotient is
        # finite; preserving d's gradient is essential for ambient derivatives.
        r=[p-a for a in x];rr=[ad.norm(a) for a in r];den=rr[0]*rr[1]*rr[2]+ad.dot(r[0],r[1])*rr[2]+ad.dot(r[1],r[2])*rr[0]+ad.dot(r[2],r[0])*rr[1];normal=d*(2*jac/den)
    else:
        h=ad.sqrt(h2);projected=p-ad.scale(n,d);virtual=projected+ad.scale(n,h);r=[virtual-a for a in x];rr=[ad.norm(a) for a in r];num=h*jac;den=rr[0]*rr[1]*rr[2]+ad.dot(r[0],r[1])*rr[2]+ad.dot(r[1],r[2])*rr[0]+ad.dot(r[2],r[0])*rr[1];omega=2*ad.atan2(num,den);normal=d/h*omega
    return out+ad.scale(n,normal),n

def far_area(p,x,eps,order):
    """Independent smooth far branch avoids boundary/sheet cancellation."""
    q,w=np.polynomial.legendre.leggauss(order);q=(q+1)/2;w=w/2;cr=ad.cross(x[1]-x[0],x[2]-x[0]);jac=ad.norm(cr);n=ad.scale(cr,1/jac);out=[[ad.Dual(0.) for _ in range(3)] for k in range(3)]
    for i,u in enumerate(q):
        for j,v in enumerate(q):
            bary=np.array([(1-u)*(1-v),u,(1-u)*v]);y=np.array([sum(bary[k]*x[k,a] for k in range(3)) for a in range(3)],object);r=p-y;r2=ad.dot(r,r);den=r2+eps*eps;rn=ad.dot(r,n);kernel=ad.scale(ad.scale(n,r2-2*eps*eps)-ad.scale(r,3*rn),jac*w[i]*w[j]*(1-u)/(den*den*ad.sqrt(den)))
            for k in range(3):
                for a in range(3):out[k][a]+=bary[k]*kernel[a]
    return out

def triangle(point,vertices,eps=0.,incident=False):
    p=np.asarray(point,float);x=np.asarray(vertices,float)
    if p.shape!=(3,) or x.shape!=(3,3) or not np.isfinite(p).all() or not np.isfinite(x).all() or not np.isfinite(eps) or eps<0:raise ValueError('invalid_induction_data')
    lengths=[np.linalg.norm(x[j]-x[i]) for i,j in [(0,1),(1,2),(2,0)]];L=max(lengths);jac=np.linalg.norm(np.cross(x[1]-x[0],x[2]-x[0]))
    if not np.isfinite(L) or not np.isfinite(jac) or L<=0 or jac<=1e-12*L*L:raise ValueError('degenerate_induction_triangle')
    if eps>0 and (eps/L)**2==0:raise ValueError('core_precision_underflow')
    distance=triangle_distance((p-x[0])/L,(x-x[0])/L);contact=distance<=1e-11
    if incident and not contact:raise ValueError('incidence_geometry_mismatch')
    if contact and not incident:raise ValueError('nonincident_surface_contact')
    if contact and eps==0:raise ValueError('ideal_surface_self_not_supported')
    origin=x[0].copy();pp=np.array([ad.Dual.variable((p[j]-origin[j])/L,j,1/L) for j in range(3)],object);xx=np.array([[ad.Dual.variable((x[i,j]-origin[j])/L,3+3*i+j,1/L) for j in range(3)] for i in range(3)],object);core=ad.Dual.variable(eps/L,12,1/L)
    if distance>4:
        low=far_area(pp,xx,core,8);high=far_area(pp,xx,core,12);v=np.array([[a.v for a in row] for row in high])/(4*np.pi*L);j=np.array([[a.g for a in row] for row in high])/(4*np.pi*L);lv=np.array([[a.v for a in row] for row in low])/(4*np.pi*L);lj=np.array([[a.g for a in row] for row in low])/(4*np.pi*L)
        if np.linalg.norm(v-lv)>1e-10+1e-8*np.linalg.norm(v) or np.linalg.norm(j-lj)*L>1e-10+1e-8*np.linalg.norm(j)*L:raise ValueError('far_quadrature_unresolved')
        return v,j,dict(branch='far_Gauss12_checked_against8',normalized_distance=distance,incident=False,eps=float(eps))
    I,n=sheet_vector(pp,xx,core);e=xx[1]-xx[0];f=xx[2]-xx[0];ee=ad.dot(e,e);ef=ad.dot(e,f);ff=ad.dot(f,f);det=ee*ff-ef*ef;values=[];derivatives=[]
    for k in range(3):
        mu=np.eye(3)[k];b1=mu[1]-mu[0];b2=mu[2]-mu[0];grad=ad.scale(e,(ff*b1-ef*b2)/det)+ad.scale(f,(ee*b2-ef*b1)/det);v=ad.cross(ad.cross(n,grad),I)
        for i,j in [(0,2),(2,1),(1,0)]:v+=segment(pp,xx[i],xx[j],core,mu[i],mu[j])
        v=ad.scale(v,1/(4*np.pi*L));values.append([a.v for a in v]);derivatives.append([a.g for a in v])
    values=np.array(values);derivatives=np.array(derivatives)
    if not np.isfinite(values).all() or not np.isfinite(derivatives).all():raise ValueError('induction_precision_unresolved')
    return values,derivatives,dict(branch='analytic_boundary_and_sheet',normalized_distance=distance,incident=bool(incident),eps=float(eps))

def patch(point,quad,eps=0.,incident=False):
    quad=np.asarray(quad);result=np.zeros((2,3));J=np.zeros((2,3,16));maps=np.array([[1.,0],[1.,0],[0,1.],[0,1.]])
    contacts=[]
    for ids in [[0,2,1],[0,3,2]]:contacts.append(triangle_distance(point,quad[ids])<=1e-11*max(np.linalg.norm(quad[i]-quad[j]) for i in range(4) for j in range(i)))
    if incident and not any(contacts):raise ValueError('incidence_geometry_mismatch')
    for ids,local_contact in zip([[0,2,1],[0,3,2]],contacts):
        v,j,_=triangle(point,quad[ids],eps,incident and local_contact);result+=maps[ids].T@v;D=np.einsum('ik,iaj->kaj',maps[ids],j);J[:,:,:3]+=D[:,:,:3];J[:,:,15]+=D[:,:,12]
        for i,k in enumerate(ids):J[:,:,3+3*k:6+3*k]+=D[:,:,3+3*i:6+3*i]
    return result,J
