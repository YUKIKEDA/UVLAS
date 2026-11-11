"""Loads and objective strain-rate dissipation on the adopted shell/beam.

Design prototype only. Rotational forces and tangents use a fixed local
exponential chart, including the chart pullback of physical moments.
"""
import numpy as np
import beam_se3_v2 as beam
from pg_ritz_mechanics import MB

def cross_matrix(v):
    out=np.zeros(np.shape(v)[:-1]+(3,3),dtype=np.result_type(v,float));out[...,0,1]=-v[...,2];out[...,0,2]=v[...,1];out[...,1,0]=v[...,2];out[...,1,2]=-v[...,0];out[...,2,0]=-v[...,1];out[...,2,1]=v[...,0];return out

def surface(N,Da,Db,q,w,kind,value,area0=None,tangent=True):
    """Oriented face chart; positive pressure is compressive.

    Da/Db are derivatives wrt oriented surface coordinates, NOT Cartesian G.
    They include the outward orientation, e.g. by reversing Da for a minus face.
    """
    expected=() if kind=='pressure' else (3,)
    if np.shape(value)!=expected or not np.isfinite(value).all() or not np.isfinite(q).all() or np.any(w<=0):raise ValueError('invalid_surface_load_data')
    a=Da@q;b=Db@q;area=np.cross(a,b);norm=np.sqrt(np.sum(area*area,axis=1))
    if np.any(norm.real<=1e-12):raise ValueError('degenerate_loaded_surface')
    if kind=='pressure':density=-np.asarray(value)*area
    elif kind=='traction_current':density=norm[:,None]*np.asarray(value)
    elif kind=='traction_reference':
        if area0 is None:raise ValueError('missing_reference_area')
        density=np.sqrt(np.sum(area0*area0,axis=1))[:,None]*np.asarray(value)
    else:raise ValueError('unknown_surface_load_measure')
    f=N.T@(w[:,None]*density);K=None
    if tangent:
        Darea=-np.einsum('pj,pab->pajb',Da,cross_matrix(b))+np.einsum('pj,pab->pajb',Db,cross_matrix(a))
        if kind=='pressure':Dforce=-np.asarray(value)*Darea
        elif kind=='traction_current':Dnorm=np.einsum('pa,pajb->pjb',area/norm[:,None],Darea);Dforce=np.einsum('a,pjb->pajb',np.asarray(value),Dnorm)
        else:Dforce=np.zeros_like(Darea)
        K=np.einsum('p,pi,pajb->iajb',w,N,Dforce,optimize=True).reshape(q.size,q.size)
    return f,K,dict(area=area,density=density,point=N@q)

def point_wrench(r,R,t,c,force,moment,force_frame='world',moment_frame='world'):
    beam.validate(r,R,1.)
    if force_frame not in ['world','body'] or moment_frame not in ['world','body']:raise ValueError('invalid_wrench_frame')
    x,Q=beam.point_jet(r,R,t,c);X=beam.float_array(x);A=beam.float_array(Q);J=np.array([a.g for a in x]);H=np.array([a.H for a in x]);JQ=np.array([[a.g for a in row] for row in Q]);HQ=np.array([[a.H for a in row] for row in Q]);W=np.empty((3,12));DW=np.empty((3,12,12))
    for j in range(12):
        W[:,j]=np.array(beam.vee(A.T@JQ[:,:,j]),float)
        for k in range(12):DW[:,j,k]=np.array(beam.vee(JQ[:,:,k].T@JQ[:,:,j]+A.T@HQ[:,:,j,k]),float)
    fworld=A@force if force_frame=='body' else np.asarray(force);mbody=A.T@moment if moment_frame=='world' else np.asarray(moment);Df=np.einsum('abj,b->aj',JQ,force) if force_frame=='body' else np.zeros((3,12));Dm=np.einsum('baj,b->aj',JQ,moment) if moment_frame=='world' else np.zeros((3,12));f=J.T@fworld+W.T@mbody;K=np.einsum('a,aij->ij',fworld,H)+J.T@Df+np.einsum('a,aij->ij',mbody,DW)+W.T@Dm
    return f,K,dict(point=X,R=A,J=J,W=W,force_world=fworld,moment_world=A@mbody)

def viscosity_factor(D,scales=None):
    D=np.asarray(D,float)
    if D.shape!=(6,6) or not np.isfinite(D).all() or np.linalg.norm(D-D.T)>1e-12*max(np.linalg.norm(D),1e-300):raise ValueError('invalid_viscosity_matrix')
    scales=np.ones(6) if scales is None else np.asarray(scales,float)
    if scales.shape!=(6,) or not np.isfinite(scales).all() or np.any(scales<=0):raise ValueError('invalid_viscosity_scale')
    A=D*scales[:,None]*scales[None,:];vals,V=np.linalg.eigh(A)
    # Fail closed; do not repair a negative coefficient by clipping it.
    if vals[0]<0:raise ValueError('negative_viscosity')
    return np.sqrt(vals)[:,None]*V.T/scales[None,:]

def shell_damping(data,u,v,factor,tangent=True):
    if np.ndim(factor)!=2 or np.shape(factor)[1]!=6 or not np.isfinite(factor).all():raise ValueError('invalid_viscosity_factor')
    G=data['G'];w=data['weights'];F=np.eye(3)+np.einsum('ia,pib->pab',u,G);V=np.einsum('ia,pib->pab',v,G);Edot=(F.transpose(0,2,1)@V+V.transpose(0,2,1)@F)/2;e=np.einsum('kab,pab->pk',MB,Edot);D=factor.T@factor;s=e@D.T;S=np.einsum('pk,kab->pab',s,MB);P=F@S;force=np.einsum('p,pib,pab->ia',w,G,P,optimize=True);dissipation=float(np.real(np.sum(w[:,None]*(e@factor.T)**2)));Kq=Kv=None
    if tangent:
        ns=len(u);Kq=np.zeros((3*ns,3*ns),dtype=np.result_type(u,v,float));Kv=Kq.copy()
        for start in range(0,len(G),128):
            sl=slice(start,start+128);B=np.einsum('kij,pai,pnj->pkna',MB,F[sl],G[sl],optimize=True).reshape(len(G[sl]),6,3*ns);Bv=np.einsum('kij,pai,pnj->pkna',MB,V[sl],G[sl],optimize=True).reshape(len(G[sl]),6,3*ns);DB=np.einsum('kl,plj->pkj',D,B);DBv=np.einsum('kl,plj->pkj',D,Bv);Kv+=np.einsum('p,pki,pkj->ij',w[sl],B,DB,optimize=True);Kq+=np.einsum('p,pki,pkj->ij',w[sl],B,DBv,optimize=True)
        geo=np.einsum('p,pni,pij,pmj->nm',w,G,S,G,optimize=True);Kq+=np.kron(geo,np.eye(3))
    return force,Kq,Kv,dict(dissipation=dissipation,strain_rate=e,viscous_second_PK=S,F=F,V=V)

def beam_damping(r,R,L,v,factor,reference=None):
    beam.validate(r,R,L)
    if np.ndim(factor)!=2 or np.shape(factor)[1]!=6 or not np.isfinite(factor).all():raise ValueError('invalid_section_viscosity_factor')
    rs,Rs=beam.jets(r,R);strain=beam.relative(rs[0],Rs[0],rs[1],Rs[1])/L;B=np.array([a.g for a in strain]);H=np.array([a.H for a in strain]);D=factor.T@factor;rate=B@v;s=D@rate;f=L*B.T@s;Kv=L*B.T@D@B;Kq=L*(np.einsum('k,kij->ij',s,H)+B.T@D@np.einsum('kij,j->ki',H,v));return f,Kq,Kv,dict(strain_rate=rate,dissipation=float(L*np.dot(factor@rate,factor@rate)),B=B,H=H)

def strain_support(e,J,H,edot,K,D,edot_q):
    """Elastic/viscous support from the SAME relative-coordinate map.

    H[k,i,j]=e_k,ij. edot=J v+e_t is supplied from that map. The state
    full state derivative of edot is supplied as edot_q, including H*v and
    the moving-base term e_tq. It must come from the same relative map.
    Returns restoring (internal-side) force, not external applied force.
    """
    viscosity_factor(K);viscosity_factor(D)
    value=K@e+D@edot;force=J.T@value;Kv=J.T@D@J;Kq=np.einsum('k,kij->ij',value,H)+J.T@K@J+J.T@D@edot_q
    return force,Kq,Kv,dict(potential=float(e@K@e/2),dissipation=float(edot@D@edot))
