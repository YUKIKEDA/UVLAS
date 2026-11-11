"""Coefficient-compatible ANCF space, independent r / polynomial-normal-weighted d / c.

The numerator dictionary contains the normal and the exact shear-free director
companions of every position basis. It is compressed as a polynomial space,
before division by g**1.5. Thickness powers are separate coordinates.
"""
from types import SimpleNamespace
from functools import lru_cache
import numpy as np
from scipy.linalg import svd
from numpy.polynomial import legendre as L
import pg_ritz_geometry as geometry
from pg_ritz_shell_v3 import native
import pg_ritz_particular_v2 as particular


def rraw(p,n,clamped):
    jets=geometry.trial(p,n)
    if not clamped:return jets
    edge=p.copy();edge[:,0]=0;boundary=geometry.trial(edge,n)
    x=p[:,0,None]
    return [jets[0]-(1-x)*boundary[0],jets[1]+boundary[0],jets[2]-(1-x)*boundary[2],
            jets[3],jets[4]+boundary[2],jets[5]-(1-x)*boundary[5]]


def nodes(n,degree):
    t=-np.cos(np.arange(degree+1)*np.pi/degree)
    return np.r_[np.concatenate([(j+(t[:-1]+1)/2)/(2*n) for j in range(2*n)]),1.]


def hermite_nodes(n):
    inner=-np.cos(np.arange(1,8)*np.pi/8)
    coordinates=[];types=[]
    for j in range(2*n):
        coordinates.extend([j/(2*n)]*2);types.extend([0,1])
        coordinates.extend((j+(inner+1)/2)/(2*n));types.extend([0]*7)
    coordinates.extend([1.,1.]);types.extend([0,1])
    return np.array(coordinates),np.array(types)


@lru_cache(None)
def local_interpolator(degree,hermite):
    if hermite:
        inner=-np.cos(np.arange(1,8)*np.pi/8)
        t=np.r_[-1.,-1.,inner,1.,1.]
        V=L.legvander(t,10)
        D=L.legder(np.eye(11),axis=0)
        V[[1,-1]]=L.legvander(t[[1,-1]],9)@D
    else:
        t=-np.cos(np.arange(degree+1)*np.pi/degree)
        V=L.legvander(t,degree)
    inverse=np.linalg.inv(V)
    return inverse,L.legder(inverse,axis=0)


def interpolation(p,n,degree,coeff,hermite=False):
    """C0 tensor polynomial from shared Chebyshev-Lobatto nodal coefficients."""
    inverse,derivative=local_interpolator(degree,hermite)
    step=degree-1 if hermite else degree
    cell=np.minimum(np.maximum((p*2*n).astype(int),0),2*n-1)
    s=2*(p*2*n-cell)-1
    out=[np.empty((len(p),coeff.shape[-1])) for _ in range(3)]
    for i,j in np.unique(cell,axis=0):
        mask=(cell[:,0]==i)&(cell[:,1]==j)
        u,v=s[mask].T
        X=L.legvander(u,degree)@inverse;Y=L.legvander(v,degree)@inverse
        DX=L.legvander(u,degree-1)@derivative*(4*n)
        DY=L.legvander(v,degree-1)@derivative*(4*n)
        local=coeff[i*step:i*step+degree+1,j*step:j*step+degree+1].reshape((degree+1)**2,-1)
        for a,(U,V) in enumerate([(X,Y),(DX,Y),(X,DY)]):
            out[a][mask]=np.einsum('pi,pj->pij',U,V).reshape(sum(mask),-1)@local
    return out


def numerator_dictionary(p,n,parameters,T,clamped,Rcoeff=None):
    """Scalar components of H^{-1}[-n.v_x,-n.v_y,0] times g**1.5."""
    r,A,normal,dn,g=geometry.graph(p,parameters)
    fx,fy=A[:,2,0],A[:,2,1]
    R=interpolation(p,n,10,Rcoeff,True) if Rcoeff is not None else [a@T for a in rraw(p,n,clamped)[:3]]
    normal_raw=np.column_stack((-fx,-fy,np.ones(len(p))))
    a=np.column_stack((-(1+fy*fy),fx*fy,-fx))
    b=np.column_stack((fx*fy,-(1+fx*fx),-fy))
    companion=(a[:,:,None,None]*normal_raw[:,None,:,None]*R[1][:,None,None,:]
               +b[:,:,None,None]*normal_raw[:,None,:,None]*R[2][:,None,None,:]).reshape(len(p),-1)
    # Independent low-order numerator functions and exact reference normal.
    result=np.column_stack((native(p,n)[1][0],normal_raw*g[:,None],companion))
    if clamped:
        edge=p.copy();edge[:,0]=0
        # Evaluate the same dictionary at x=0; R itself already has the clamp.
        r,A,normal,dn,ge=geometry.graph(edge,parameters);ex,ey=A[:,2,0],A[:,2,1]
        E=interpolation(edge,n,10,Rcoeff,True) if Rcoeff is not None else [a@T for a in rraw(edge,n,clamped)[:3]]
        en=np.column_stack((-ex,-ey,np.ones(len(p))))
        ea=np.column_stack((-(1+ey*ey),ex*ey,-ex));eb=np.column_stack((ex*ey,-(1+ex*ex),-ey))
        ec=(ea[:,:,None,None]*en[:,None,:,None]*E[1][:,None,None,:]
            +eb[:,:,None,None]*en[:,None,:,None]*E[2][:,None,None,:]).reshape(len(p),-1)
        boundary=np.column_stack((native(edge,n)[1][0],en*ge[:,None],ec))
        result-=(1-p[:,0,None])*boundary
    return result


def independent(A,weight):
    norm=np.linalg.norm(A*weight[:,None],axis=0)
    keep=norm>1e-11*max(norm)
    scaled=A[:,keep]*weight[:,None]/norm[keep]
    U,s,V=svd(scaled,full_matrices=False,check_finite=False)
    ranks=[int(sum(s>tol*max(1.,np.linalg.norm(scaled)))) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('coefficient_dictionary_rank:'+str(ranks))
    rank=ranks[0];T=np.zeros((A.shape[1],rank));T[keep]=V[:rank].T/norm[keep,None]/s[:rank]
    return T,ranks,float(s[rank-1]/s[0])


def prepare(n,parameters,order=14,clamped=True):
    p,w=geometry.quadrature(n,order)
    T,rr,rc=independent(rraw(p,n,clamped)[0],np.sqrt(w))
    hx,ht=hermite_nodes(n);X,Y=np.meshgrid(hx,hx,indexing='ij')
    dx,dy=np.meshgrid(ht,ht,indexing='ij');jet_index=np.where(dx.ravel(),np.where(dy.ravel(),4,1),np.where(dy.ravel(),2,0))
    hp=np.column_stack((X.ravel(),Y.ravel()));hjets=rraw(hp,n,clamped)
    controls=np.zeros_like(hjets[0])
    for j in [0,1,2,4]:controls[jet_index==j]=hjets[j][jet_index==j]
    controls*=((1/(4*n))**(dx+dy)).ravel()[:,None]
    Rcoeff=(controls@T).reshape(len(hx),len(hx),-1)
    raw=numerator_dictionary(p,n,parameters,T,clamped,Rcoeff)
    D,dr,dc=independent(raw,np.sqrt(w))
    degree=13;grid=nodes(n,degree);X,Y=np.meshgrid(grid,grid,indexing='ij')
    points=np.column_stack((X.ravel(),Y.ravel()))
    coefficients=(numerator_dictionary(points,n,parameters,T,clamped,Rcoeff)@D).reshape(len(grid),len(grid),-1)
    native_c=native(p,n)[2][0];keep_c=np.arange(native_c.shape[1])
    if clamped:keep_c=keep_c[n+1:]
    C=native_c[:,keep_c];CT,cr,cc=independent(C,np.sqrt(w))
    model=dict(n=n,parameters=parameters,clamped=clamped,Tr=T,Rcoeff=Rcoeff,Dcoeff=coefficients,
        degree=degree,Cindices=keep_c,Tc=CT,nr=T.shape[1],nd=D.shape[1],nc=CT.shape[1],
        lift=dict(n=n,parameters=parameters),ranks=dict(r=rr,d=dr,c=cr),
        rcond=dict(r=rc,d=dc,c=cc))
    model['ns']=model['nr']+model['nd']+model['nc']
    # Independently evaluate the polynomial representation on the preparation nodes.
    interpolated=interpolation(p,n,degree,coefficients)[0]
    model['numerator_interpolation']=float(np.linalg.norm(interpolated-raw@D)/np.linalg.norm(raw@D))
    if model['numerator_interpolation']>1e-9:raise ValueError('numerator_polynomial_representation')
    return model


def fields(model,p):
    n=model['n'];R=interpolation(p,n,10,model['Rcoeff'],True)
    P=interpolation(p,n,model['degree'],model['Dcoeff'])
    _,A,_,_,g=geometry.graph(p,model['parameters']);kx,ky,tau=model['parameters']
    fx,fy=A[:,2,0],A[:,2,1];dg=[2*(kx*fx+tau*fy),2*(tau*fx+ky*fy)]
    D=[P[0]/g[:,None]**1.5]+[P[j+1]/g[:,None]**1.5-1.5*P[0]*dg[j][:,None]/g[:,None]**2.5 for j in range(2)]
    C=[a[:,model['Cindices']]@model['Tc'] for a in native(p,n)[2]]
    out=[]
    for j in range(3):
        a=np.zeros((3,len(p),model['ns']))
        a[0,:,:model['nr']]=R[j]
        a[1,:,model['nr']:model['nr']+model['nd']]=D[j]
        a[2,:,model['nr']+model['nd']:]=C[j]
        out.append(a)
    return out


def compatible_particular(model,p):
    """Reference-only polynomial fields for the stable coarse calculation."""
    v=particular.jets(p,model['n'],model['parameters'])
    _,A,normal,dn,g=geometry.graph(p,model['parameters'])
    H=np.concatenate((A.transpose(0,2,1),normal[:,None,:]),axis=1)
    zeros=np.zeros((len(p),32*model['n']))
    rhs=np.stack((-np.einsum('pa,pan->pn',normal,v[1]),-np.einsum('pa,pan->pn',normal,v[2]),zeros),axis=1)
    d=np.linalg.solve(H,rhs);ds=[];kx,ky,tau=model['parameters']
    for axis,(i,j) in enumerate([(3,4),(4,5)]):
        db=np.stack((-np.einsum('pa,pan->pn',dn[:,:,axis],v[1])-np.einsum('pa,pan->pn',normal,v[i]),-np.einsum('pa,pan->pn',dn[:,:,axis],v[2])-np.einsum('pa,pan->pn',normal,v[j]),zeros),axis=1)
        dH=np.zeros_like(H);dH[:,0,2]=kx if axis==0 else tau;dH[:,1,2]=tau if axis==0 else ky;dH[:,2]=dn[:,:,axis]
        ds.append(np.linalg.solve(H,db-dH@d))
    return [np.concatenate((a,a),axis=2) for a in v],np.concatenate((d,d),axis=2),[np.concatenate((a,a),axis=2) for a in ds]


lift=SimpleNamespace(geometry=geometry,evaluate=compatible_particular)


def kinematics(model,p3,h,variable=False):
    p,back=np.unique(p3[:,:2],axis=0,return_inverse=True);fs=fields(model,p);X,J,det,H,Ha=lift.geometry.reference(p3[:,:2],p3[:,2],h,model['parameters'],variable);z=H*p3[:,2]/2;za=Ha*p3[:,2,None]/2;parts=[v[:,back] for v in fs];N=parts[0][0]+z[:,None]*parts[0][1]+z[:,None]**2/2*parts[0][2];G=np.empty((len(p3),model['ns'],3))
    for a in range(2):G[:,:,a]=parts[a+1][0]+z[:,None]*parts[a+1][1]+z[:,None]**2/2*parts[a+1][2]+za[:,a,None]*(parts[0][1]+z[:,None]*parts[0][2])
    G[:,:,2]=H[:,None]/2*(parts[0][1]+z[:,None]*parts[0][2]);G=G@np.linalg.inv(J);MB=np.zeros((6,3,3));MB[0,0,0]=MB[1,1,1]=MB[2,2,2]=1
    for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2)
    B=np.einsum('kij,pnj->pkni',MB,G).reshape(len(p3),6,3*model['ns']);return N,G,B,det,X
