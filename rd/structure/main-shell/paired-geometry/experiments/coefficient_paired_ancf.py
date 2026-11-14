"""PG pair space represented and decomposed in polynomial coefficient space.

Native d is numerator-polynomial / g**1.5, so compression never tries to
distinguish a near-polynomial normal from a Cartesian polynomial d. The r
image and its pure-d kernel are separated BEFORE thickness scaling.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import svd
import coefficient_compatible_ancf as polynomial
import pg_ritz_lift_v8 as reference_lift
geometry=polynomial.geometry
lift=polynomial.lift
compatible_particular=polynomial.compatible_particular


def r_raw(p,n,ref):
    v,_,_=reference_lift.evaluate(ref,p)
    count=(n+2)**2
    native=[a[:,:count] for a in geometry.trial(p,n)]
    return [np.column_stack((a,b.reshape(len(p),-1))) for a,b in zip(native,v)]


def r_clamped(p,n,ref):
    a=r_raw(p,n,ref);edge=p.copy();edge[:,0]=0;b=r_raw(edge,n,ref);x=p[:,0,None]
    return [a[0]-(1-x)*b[0],a[1]+b[0],a[2]-(1-x)*b[2],a[3],a[4]+b[2],a[5]-(1-x)*b[5]]


def numerator_raw(p,n,parameters,ref):
    def raw(p):
        _,A,_,_,g=geometry.graph(p,parameters);fx,fy=A[:,2,0],A[:,2,1]
        v,_,_=reference_lift.evaluate(ref,p)
        t1=-fx[:,None]*v[1][:,0]-fy[:,None]*v[1][:,1]+v[1][:,2]
        t2=-fx[:,None]*v[2][:,0]-fy[:,None]*v[2][:,1]+v[2][:,2]
        num=np.stack((-(1+fy[:,None]**2)*t1+fx[:,None]*fy[:,None]*t2,
                      fx[:,None]*fy[:,None]*t1-(1+fx[:,None]**2)*t2,
                      -fx[:,None]*t1-fy[:,None]*t2),axis=1).reshape(len(p),-1)
        basic=np.column_stack((polynomial.native(p,n)[1][0],np.column_stack((-fx,-fy,np.ones(len(p))))*g[:,None]))
        return basic,num
    basic,num=raw(p);edge=p.copy();edge[:,0]=0;eb,en=raw(edge)
    return basic-(1-p[:,0,None])*eb,num-(1-p[:,0,None])*en


def rank_of(s,A,label):
    ranks=[int(sum(s>tol*max(1.,np.linalg.norm(A)))) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError(label+':'+str(ranks))
    return ranks


def prepare(n,parameters,order=15):
    ref=reference_lift.prepare(n,parameters,13)
    # Canonical C1 position coefficients: endpoint values and derivatives are shared.
    hx,ht=polynomial.hermite_nodes(n);X,Y=np.meshgrid(hx,hx,indexing='ij');dx,dy=np.meshgrid(ht,ht,indexing='ij')
    hp=np.column_stack((X.ravel(),Y.ravel()));jets=r_clamped(hp,n,ref)
    ji=np.where(dx.ravel(),np.where(dy.ravel(),4,1),np.where(dy.ravel(),2,0))
    control=np.zeros_like(jets[0])
    for j in [0,1,2,4]:control[ji==j]=jets[j][ji==j]
    control*=((1/(4*n))**(dx+dy)).ravel()[:,None]
    nr0=(n+2)**2;nd0=(2*n+1)**2+3;npair=control.shape[1]-nr0
    total=nr0+nd0+npair
    Rc=np.zeros((len(hx),len(hx),total));Rc[:,:,:nr0]=control[:,:nr0].reshape(len(hx),len(hx),nr0)
    Rc[:,:,nr0+nd0:]=control[:,nr0:].reshape(len(hx),len(hx),npair)
    # C0 director numerator coefficients, shared at every interface node.
    degree=13;gx=polynomial.nodes(n,degree);X,Y=np.meshgrid(gx,gx,indexing='ij');dp=np.column_stack((X.ravel(),Y.ravel()))
    basic,pair=numerator_raw(dp,n,parameters,ref)
    Dc=np.zeros((len(gx),len(gx),total));Dc[:,:,nr0:nr0+nd0]=basic.reshape(len(gx),len(gx),nd0)
    Dc[:,:,nr0+nd0:]=pair.reshape(len(gx),len(gx),npair)
    p,w=geometry.quadrature(n,order);sqrt=np.sqrt(w)[:,None]
    R=polynomial.interpolation(p,n,10,Rc,True)[0]
    D=polynomial.interpolation(p,n,degree,Dc)[0]
    joint=np.vstack((R*sqrt,D*sqrt));norm=np.linalg.norm(joint,axis=0)
    keep=norm>1e-11*max(norm);A=joint[:,keep]/norm[keep]
    U,s,V=svd(A,full_matrices=False,check_finite=False);ranks=rank_of(s,A,'paired_coefficient_rank');k=ranks[0]
    T=np.zeros((total,k));T[keep]=V[:k].T/norm[keep,None]/s[:k]
    Rc=Rc@T;Dc=Dc@T;R=R@T;D=D@T
    # Hierarchical coordinates: r image, d-only kernel. No thickness enters here.
    U,s,V=svd(R*sqrt,full_matrices=False,check_finite=False);rr=rank_of(s,R*sqrt,'position_image_rank');nr=rr[0]
    Vr=V[:nr].T/s[:nr];Vd=V[nr:].T
    RI=Rc@Vr;DI=Dc@Vr;DK=Dc@Vd;DKvalues=D@Vd
    U,sd,V=svd(DKvalues*sqrt,full_matrices=False,check_finite=False)
    dr=rank_of(sd,DKvalues*sqrt,'director_kernel_rank');nd=dr[0]
    if nd!=k-nr:raise ValueError('director_kernel_dimension')
    W=V[:nd].T/sd[:nd];DK=DK@W;DKvalues=DKvalues@W
    cross=DKvalues.T@(w[:,None]*(D@Vr));DI-=DK@cross
    # RI has full r rank; the second block has exactly zero r by definition.
    C=polynomial.native(p,n)[2][0];indices=np.arange(C.shape[1])[n+1:]
    CT,cr,cc=polynomial.independent(C[:,indices],np.sqrt(w))
    discarded=float(np.linalg.norm((R*sqrt)@Vd)/max(1.,np.linalg.norm(R*sqrt)))
    if discarded>1e-9:raise ValueError('nonzero_position_kernel:'+str(discarded))
    return dict(n=n,parameters=parameters,Rcoeff=RI,Dcoupled=DI,Dcoeff=DK,degree=degree,
        Cindices=indices,Tc=CT,nr=nr,nd=nd,nc=CT.shape[1],ns=nr+nd+CT.shape[1],
        lift=dict(n=n,parameters=parameters),ranks=dict(joint=ranks,r=rr,d=dr,c=cr),
        position_kernel_defect=discarded,reference_ritz=ref)


def fields(model,p):
    n=model['n'];R=polynomial.interpolation(p,n,10,model['Rcoeff'],True)
    P=polynomial.interpolation(p,n,model['degree'],np.concatenate((model['Dcoupled'],model['Dcoeff']),axis=2))
    _,A,_,_,g=geometry.graph(p,model['parameters']);kx,ky,tau=model['parameters'];fx,fy=A[:,2,0],A[:,2,1]
    dg=[2*(kx*fx+tau*fy),2*(tau*fx+ky*fy)]
    D=[P[0]/g[:,None]**1.5]+[P[j+1]/g[:,None]**1.5-1.5*P[0]*dg[j][:,None]/g[:,None]**2.5 for j in range(2)]
    C=[a[:,model['Cindices']]@model['Tc'] for a in polynomial.native(p,n)[2]]
    out=[]
    for j in range(3):
        a=np.zeros((3,len(p),model['ns']));a[0,:,:model['nr']]=R[j]
        a[1,:,:model['nr']+model['nd']]=D[j];a[2,:,model['nr']+model['nd']:]=C[j];out.append(a)
    return out


def kinematics(model,p3,h,variable=False):
    p,back=np.unique(p3[:,:2],axis=0,return_inverse=True);fs=fields(model,p);X,J,det,H,Ha=lift.geometry.reference(p3[:,:2],p3[:,2],h,model['parameters'],variable);z=H*p3[:,2]/2;za=Ha*p3[:,2,None]/2;parts=[v[:,back] for v in fs];N=parts[0][0]+z[:,None]*parts[0][1]+z[:,None]**2/2*parts[0][2];G=np.empty((len(p3),model['ns'],3))
    for a in range(2):G[:,:,a]=parts[a+1][0]+z[:,None]*parts[a+1][1]+z[:,None]**2/2*parts[a+1][2]+za[:,a,None]*(parts[0][1]+z[:,None]*parts[0][2])
    G[:,:,2]=H[:,None]/2*(parts[0][1]+z[:,None]*parts[0][2]);G=G@np.linalg.inv(J);MB=np.zeros((6,3,3));MB[0,0,0]=MB[1,1,1]=MB[2,2,2]=1
    for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2)
    B=np.einsum('kij,pnj->pkni',MB,G).reshape(len(p3),6,3*model['ns']);return N,G,B,det,X
