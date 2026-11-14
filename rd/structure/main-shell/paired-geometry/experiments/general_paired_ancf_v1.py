"""Coefficient-direct PG spaces on nonquadratic polynomial graph patches.

The exact integral lift replaces the rejected membrane-Ritz binding. Full
and homogeneous spaces use the same explicit trace lifting: V=S+L(TS),
H=(I-LT)S. In particular H is exactly ker(T) in V over Q.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from functools import lru_cache
import hashlib,json,time,os
import numpy as np
from scipy.linalg import qr,svd
from numpy.polynomial import legendre as L
from flint import fmpq_mat,arb_mat,ctx
import general_paired_polynomials_v1 as exact
import polynomial_reference_surface as geometry
import paired_exact_polynomials as pexact
import coefficient_compatible_ancf as base
from exact_paired_ancf import rank_check,evaluate_legendre
ROOT=rd_paths.ROOT

def canonical(n,parameters,boundary='clamped',arithmetic='arb192'):
    parameters=tuple(float(v) for v in parameters)
    version=hashlib.sha256(Path(exact.__file__).read_bytes()).hexdigest()[:12]
    key=hashlib.sha256(repr((n,tuple(parameters),version,boundary,arithmetic)).encode()).hexdigest()[:16]
    directory=ROOT/'build/analysis/general-paired-space';directory.mkdir(exist_ok=True)
    path=directory/(key+'.npz');meta=directory/(key+'.json')
    if path.exists() and meta.exists():
        record=json.loads(meta.read_text())
        if hashlib.sha256(path.read_bytes()).hexdigest()==record['artifact_sha256']:
            a=np.load(path);return a['coefficients'],dict(record,runtime_cache_hit=True)
    tick=time.monotonic();source=exact.build(n,parameters,boundary);degree=source.shape[3]-1;A=source.reshape(-1,source.shape[-1])
    norm=np.linalg.norm(A,axis=1);rows=np.flatnonzero(norm);E=A[rows]/norm[rows,None]
    ss=svd(E,compute_uv=False,check_finite=False)
    ranks=[int(sum(ss>tol*max(1.,np.linalg.norm(E)))) for tol in [1e-9,1e-10,1e-11]]
    rsize=source[0].size//source.shape[-1];rrows=np.flatnonzero(norm[:rsize])
    ER=A[rrows]/norm[rrows,None];ss=svd(ER,compute_uv=False,check_finite=False)
    rr=[int(sum(ss>tol*max(1.,np.linalg.norm(ER)))) for tol in [1e-9,1e-10,1e-11]]
    rational=exact.rational_matrix(n,parameters,rows,np.arange(source.shape[-1]),boundary)
    exact_rank=int(rational.rank())
    position=exact.rational_matrix(n,parameters,rrows,np.arange(source.shape[-1]),boundary)
    exact_r_rank=int(position.rank())
    k=exact_rank;nr=exact_r_rank
    _,_,column_order=qr(E,mode='economic',pivoting=True,check_finite=False);columns=column_order[:k]
    # Force position pivots first so the remaining coordinates have exactly zero r.
    _,_,order=qr(A[rrows][:,columns].T,mode='economic',pivoting=True,check_finite=False)
    r_pivots=rrows[order[:nr]]
    Qp,_=qr((A[r_pivots][:,columns]/norm[r_pivots,None]).T,mode='economic',check_finite=False)
    drows=np.flatnonzero(norm[rsize:])+rsize;D=A[drows][:,columns]
    residual=D-(D@Qp)@Qp.T
    _,_,order=qr(residual.T,mode='economic',pivoting=True,check_finite=False)
    pivots=np.r_[r_pivots,drows[order[:k-nr]]]
    print('EXACT_RANK',n,k,nr,'seconds',time.monotonic()-tick,flush=True)
    lookup={int(v):i for i,v in enumerate(rows)}
    pivot=fmpq_mat([[rational[lookup[int(i)],int(j)] for j in columns] for i in pivots])
    selected=fmpq_mat([[rational[i,int(j)] for j in columns] for i in range(len(rows))])
    attempts=[]
    for bits in ([192,256,320,384,512] if arithmetic=='arb192' else [0]):
        attempts.append(bits)
        if bits:
            ctx.prec=bits;ctx.threads=4
            try:basis=arb_mat(selected)*arb_mat(pivot).inv()
            except (ValueError,ZeroDivisionError):continue
        else:basis=selected*pivot.inv()
        print('CERTIFIED_BASIS',n,bits,time.monotonic()-tick,flush=True)
        output=np.zeros((A.shape[0],k));zero_r=True;radius=0.;finite=True
        for i,row in enumerate(rows):
            for j in range(k):
                value=basis[i,j]
                if bits:
                    if not value.is_finite():finite=False;continue
                    radius=max(radius,float(value.rad()))
                    if row<rsize and j>=nr:
                        if not value.contains(0):zero_r=False
                        output[int(row),j]=0.
                    else:output[int(row),j]=float(value)
                else:
                    if row<rsize and j>=nr and value!=0:zero_r=False
                    output[int(row),j]=float(value)
        if finite and zero_r and radius<=1e-25:break
    else:raise ValueError('uncertified_coefficient_basis')
    if not zero_r:raise ValueError('nonzero_exact_position_kernel')
    if radius>1e-25:raise ValueError('coefficient_interval_radius:'+str(radius))
    C=output.reshape(2,2*n,2*n,degree+1,degree+1,k)
    record=dict(n=n,parameters=parameters,boundary=boundary,r_degree=exact.specification(parameters)[1],d_degree=degree,exact_rank=k,exact_r_rank=nr,exact_d_kernel=k-nr,
        numeric_guide_ranks=ranks,numeric_position_guide=rr,pure_d_r_exactly_zero=True,
        arithmetic=arithmetic,max_coefficient_radius=radius,
        precision_bits=bits,attempted_precisions=attempts,
        seconds=time.monotonic()-tick,coefficient_source_sha256=hashlib.sha256(Path(exact.__file__).read_bytes()).hexdigest())
    temporary=path.with_name(path.stem+'.'+str(os.getpid())+'.npz')
    np.savez_compressed(temporary,coefficients=C);record['artifact']=path.relative_to(ROOT).as_posix()
    record['artifact_sha256']=hashlib.sha256(temporary.read_bytes()).hexdigest()
    os.replace(temporary,path);temporary_meta=meta.with_name(meta.stem+'.'+str(os.getpid())+'.json')
    temporary_meta.write_text(json.dumps(record,indent=2),encoding='utf-8');os.replace(temporary_meta,meta)
    print('EXACT_CANONICAL',record,flush=True)
    return C,dict(record,runtime_cache_hit=False)


def hermite_nodes(n,degree):
    inner=-np.cos(np.arange(1,degree-2)*np.pi/(degree-2));coords=[];types=[]
    for j in range(2*n):
        coords.extend([j/(2*n)]*2);types.extend([0,1])
        coords.extend((j+(inner+1)/2)/(2*n));types.extend([0]*(degree-3))
    coords.extend([1.,1.]);types.extend([0,1])
    return np.array(coords),np.array(types)


@lru_cache(None)
def local_interpolator(degree,hermite):
    if hermite:
        inner=-np.cos(np.arange(1,degree-2)*np.pi/(degree-2))
        t=np.r_[-1.,-1.,inner,1.,1.];V=L.legvander(t,degree)
        V[[1,-1]]=L.legvander(t[[1,-1]],degree-1)@L.legder(np.eye(degree+1),axis=0)
    else:V=L.legvander(-np.cos(np.arange(degree+1)*np.pi/degree),degree)
    inverse=np.linalg.inv(V)
    return inverse,L.legder(inverse,axis=0)


def interpolation(points,n,degree,coeff,hermite=False):
    inverse,derivative=local_interpolator(degree,hermite);step=degree-1 if hermite else degree
    cell=np.minimum(np.maximum((points.real*2*n).astype(int),0),2*n-1);t=2*(points*2*n-cell)-1
    out=[np.empty((len(points),coeff.shape[-1]),dtype=np.result_type(points,coeff)) for _ in range(3)]
    for i,j in np.unique(cell,axis=0):
        mask=(cell[:,0]==i)&(cell[:,1]==j);x,y=t[mask].T
        X=L.legvander(x,degree)@inverse;Y=L.legvander(y,degree)@inverse
        DX=L.legvander(x,degree-1)@derivative*(4*n);DY=L.legvander(y,degree-1)@derivative*(4*n)
        local=coeff[i*step:i*step+degree+1,j*step:j*step+degree+1].reshape((degree+1)**2,-1)
        for a,(U,V) in enumerate([(X,Y),(DX,Y),(X,DY)]):out[a][mask]=np.einsum('pi,pj->pij',U,V).reshape(sum(mask),-1)@local
    return out


def reference_normal_coefficients(n,parameters,degree):
    coeffs,rd,dd=exact.specification(parameters);out=np.zeros((2*n,2*n,degree+1,degree+1,3))
    width=geometry.chart(parameters)[1]
    for cx in range(2*n):
        for cy in range(2*n):
            x,y,f=exact.graph_polynomial(n,cx,cy,coeffs);fx=pexact.scale(pexact.derivative(f,0,n),1/width);fy=pexact.derivative(f,1,n)
            g=pexact.add({(0,0):pexact.F(1)},pexact.mul(fx,fx),pexact.mul(fy,fy))
            for j,a in enumerate([pexact.scale(pexact.mul(g,fx),-1),pexact.scale(pexact.mul(g,fy),-1),g]):
                for (i,k),v in pexact.to_legendre(a,degree).items():out[cx,cy,i,k,j]=float(v)
    return out


def prepare(n,parameters,order=None,boundary='clamped'):
    parameters=tuple(float(v) for v in parameters);C,record=canonical(n,parameters,boundary)
    nr=record['exact_r_rank'];nd=record['exact_d_kernel'];rd=record['r_degree'];dd=record['d_degree']
    if order is None:order=dd+2
    def weights(degree):
        a,b=np.meshgrid(np.arange(degree+1),np.arange(degree+1),indexing='ij')
        return np.broadcast_to(1/(2*n)/np.sqrt((2*a+1)*(2*b+1)),(2*n,2*n,degree+1,degree+1)).reshape(-1,1)
    WR,WD=weights(rd),weights(dd)
    R=C[0,:,:,:rd+1,:rd+1,:nr];DI=C[1,:,:,:,:,:nr];DK=C[1,:,:,:,:,nr:]
    RA=R.reshape(-1,nr)*WR;DA=DK.reshape(-1,nd)*WD
    U,s,V=svd(RA,full_matrices=False,check_finite=False);rr=rank_check(RA,'generic_scaled_r_rank')
    if rr!=[nr]*3:raise ValueError('r_image_dimension')
    R=R@(V.T/s);DI=DI@(V.T/s)
    U,s,V=svd(DA,full_matrices=False,check_finite=False);dr=rank_check(DA,'generic_scaled_d_rank')
    if dr!=[nd]*3:raise ValueError('d_kernel_dimension')
    DK=DK@(V.T/s);DKvalues=DK.reshape(-1,nd)*WD
    DI-=DK@(DKvalues.T@(DI.reshape(-1,nr)*WD))
    raw=exact.coarse(n,parameters,boundary);RC=raw[0,:,:,:rd+1,:rd+1];DC=raw[1]
    RV=RC.reshape(-1,27)*WR;U,sc,Vc=svd(RV,full_matrices=False,check_finite=False)
    cranks=rank_check(RV,'generic_protected_r_rank');nc=cranks[0]
    TC=Vc[:nc].T/sc[:nc];coarse_R=RC@TC;coarse_D=DC@TC
    coarse_D-=DK@(DKvalues.T@(coarse_D.reshape(-1,nc)*WD))
    overlap=(coarse_R.reshape(-1,nc)*WR).T@(R.reshape(-1,nr)*WR)
    RF=R-coarse_R@overlap;DF=DI-coarse_D@overlap
    RFvalues=RF.reshape(-1,nr)*WR;U,s,V=svd(RFvalues,full_matrices=False,check_finite=False)
    franks=rank_check(RFvalues,'generic_fine_r_rank');nf=franks[0]
    if nf!=nr-nc:raise ValueError('generic_position_split')
    TF=V[:nf].T/s[:nf];R=np.concatenate((coarse_R,RF@TF),axis=4);DI=np.concatenate((coarse_D,DF@TF),axis=4)
    coordinates=np.zeros((nr+nd,27));coordinates[:nc]=sc[:nc,None]*Vc[:nc]
    coordinates[nr:]=DKvalues.T@((DC-coarse_D@coordinates[:nc]).reshape(-1,27)*WD)
    hx,ht=hermite_nodes(n,rd);X,Y=np.meshgrid(hx,hx,indexing='ij');dx,dy=np.meshgrid(ht,ht,indexing='ij')
    hp=np.column_stack((X.ravel(),Y.ravel()));jets=evaluate_legendre(hp,n,R,True)
    which=np.where(dx.ravel(),np.where(dy.ravel(),4,1),np.where(dy.ravel(),2,0));control=np.zeros_like(jets[0])
    for j in [0,1,2,4]:control[which==j]=jets[j][which==j]
    control*=((1/(4*n))**(dx+dy)).ravel()[:,None]
    Rc=control.reshape(len(hx),len(hx),nr)
    gx=base.nodes(n,dd);X,Y=np.meshgrid(gx,gx,indexing='ij');dp=np.column_stack((X.ravel(),Y.ravel()))
    Dc=evaluate_legendre(dp,n,np.concatenate((DI,DK),axis=4))[0].reshape(len(gx),len(gx),nr+nd)
    if boundary=='clamped':Rc[0]=0;Dc[0]=0
    q,w=geometry.quadrature(n,order);CC=base.native(q,n)[2][0];indices=np.arange(CC.shape[1])
    if boundary=='clamped':indices=indices[n+1:]
    Tc,cc,ccond=base.independent(CC[:,indices],np.sqrt(w));nc2=Tc.shape[1]
    coordinates=np.vstack((coordinates,np.zeros((nc2,27))))
    q0=None
    if boundary=='extended':
        q0=coordinates[:,:3].copy()
        normal=reference_normal_coefficients(n,parameters,dd)
        q0[nr:nr+nd]+=DKvalues.T@(normal.reshape(-1,3)*WD)
    return dict(n=n,parameters=parameters,boundary=boundary,r_degree=rd,degree=dd,Rcoeff=Rc,Dcoeff=Dc,
        Cindices=indices,Tc=Tc,nr=nr,nd=nd,nc=nc2,ns=nr+nd+nc2,canonical=record,
        ranks=dict(r=rr,d=dr,c=cc),nr_coarse=nc,coarse_coordinates=coordinates[:,3:],reference_coordinates=q0)


def fields(model,p):
    n=model['n'];R=interpolation(p,n,model['r_degree'],model['Rcoeff'],True)
    P=interpolation(p,n,model['degree'],model['Dcoeff'])
    g=geometry.graph(p,model['parameters'])[4];dg=geometry.metric_derivatives(p,model['parameters'])
    D=[P[0]/g[:,None]**1.5]+[P[j+1]/g[:,None]**1.5-1.5*P[0]*dg[j][:,None]/g[:,None]**2.5 for j in range(2)]
    CC=[a[:,model['Cindices']]@model['Tc'] for a in c_shapes(p,n)]
    out=[]
    for j in range(3):
        a=np.zeros((3,len(p),model['ns']),dtype=np.result_type(p,float));a[0,:,:model['nr']]=R[j]
        a[1,:,:model['nr']+model['nd']]=D[j];a[2,:,model['nr']+model['nd']:]=CC[j];out.append(a)
    return out


def c_shapes(points,n):
    cell=np.minimum(np.maximum((points.real*n).astype(int),0),n-1);t=points*n-cell
    result=[np.zeros((len(points),(n+1)**2),dtype=np.result_type(points,float)) for _ in range(3)]
    for a in [0,1]:
        for b in [0,1]:
            ix=(cell[:,0]+a)*(n+1)+cell[:,1]+b
            x=t[:,0] if a else 1-t[:,0];y=t[:,1] if b else 1-t[:,1]
            result[0][np.arange(len(points)),ix]=x*y
            result[1][np.arange(len(points)),ix]=(n if a else -n)*y
            result[2][np.arange(len(points)),ix]=(n if b else -n)*x
    return result



def kinematics(model,p3,h,variable=False):
    p,back=np.unique(p3[:,:2],axis=0,return_inverse=True);fs=fields(model,p);X,J,det,H,Ha=geometry.reference(p3[:,:2],p3[:,2],h,model['parameters'],variable);z=H*p3[:,2]/2;za=Ha*p3[:,2,None]/2;parts=[v[:,back] for v in fs];N=parts[0][0]+z[:,None]*parts[0][1]+z[:,None]**2/2*parts[0][2];G=np.empty((len(p3),model['ns'],3))
    for a in range(2):G[:,:,a]=parts[a+1][0]+z[:,None]*parts[a+1][1]+z[:,None]**2/2*parts[a+1][2]+za[:,a,None]*(parts[0][1]+z[:,None]*parts[0][2])
    G[:,:,2]=H[:,None]/2*(parts[0][1]+z[:,None]*parts[0][2]);G=G@np.linalg.inv(J);MB=np.zeros((6,3,3));MB[0,0,0]=MB[1,1,1]=MB[2,2,2]=1
    for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2)
    B=np.einsum('kij,pnj->pkni',MB,G).reshape(len(p3),6,3*model['ns']);return N,G,B,det,X
