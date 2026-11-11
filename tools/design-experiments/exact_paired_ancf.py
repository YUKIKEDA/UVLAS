"""Exact coefficient-space PG assembly followed by thickness-hierarchical coordinates.

FLINT rational arithmetic constructs the r image / pure-d kernel exactly.
Only the final, well-scaled polynomial basis is converted to float64. Rank
decisions are never made from cancellation-prone physical field samples.
"""
from pathlib import Path
import hashlib,json,time,os
import numpy as np
from scipy.linalg import qr,svd
from numpy.polynomial import legendre as L
from flint import fmpq_mat,arb_mat,ctx
import paired_exact_polynomials as exact
import coefficient_compatible_ancf as polynomial
geometry=polynomial.geometry
lift=polynomial.lift
compatible_particular=polynomial.compatible_particular
ROOT=Path(__file__).resolve().parents[2]


def rank_check(A,label):
    s=svd(A,compute_uv=False,check_finite=False)
    ranks=[int(sum(s>tol*max(1.,np.linalg.norm(A)))) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError(label+':'+str(ranks))
    return ranks


def canonical(n,parameters,arithmetic='arb192'):
    version=hashlib.sha256(Path(exact.__file__).read_bytes()).hexdigest()[:12]
    key=hashlib.sha256(repr((n,tuple(parameters),version,arithmetic)).encode()).hexdigest()[:16]
    directory=ROOT/'build/analysis/exact-paired-space';directory.mkdir(exist_ok=True)
    path=directory/(key+'.npz');meta=directory/(key+'.json')
    if path.exists() and meta.exists():
        record=json.loads(meta.read_text())
        if hashlib.sha256(path.read_bytes()).hexdigest()==record['artifact_sha256']:
            a=np.load(path);return a['coefficients'],dict(record,runtime_cache_hit=True)
    tick=time.monotonic();source=exact.build(n,parameters);A=source.reshape(-1,source.shape[-1])
    norm=np.linalg.norm(A,axis=1);rows=np.flatnonzero(norm);E=A[rows]/norm[rows,None]
    ss=svd(E,compute_uv=False,check_finite=False)
    ranks=[int(sum(ss>tol*max(1.,np.linalg.norm(E)))) for tol in [1e-9,1e-10,1e-11]]
    rsize=source[0].size//source.shape[-1];rrows=np.flatnonzero(norm[:rsize])
    ER=A[rrows]/norm[rrows,None];ss=svd(ER,compute_uv=False,check_finite=False)
    rr=[int(sum(ss>tol*max(1.,np.linalg.norm(ER)))) for tol in [1e-9,1e-10,1e-11]]
    rational=exact.rational_matrix(n,parameters,rows,np.arange(source.shape[-1]))
    exact_rank=int(rational.rank())
    position=exact.rational_matrix(n,parameters,rrows,np.arange(source.shape[-1]))
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
    for bits in ([192,256,320] if arithmetic=='arb192' else [0]):
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
    C=output.reshape(2,2*n,2*n,14,14,k)
    record=dict(n=n,parameters=parameters,exact_rank=k,exact_r_rank=nr,exact_d_kernel=k-nr,
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


def evaluate_legendre(p,n,coeff,second=False):
    degree=coeff.shape[2]-1;cell=np.minimum(np.maximum((p*2*n).astype(int),0),2*n-1)
    t=2*(p*2*n-cell)-1
    pairs=[(0,0),(1,0),(0,1)]+([(2,0),(1,1),(0,2)] if second else [])
    out=[np.empty((len(p),coeff.shape[-1])) for _ in pairs]
    derivatives=[L.legder(np.eye(degree+1),m=j,axis=0) for j in range(3)]
    for i,j in np.unique(cell,axis=0):
        mask=(cell[:,0]==i)&(cell[:,1]==j);x,y=t[mask].T
        X=[L.legvander(x,degree-d)@derivatives[d]*(4*n)**d for d in range(3)]
        Y=[L.legvander(y,degree-d)@derivatives[d]*(4*n)**d for d in range(3)]
        local=coeff[i,j].reshape((degree+1)**2,-1)
        for k,(a,b) in enumerate(pairs):out[k][mask]=np.einsum('pi,pj->pij',X[a],Y[b]).reshape(sum(mask),-1)@local
    return out


def prepare(n,parameters,order=15):
    C,record=canonical(n,parameters);nr=record['exact_r_rank'];nd=record['exact_d_kernel']
    a,b=np.meshgrid(np.arange(14),np.arange(14),indexing='ij')
    weight=np.broadcast_to(1/(2*n)/np.sqrt((2*a+1)*(2*b+1)),(2*n,2*n,14,14)).reshape(-1,1)
    R=C[0,:,:,:,:,:nr];DI=C[1,:,:,:,:,:nr];DK=C[1,:,:,:,:,nr:]
    RA=R.reshape(-1,nr)*weight;DA=DK.reshape(-1,nd)*weight
    U,s,V=svd(RA,full_matrices=False,check_finite=False);rr=rank_check(RA,'scaled_r_image_rank')
    if rr!=[nr]*3:raise ValueError('r_image_dimension')
    Tr=V.T/s;R=R@Tr;DI=DI@Tr
    U,s,V=svd(DA,full_matrices=False,check_finite=False);dr=rank_check(DA,'scaled_d_kernel_rank')
    if dr!=[nd]*3:raise ValueError('d_kernel_dimension')
    Td=V.T/s;DK=DK@Td
    DKvalues=DK.reshape(-1,nd)*weight;DIvalues=DI.reshape(-1,nr)*weight
    cross=DKvalues.T@DIvalues;DI-=DK@cross
    # Convert the completed polynomial basis to shared C1/C0 controls for evaluation.
    hx,ht=polynomial.hermite_nodes(n);X,Y=np.meshgrid(hx,hx,indexing='ij');dx,dy=np.meshgrid(ht,ht,indexing='ij')
    hp=np.column_stack((X.ravel(),Y.ravel()));jets=evaluate_legendre(hp,n,R,True)
    index=np.where(dx.ravel(),np.where(dy.ravel(),4,1),np.where(dy.ravel(),2,0));controls=np.zeros_like(jets[0])
    for j in [0,1,2,4]:controls[index==j]=jets[j][index==j]
    controls*=((1/(4*n))**(dx+dy)).ravel()[:,None]
    Rc=controls.reshape(len(hx),len(hx),nr);Rc[0]=0
    gx=polynomial.nodes(n,13);X,Y=np.meshgrid(gx,gx,indexing='ij');p=np.column_stack((X.ravel(),Y.ravel()))
    Dc=evaluate_legendre(p,n,np.concatenate((DI,DK),axis=4))[0].reshape(len(gx),len(gx),nr+nd);Dc[0]=0
    p,w=geometry.quadrature(n,order);nativeC=polynomial.native(p,n)[2][0];indices=np.arange(nativeC.shape[1])[n+1:]
    Tc,cr,cc=polynomial.independent(nativeC[:,indices],np.sqrt(w))
    return dict(n=n,parameters=parameters,Rcoeff=Rc,Dcoeff=Dc,degree=13,Cindices=indices,Tc=Tc,
        nr=nr,nd=nd,nc=Tc.shape[1],ns=nr+nd+Tc.shape[1],lift=dict(n=n,parameters=parameters),
        ranks=dict(r=rr,d=dr,c=cr),canonical=record)


def fields(model,p):
    n=model['n'];R=polynomial.interpolation(p,n,10,model['Rcoeff'],True)
    P=polynomial.interpolation(p,n,model['degree'],model['Dcoeff'])
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
