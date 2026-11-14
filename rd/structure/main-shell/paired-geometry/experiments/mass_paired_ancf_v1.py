"""Verified geometric L2 coordinates for the exact polynomial PG space.

Exact RREF chooses a basis without floating pivot/rank ambiguity. A geometric
reference-density mass inner product then conditions the SAME space. All
coefficient combinations are completed with ball arithmetic before float64
evaluation. The global smooth subspace is ordered first and preserved.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from functools import lru_cache
import hashlib,json,time
import numpy as np
from flint import fmpq_mat,arb_mat,arb,fmpq,ctx
import general_paired_polynomials as exact
import polynomial_reference_surface as geometry
import general_paired_ancf as polynomial
ROOT=rd_paths.ROOT


def fraction(x):return fmpq(x.numerator,x.denominator)


def pivots(A,rank):
    return [next(j for j in range(A.ncols()) if A[i,j]!=0) for i in range(rank)]


@lru_cache(maxsize=2)
def source_basis(n,parameters,boundary):
    tick=time.monotonic();C=exact.build(n,parameters,boundary);flat=C.reshape(-1,C.shape[-1])
    rows=np.flatnonzero(np.any(flat!=0,axis=1));Q=exact.rational_matrix(n,parameters,rows,np.arange(flat.shape[1]),boundary)
    RR,rank=Q.transpose().rref();pr=pivots(RR,rank)
    source=fmpq_mat([[RR[j,i] for j in range(rank)] for i in range(len(rows))])
    raw_coarse=exact.coarse_source_map(n,parameters,boundary)
    Qp=fmpq_mat([[Q[int(i),j] for j in range(Q.ncols())] for i in pr])
    coordinates=Qp*raw_coarse
    augmented=fmpq_mat(rank,30+rank)
    for i in range(rank):
        for j in range(30):augmented[i,j]=coordinates[i,j]
        augmented[i,30+i]=1
    reduced,full_rank=augmented.rref();pc=pivots(reduced,full_rank)
    if full_rank!=rank:raise ValueError('coarse_extension_rank')
    coarse=[j for j in pc if j<30];fine=[j-30 for j in pc if j>=30];nc=len(coarse)
    missing=sorted(set(range(rank))-set(fine))
    G=fmpq_mat([[coordinates[i,j] for j in coarse] for i in missing])
    coarse_coordinates=G.inv()*fmpq_mat([[coordinates[i,j] for j in range(30)] for i in missing])
    source_coarse=source*fmpq_mat([[coordinates[i,j] for j in coarse] for i in range(rank)])
    basis=fmpq_mat(len(rows),rank)
    for i in range(len(rows)):
        for j in range(nc):basis[i,j]=source_coarse[i,j]
        for j,k in enumerate(fine):basis[i,nc+j]=source[i,k]
    # Exact check: every mandatory field has coordinates only in the first block.
    Z=fmpq_mat(rank,30)
    for i in range(nc):
        for j in range(30):Z[i,j]=coarse_coordinates[i,j]
    if basis*Z!=Q*raw_coarse:raise ValueError('exact_coarse_inclusion')
    rd,dd=exact.specification(parameters)[1:]
    record=dict(n=n,parameters=parameters,boundary=boundary,exact_rank=int(rank),coarse_rank=nc,
                exact_coarse_inclusion=True,source_seconds=time.monotonic()-tick,r_degree=rd,d_degree=dd)
    print('EXACT_RREF',record,flush=True)
    return rows,basis,Z,record


def legendre_values(x,degree):
    result=[arb(1)]
    if degree:result.append(x)
    for j in range(1,degree):result.append(((2*j+1)*x*result[-1]-j*result[-2])/(j+1))
    return result


def evaluate_geometry(x,y,parameters):
    coefficients=geometry.coefficients(parameters);width=arb(fraction(geometry.chart(parameters)[1]))
    def evaluate(dx=0,dy=0):
        out=arb(0)
        for (i,j),a in coefficients.items():
            if i<dx or j<dy:continue
            multiplier=fraction(a)
            for k in range(dx):multiplier*=i-k
            for k in range(dy):multiplier*=j-k
            out+=arb(multiplier)*x**(i-dx)*y**(j-dy)
        return out
    skew=arb(fraction(geometry.shear(parameters)));fx=evaluate(1,0)/width;fy=evaluate(0,1)-skew*fx;g=1+fx*fx+fy*fy
    fxx=evaluate(2,0)/(width*width);fxy=evaluate(1,1)/width-skew*fxx
    fyy=evaluate(0,2)-2*skew*evaluate(1,1)/width+skew*skew*fxx
    trace=((1+fy*fy)*fxx-2*fx*fy*fxy+(1+fx*fx)*fyy)/(g*g.sqrt())
    gauss=(fxx*fyy-fxy*fxy)/(g*g)
    return g,width*g.sqrt(),trace,gauss


def c_polynomials(n,cx,cy,indices):
    from fractions import Fraction as F
    x={(0,0):F(2*cx+1,4*n),(1,0):F(1,4*n)}
    y={(0,0):F(2*cy+1,4*n),(0,1):F(1,4*n)}
    cellx=cx//2;celly=cy//2
    X=[exact.p.add({(0,0):F(cellx+1)},exact.p.scale(x,-n)),exact.p.add(exact.p.scale(x,n),{(0,0):F(-cellx)})]
    Y=[exact.p.add({(0,0):F(celly+1)},exact.p.scale(y,-n)),exact.p.add(exact.p.scale(y,n),{(0,0):F(-celly)})]
    lookup={int(j):k for k,j in enumerate(indices)};out={}
    for a in [0,1]:
        for b in [0,1]:
            index=(cellx+a)*(n+1)+celly+b
            if index not in lookup:continue
            for key,value in exact.p.to_legendre(exact.p.mul(X[a],Y[b]),1).items():
                out.setdefault(key,{})[lookup[index]]=fraction(value)
    return out


def field_coefficients(n,rows,basis,degree,cx,cy,field,columns):
    result={};size=degree+1
    for i,index in enumerate(rows):
        b=int(index)%size;q=int(index)//size;a=q%size;q//=size;y=q%(2*n);q//=2*n;x=q%(2*n);kind=q//(2*n)
        if kind!=field or x!=cx or y!=cy:continue
        values=[basis[i,j] for j in range(columns)]
        if any(values):result[(a,b)]=values
    return result


def eval_matrix(points,cell,n,degree,coeffs,columns):
    cx,cy=cell;keys=list(coeffs)
    if not keys:return arb_mat(len(points),columns)
    values=arb_mat(len(points),len(keys));coefficient=arb_mat([[coeffs[key][j] for j in range(columns)] for key in keys])
    for row,(x,y) in enumerate(points):
        tx=4*n*x-(2*cx+1);ty=4*n*y-(2*cy+1)
        X=legendre_values(tx,degree);Y=legendre_values(ty,degree)
        for j,(a,b) in enumerate(keys):values[row,j]=X[a]*Y[b]
    return values*coefficient


def weighted(matrix,weights):
    return arb_mat([[matrix[i,j]*weights[i] for j in range(matrix.ncols())] for i in range(matrix.nrows())])


def cholesky(M):
    count=M.nrows();L=[[arb(0) for _ in range(count)] for _ in range(count)]
    for i in range(count):
        for j in range(i+1):
            value=M[i,j]-sum((L[i][k]*L[j][k] for k in range(j)),arb(0))
            if i==j:
                if not value>0:raise ValueError('uncertified_positive_mass_pivot')
                L[i][j]=value.sqrt()
            else:L[i][j]=value/L[j][j]
    return arb_mat(L)


def geometric_mass(n,parameters,h,variable,boundary,rows,basis,degree,indices,order):
    k=basis.ncols();nc=len(indices);M=arb_mat(k+nc,k+nc)
    nodes,weights=np.polynomial.legendre.leggauss(order);active_counts=[]
    for cx in range(2*n):
        for cy in range(2*n):
            points=[];wp=[]
            for i in range(order):
                for j in range(order):
                    points.append(((arb(cx)+(arb(float(nodes[i]))+1)/2)/(2*n),(arb(cy)+(arb(float(nodes[j]))+1)/2)/(2*n)))
                    wp.append(arb(float(weights[i]))*arb(float(weights[j]))/(16*n*n))
            rc=field_coefficients(n,rows,basis,degree,cx,cy,0,k)
            dc=field_coefficients(n,rows,basis,degree,cx,cy,1,k)
            active=[j for j in range(k) if any(a[j]!=0 for a in rc.values()) or any(a[j]!=0 for a in dc.values())]
            na=len(active);rc={key:[a[j] for j in active] for key,a in rc.items()};dc={key:[a[j] for j in active] for key,a in dc.items()}
            R=eval_matrix(points,(cx,cy),n,degree,rc,na);D=eval_matrix(points,(cx,cy),n,degree,dc,na)
            cp=c_polynomials(n,cx,cy,indices)
            active_c=sorted({j for values in cp.values() for j in values});nac=len(active_c)
            cdict={key:[values.get(j,fmpq(0)) for j in active_c] for key,values in cp.items()}
            C=eval_matrix(points,(cx,cy),n,1,cdict,nac);active_counts.append(na+nac)
            moments=[[] for _ in range(5)]
            for i,(x,y) in enumerate(points):
                g,det0,trace,gauss=evaluate_geometry(x,y,parameters)
                H=arb(h)*(1+arb(.1)*x+arb(.05)*y) if variable else arb(h)
                def integral(power):return arb(0) if power%2 else H**(power+1)/(2**power*(power+1))
                for power in range(5):moments[power].append(wp[i]*det0*(integral(power)-trace*integral(power+1)+gauss*integral(power+2)))
                divisor=g*g.sqrt()
                for j in range(na):D[i,j]/=divisor
            RD=R.transpose()*weighted(D,moments[1])
            block=R.transpose()*weighted(R,moments[0])+D.transpose()*weighted(D,moments[2])+RD+RD.transpose()
            cross=R.transpose()*weighted(C,[a/2 for a in moments[2]])+D.transpose()*weighted(C,[a/2 for a in moments[3]])
            cc=C.transpose()*weighted(C,[a/4 for a in moments[4]])
            for i,gi in enumerate(active):
                for j,gj in enumerate(active):M[gi,gj]+=block[i,j]
                for j,cj in enumerate(active_c):M[gi,k+cj]+=cross[i,j];M[k+cj,gi]+=cross[i,j]
            for i,ci in enumerate(active_c):
                for j,cj in enumerate(active_c):M[k+ci,k+cj]+=cc[i,j]
    return M,dict(max_active=max(active_counts),mean_active=sum(active_counts)/len(active_counts),total_columns=k+nc)


def prepare(n,parameters,h=.01,variable=False,boundary='clamped',order=None):
    parameters=tuple(float(v) for v in parameters);tick=time.monotonic()
    rows,basis,coordinates,record=source_basis(n,parameters,boundary)
    rd=record['r_degree'];dd=record['d_degree'];k=basis.ncols();nc0=record['coarse_rank']
    indices=np.arange((n+1)**2)
    if boundary=='clamped':indices=indices[n+1:]
    nc=len(indices);ns=k+nc
    if order is None:order=dd+2
    ctx.prec=768;ctx.threads=4
    M,locality=geometric_mass(n,parameters,h,variable,boundary,rows,basis,dd,indices,order)
    print('VERIFIED_MASS',n,boundary,ns,time.monotonic()-tick,flush=True)
    Lc=cholesky(M);T=Lc.inv().transpose()
    protected=fmpq_mat(ns,30)
    for i in range(k):
        for j in range(30):protected[i,j]=coordinates[i,j]
    Q= Lc.transpose()*arb_mat(protected)
    coarse_coordinates=np.array([[float(Q[i,j]) for j in range(30)] for i in range(ns)])
    R=np.zeros((2*n,2*n,rd+1,rd+1,ns));D=np.zeros((2*n,2*n,dd+1,dd+1,ns));C=np.zeros((2*n,2*n,2,2,ns));radius=0.
    for cx in range(2*n):
        for cy in range(2*n):
            for field,output,degree in [(0,R,rd),(1,D,dd),(2,C,1)]:
                if field<2:
                    source=field_coefficients(n,rows,basis,dd,cx,cy,field,k)
                    coeff={key:values+[fmpq(0)]*nc for key,values in source.items()}
                else:
                    source=c_polynomials(n,cx,cy,indices)
                    coeff={key:[fmpq(0)]*k+[values.get(j,fmpq(0)) for j in range(nc)] for key,values in source.items()}
                keys=list(coeff)
                if not keys:continue
                transformed=arb_mat([coeff[key] for key in keys])*T
                for i,(a,b) in enumerate(keys):
                    if a>degree or b>degree:raise ValueError('polynomial_degree')
                    for j in range(ns):
                        value=transformed[i,j]
                        if not value.is_finite():raise ValueError('nonfinite_normalized_coefficient')
                        radius=max(radius,float(value.rad()));output[cx,cy,a,b,j]=float(value)
    if radius>1e-25:raise ValueError('normalized_coefficient_interval:'+str(radius))
    # Shared controls retain exact C1/C0 conformity after the float64 conversion.
    hx,ht=polynomial.hermite_nodes(n,rd);X,Y=np.meshgrid(hx,hx,indexing='ij');dx,dy=np.meshgrid(ht,ht,indexing='ij')
    points=np.column_stack((X.ravel(),Y.ravel()));jets=polynomial.evaluate_legendre(points,n,R,True)
    which=np.where(dx.ravel(),np.where(dy.ravel(),4,1),np.where(dy.ravel(),2,0));control=np.zeros_like(jets[0])
    for j in [0,1,2,4]:control[which==j]=jets[j][which==j]
    control*=((1/(4*n))**(dx+dy)).ravel()[:,None];Rc=control.reshape(len(hx),len(hx),ns)
    def nodal(coeff,degree):
        nodes=polynomial.base.nodes(n,degree);X,Y=np.meshgrid(nodes,nodes,indexing='ij')
        points=np.column_stack((X.ravel(),Y.ravel()))
        if degree>=2:values=polynomial.evaluate_legendre(points,n,coeff)[0]
        else:
            cells=np.minimum((points*2*n).astype(int),2*n-1);t=2*(points*2*n-cells)-1
            values=np.empty((len(points),ns))
            for i,j in np.unique(cells,axis=0):
                mask=(cells[:,0]==i)&(cells[:,1]==j)
                values[mask]=np.einsum('pa,pb,abk->pk',np.polynomial.legendre.legvander(t[mask,0],degree),np.polynomial.legendre.legvander(t[mask,1],degree),coeff[i,j])
        return values.reshape(len(nodes),len(nodes),ns)
    Dc=nodal(D,dd);Cc=nodal(C,1)
    if boundary=='clamped':Rc[0]=0;Dc[0]=0;Cc[0]=0
    q0=coarse_coordinates[:,:3]+coarse_coordinates[:,27:30] if boundary=='extended' else None
    return dict(n=n,parameters=parameters,h=h,variable=variable,boundary=boundary,ns=ns,r_degree=rd,degree=dd,
        Rcoeff=Rc,Dcoeff=Dc,Ccoeff=Cc,coarse_coordinates=coarse_coordinates[:,3:27],reference_coordinates=q0,
        coarse_rank=nc0,exact_source=record,normalization=dict(precision_bits=768,max_radius=radius,quadrature=order,
        seconds=time.monotonic()-tick,verified_positive_pivots=ns,locality=locality),source_mass_cholesky=np.array([[float(Lc[i,j]) for j in range(ns)] for i in range(ns)]))


def fields(model,points):
    n=model['n'];R=polynomial.interpolation(points,n,model['r_degree'],model['Rcoeff'],True)
    P=polynomial.interpolation(points,n,model['degree'],model['Dcoeff']);C=polynomial.interpolation(points,n,1,model['Ccoeff'])
    g=geometry.graph(points,model['parameters'])[4];dg=geometry.metric_derivatives(points,model['parameters'])
    D=[P[0]/g[:,None]**1.5]+[P[j+1]/g[:,None]**1.5-1.5*P[0]*dg[j][:,None]/g[:,None]**2.5 for j in range(2)]
    return [np.stack((R[j],D[j],C[j])) for j in range(3)]


def kinematics(model,p3,h,variable=False):
    if h!=model['h'] or variable!=model['variable']:raise ValueError('reference_thickness_identity')
    p,back=np.unique(p3[:,:2],axis=0,return_inverse=True);F=fields(model,p)
    X,J,det,H,Ha=geometry.reference(p3[:,:2],p3[:,2],h,model['parameters'],variable)
    z=H*p3[:,2]/2;za=Ha*p3[:,2,None]/2;parts=[a[:,back] for a in F]
    N=parts[0][0]+z[:,None]*parts[0][1]+z[:,None]**2/2*parts[0][2]
    G=np.empty((len(p3),model['ns'],3),dtype=np.result_type(p3,float))
    for a in range(2):G[:,:,a]=parts[a+1][0]+z[:,None]*parts[a+1][1]+z[:,None]**2/2*parts[a+1][2]+za[:,a,None]*(parts[0][1]+z[:,None]*parts[0][2])
    G[:,:,2]=H[:,None]/2*(parts[0][1]+z[:,None]*parts[0][2]);G=G@np.linalg.inv(J)
    MB=np.zeros((6,3,3));MB[0,0,0]=MB[1,1,1]=MB[2,2,2]=1
    for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):MB[k,i,j]=MB[k,j,i]=1/np.sqrt(2)
    B=np.einsum('kij,pnj->pkni',MB,G).reshape(len(p3),6,3*model['ns'])
    return N,G,B,det,X
