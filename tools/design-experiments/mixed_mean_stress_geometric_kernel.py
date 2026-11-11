"""Complete analytic point-isometry subspaces used as known null identities."""
from fractions import Fraction as F
import numpy as np
from scipy.linalg import qr,solve_triangular
import moment_enriched_exact_isometry_v2 as exact
p=exact;thin=exact.thin;c=thin.c
def local_known(q,k,R):
    poly,director,cert=exact.certificate(np.asarray(q,dtype=np.float64).tobytes(),k)
    plane=[c.LOOK[i,j,0] for i in range(6) for j in range(6)];dirrows=[c.LOOK[i,j,0] for i in range(6) for j in range(6) if (i,j)!=(5,5)];n=poly.shape[-1]
    latent=np.kron(R,np.eye(3))@poly[plane].reshape(108,n);Q,T=qr(latent,mode='full',check_finite=False);d=director[dirrows].reshape(105,n);d=solve_triangular(T[:n].T,d.T,lower=True,check_finite=False).T
    if not cert['all_rational_substitutions_zero']:raise ValueError('known_geometric_kernel_certificate')
    return Q[:,:n],Q[:,n:],d,cert
def sc(a,s):return {key:value*s for key,value in a.items() if value*s}
def sub(a,b):return p.add(a,b,F(-1))
def power(a,n):
    out={(0,0):F(1)}
    for _ in range(n):out=p.mul(out,a)
    return out
def compose(poly,x,y):
    out={}
    for (i,j),value in poly.items():out=p.add(out,p.mul(power(x,i),power(y,j)),value)
    return out
def physical_fields(k):
    x={(1,0):F(1)};y={(0,1):F(1)};out=[];names=[]
    if not k:
        for a in range(2,6):
            for b in range(6-a):
                w={(a,b):F(1)};out.append(([{}, {},w],[sc(p.deriv(w,0),F(-1,2)),sc(p.deriv(w,1),F(-1,2)),{}]));names.append(f'flat_P5_normal_{a}_{b}')
        return out,names
    kk=F(str(k))
    for family,degree in [('b',2),('b',3),('b',4),('a',4)]:
        a=power(x,degree) if family=='a' else {};b=power(x,degree) if family=='b' else {};ap=p.deriv(a,0);app=p.deriv(ap,0);appp=p.deriv(app,0);bp=p.deriv(b,0);cx={}
        for (i,j),value in p.mul(x,bp).items():cx[i+1,j]=-kk*value/(i+1)
        u=[p.add(p.mul(y,sub(p.mul(x,app),ap)),cx),a,sub(b,sc(p.mul(y,app),1/kk))]
        d=[sub(sc(p.mul(y,appp),1/kk),bp),p.add(sub(sc(p.mul(power(x,2),app),kk),sc(p.mul(x,ap),kk)),sc(app,1/kk)),sub(p.mul(p.mul(x,y),appp),sc(p.mul(x,bp),kk))]
        out.append((u,[sc(v,F(1,2)) for v in d]));names.append(f'curved_{family}_{degree}')
    return out,names
def native_global_fields(q,k,state):
    q=q.reshape(75,3);x=p.bilinear(q[[0,4,20,24],0]);y=p.bilinear(q[[0,4,20,24],1]);fields,names=physical_fields(k);coeff=[np.zeros((len(c.IDS),3,len(fields))) for _ in range(2)]
    for j,(u,d) in enumerate(fields):
        for order,parts in enumerate([u,d]):
            for cart,poly in enumerate(parts):
                for (a,b),value in compose(poly,x,y).items():
                    for i,vi in c.powers(a):
                        for l,vl in c.powers(b):
                            if i>5 or l>5:raise ValueError('global_polynomial_outside_position_space')
                            row=c.LOOK[i,l,order];coeff[order][row,cart,j]+=float(value*vi*vl)/c.NORM[row]
    native=[]
    for coef in coeff:
        v=np.zeros((state['nscalar'],3,len(fields)));v[:76]=coef[thin.f.VI]
        for face in state['faces']:v[face['offset']:face['offset']+face['n']]=np.einsum('mn,maj->naj',face['Mcoeff'],c.restrict(coef,face['axis'],face['sign']))
        native.append(v.reshape(-1,len(fields)))
    return np.vstack(native),names
