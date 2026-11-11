"""MMS-1 full q0 and full odd TB-kernel q1 necessary constraints."""
import numpy as np
from scipy.linalg import cholesky,solve_triangular,qr
import moment_enriched_thin_v2 as old
import mixed_mean_stress_space_v2 as stress
f=old.f;c=f.c;m=f.m;fit=old.fit;nullspace=old.nullspace;rank_info=old.rank_info
def scalar_position(s):
    ids=[v for v in f.PID if v[2]==0];coeff=np.zeros((len(c.IDS),36))
    for j,key in enumerate(ids):coeff[c.LOOK[key],j]=1.
    I=np.zeros((s['nscalar'],36));I[:76]=coeff[f.VI]
    for face in s['faces']:I[face['offset']:face['offset']+face['n']]=face['Mcoeff'].T@c.restrict(coeff,face['axis'],face['sign'])
    Q,R=qr(I,mode='economic',check_finite=False);return Q,R,coeff
def odd_TB(s,independent=False):
    cols=[];ns=s['nscalar']
    for i,key in enumerate(f.VID):
        if key[2]==1:v=np.zeros(ns);v[i]=1;cols.append(v)
    bot,top=s['faces'][:2]
    for j in range(25):
        v=np.zeros(ns);v[bot['offset']+j]=-1/np.sqrt(2);v[top['offset']+j]=1/np.sqrt(2);cols.append(v)
    for face in s['faces'][2:]:
        for j in range(face['n']):
            if j%3==1:v=np.zeros(ns);v[face['offset']+j]=1;cols.append(v)
    O=np.column_stack(cols);J=np.vstack([face['hat'].T@((face['w']/4)[:,None]*face['jump']) for face in s['faces'][:2]]);Z,info=nullspace(J@O,independent=independent)
    return O@Z,info
def local(q,kappa,order=13,zorder=7,independent=False):
    op=old.operators(q,order,zorder,independent);s=op['state'];ns=s['nscalar'];nc=3*ns
    if ns!=186:raise ValueError('MMS1_first_stage_requires_basic_benchmark_faces')
    Q,R,poly=scalar_position(s);O,oddinfo=odd_TB(s,independent);P0=np.kron(Q,np.eye(3));P1=np.kron(O,np.eye(3));count0=P0.shape[1]
    B=np.concatenate(((op['A'][1]@P0),(op['A'][0]@P1)),axis=2)/op['d0'][:,None,None]
    plane,wp=m.quadrature(order,2);z,wz=np.polynomial.legendre.leggauss(zorder);B=B.reshape(len(plane),zorder,6,-1);mean=np.einsum('z,pzkn->pkn',wz/2,B);fluct=float(np.linalg.norm(B-mean[:,None])/max(1.,np.linalg.norm(B)))
    L,meta=stress.evaluate(q,kappa,plane,independent);d0=op['d0'].reshape(len(plane),zorder)[:,0];weight=2*wp*d0;Ci=np.linalg.inv(m.C)
    gram=np.einsum('p,pir,ij,pjs->rs',weight,L,Ci,L,optimize=True);chol=cholesky(gram,lower=True,check_finite=False)
    moment=np.einsum('p,pkr,pkn->rn',weight,L,mean,optimize=True);A=solve_triangular(chol,moment,lower=True,check_finite=False);Zlatent,info=nullspace(A,independent=independent)
    if fluct>1e-9:raise ValueError('nonzero_thickness_fluctuation_constraint:'+str(fluct))
    Z=np.zeros((2*nc,Zlatent.shape[1]));Z[:nc]=P0@Zlatent[:count0];Z[nc:]=P1@Zlatent[count0:]
    native0=np.kron(solve_triangular(R,np.eye(36),check_finite=False),np.eye(3))@Zlatent[:count0]
    residual=float(np.linalg.norm(A@Zlatent)/max(1.,np.linalg.norm(A)))
    # Test all six physical constant stresses in the complete selected space.
    C0=np.tile(np.eye(6),(len(plane),1)).reshape(len(plane),6,6);F=L.reshape(-1,L.shape[-1]);coef=np.linalg.lstsq(F,C0.reshape(-1,6),rcond=1e-12)[0];constant=float(np.linalg.norm(F@coef-C0.reshape(-1,6))/np.sqrt(6*len(plane)))
    if constant>1e-9:raise ValueError('constant_physical_stress_not_preserved:'+str(constant))
    return dict(q=q,Z=Z,poly0=native0,polycoeff=poly,nc=nc,ni=3*s['ni'],nb=3*(ns-s['ni']),ny=Z.shape[1],nscalar=ns,info=info,stress=meta,odd_TB=oddinfo,fluctuation=fluct,residual=residual,constant_stress=constant,gram_condition=float(np.linalg.cond(gram)),faces=[dict(axis=x['axis'],sign=x['sign'],offset=x['offset'],n=x['n'],Mcoeff=x['Mcoeff']) for x in s['faces'][2:]])
def target(q,p,k):
    N,_=m.shape(p,0);x=(N@q.reshape(75,3)[:25])[:,0]
    b=np.column_stack((x**4/12,x**5,np.exp(x)-1-x));u=np.column_stack((-k*x**5/15,-5*k*x**6/6,-k*((x-1)*np.exp(x)+1-x*x/2)))
    out=np.zeros((len(p),3,3));out[:,0,:]=u;out[:,2,:]=b;return out
def sample(cell,k,order):
    plane,wp=m.quadrature(order,2);p=np.column_stack((plane,np.zeros(len(plane))));L,_=c.scalar_basis(p,order==15);field=L@cell['polycoeff'];N,g=m.shape(p,0,order==15);J=np.einsum('na,pnb->pab',cell['q'].reshape(75,3)[:25],g);area=np.linalg.norm(np.cross(J[:,:,0],J[:,:,1]),axis=1);weight=np.sqrt(wp*area)
    E=np.kron(field,np.eye(3))@cell['poly0'];return E* np.repeat(weight,3)[:,None],(target(cell['q'],p,k)*weight[:,None,None]).reshape(-1,3)
