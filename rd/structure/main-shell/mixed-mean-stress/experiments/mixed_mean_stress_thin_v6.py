"""Analytic director/gauge separation and global position-kernel elimination."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import qr,svd,cholesky,solve_triangular
import mixed_mean_stress_thin_v5 as prior
from mixed_mean_stress_constraints import nullspace as original_constraints
old=prior.old;f=prior.f;c=prior.c;m=prior.m;fit=prior.fit;stress=prior.stress;nullspace=prior.nullspace;rank_info=prior.rank_info
scalar_position=prior.scalar_position;target=prior.target;sample=prior.sample
def reduced_nullspace(A,independent=False):
    sv=svd(A,compute_uv=False,check_finite=False);raw=rank_info(sv)
    if len(set(raw))==1:
        Z,info=nullspace(A,independent=independent);info.update(raw_ranks=raw,block_equilibrated=False);return Z,info
    if A.shape[0]!=72:raise ValueError('MMS1_unexpected_reduced_rank:'+str(raw))
    _,_,p=qr(A,mode='economic',pivoting=True,check_finite=False);strong=p[:69];free=p[69:];G=A[:,strong];sr=rank_info(svd(G,compute_uv=False,check_finite=False))
    if sr!=[69]*3:raise ValueError('MMS1_membrane_basic_rank:'+str(sr))
    Q,R=qr(G,mode='full',check_finite=False);top=Q[:,:69].T@A[:,free];weak=Q[:,69:].T@A[:,free];scale=np.linalg.norm(weak);Zf,wi=nullspace(weak/scale,independent=independent);Z=np.zeros((A.shape[1],Zf.shape[1]));Z[free]=Zf;Z[strong]=-solve_triangular(R[:69],top@Zf,check_finite=False);Z=qr(Z,mode='economic',check_finite=False)[0]
    residual=float(np.linalg.norm(A@Z)/max(1.,np.linalg.norm(A)));wr=float(np.linalg.norm(weak@Z[free])/scale)
    if max(residual,wr)>1e-9:raise ValueError('MMS1_reduced_block_residual')
    rank=69+wi['ranks'][0];return Z,dict(raw_ranks=raw,ranks=[rank]*3,columns=A.shape[1],nullity=Z.shape[1],residual=residual,block_equilibrated=True,strong_ranks=sr,completion=wi,completion_scale=float(scale),completion_residual=wr)
def local(q,kappa,order=13,zorder=7,independent=False):
    op=old.operators(q,order,zorder,independent);s=op['state'];ns=s['nscalar'];nc=3*ns
    if ns!=186:raise ValueError('MMS1_basic_faces_required')
    Q0,R,poly=scalar_position(s);P0=np.kron(Q0,np.eye(3));rinv=solve_triangular(R,np.eye(36),check_finite=False)
    # Complete polynomial director lift Ih(zeta*psi) and its exact zero image.
    keys=[(i,j,0) for i in range(6) for j in range(6) if (i,j)!=(5,5)];coef=np.zeros((len(c.IDS),35));plane_coef=np.zeros_like(coef)
    for j,(a,b,z) in enumerate(keys):coef[c.LOOK[a,b,1],j]=1/np.sqrt(3);plane_coef[c.LOOK[a,b,0],j]=1.
    Idir=np.zeros((ns,35));Idir[:76]=coef[f.VI]
    for face in s['faces']:Idir[face['offset']:face['offset']+face['n']]=face['Mcoeff'].T@c.restrict(coef,face['axis'],face['sign'])
    Pd=np.kron(Idir,np.eye(3));O,odd=prior.odd_TB(s,independent);Gnormal=s['bodycoeff'][:,2,:]@s['D'];Zh,hominfo=nullspace(Gnormal@O,independent=independent);Ph=np.kron(O@Zh,np.eye(3))
    if hominfo['nullity']!=10:raise ValueError('MMS1_odd_director_homogeneous_dimension')
    # Analytic cylinder surface gradients. No subtraction of a q1-zero image
    # from nearly dependent position columns is used in the constraint matrix.
    plane,wp=m.quadrature(order,2);p3=np.column_stack((plane,np.zeros(len(plane))));L,dL=c.scalar_basis(p3,independent);phi=L@plane_coef;gp=dL@poly
    _,_,den,meta_dummy=stress.coefficients(np.asarray(q,dtype=np.float64).tobytes(),kappa,independent) if hasattr(stress,'coefficients') else (None,None,None,None)
    corners=q.reshape(75,3)[[0,4,20,24]];xp=stress.p.bilinear(corners[:,0]);yp=stress.p.bilinear(corners[:,1])
    def evaluate(poly):return sum(float(v)*plane[:,0]**i*plane[:,1]**j for (i,j),v in poly.items()) if poly else np.zeros(len(plane))
    x=evaluate(xp);xx=evaluate(stress.p.deriv(xp,0));xy=evaluate(stress.p.deriv(xp,1));yx=evaluate(stress.p.deriv(yp,0));yy=evaluate(stress.p.deriv(yp,1));jac=xx*yy-xy*yx;g=1+kappa*kappa*x*x;d=np.column_stack((-kappa*x,np.zeros(len(x)),np.ones(len(x))));tx=np.column_stack((np.ones(len(x)),np.zeros(len(x)),kappa*x))
    gx=(yy[:,None]*gp[:,0]-yx[:,None]*gp[:,1])/jac[:,None];gy=(-xy[:,None]*gp[:,0]+xx[:,None]*gp[:,1])/jac[:,None];grad=tx[:,:,None]*gx[:,None,:]/g[:,None,None];grad[:,1]+=gy;grad=grad@rinv
    B0=np.einsum('kab,pbn->pkna',m.MB,grad).reshape(len(plane),6,108);Bd=2*np.einsum('kab,pb,pn->pkna',m.MB,d/g[:,None],phi).reshape(len(plane),6,105)
    Ls,meta=stress.evaluate(q,kappa,plane,independent);weight=wp*g*jac;gram=np.einsum('p,pir,ij,pjs->rs',weight,Ls,np.linalg.inv(m.C),Ls,optimize=True);chol=cholesky(gram,lower=True,check_finite=False)
    C0=solve_triangular(chol,np.einsum('p,pkr,pkn->rn',weight,Ls,B0,optimize=True),lower=True,check_finite=False);Cd=solve_triangular(chol,np.einsum('p,pkr,pkn->rn',weight,Ls,Bd,optimize=True),lower=True,check_finite=False)
    dr=rank_info(svd(Cd,compute_uv=False,check_finite=False))
    if dr!=[105]*3:raise ValueError('MMS1_director_work_rank:'+str(dr))
    Qd,Rd=qr(Cd,mode='full',check_finite=False);K=Qd[:,105:].T@C0;Z0,info=reduced_nullspace(K,independent);director=-solve_triangular(Rd[:105],Qd[:,:105].T@C0@Z0,check_finite=False);n0=Z0.shape[1];nh=Ph.shape[1]
    Z=np.zeros((2*nc,n0+nh));Z[:nc,:n0]=P0@Z0;Z[nc:,:n0]=Pd@director;Z[nc:,n0:]=Ph
    poly0=np.zeros((108,n0+nh));poly0[:,:n0]=np.kron(rinv,np.eye(3))@Z0;residual=float(np.linalg.norm(C0@Z0+Cd@director)/max(1.,np.linalg.norm(C0),np.linalg.norm(Cd)))
    # Original float64 weak-gradient equations independently observe the lift.
    original=(op['A'][1]@Z[:nc]+op['A'][0]@Z[nc:])/op['d0'][:,None,None];nz=zorder;_,wz=np.polynomial.legendre.leggauss(nz);orig=original.reshape(len(plane),nz,6,-1);mean=np.einsum('z,pzkn->pkn',wz/2,orig);fluct=float(np.linalg.norm(orig-mean[:,None])/max(1.,np.linalg.norm(orig)));work=solve_triangular(chol,np.einsum('p,pkr,pkn->rn',weight,Ls,mean,optimize=True),lower=True);original_res=float(np.linalg.norm(work)/max(1.,np.linalg.norm(Z)))
    F=Ls.reshape(-1,Ls.shape[-1]);const=np.tile(np.eye(6),(len(plane),1)).reshape(-1,6);constant=float(np.linalg.norm(F@np.linalg.lstsq(F,const,rcond=1e-12)[0]-const)/np.sqrt(6*len(plane)))
    if max(residual,fluct,original_res,constant)>1e-9:raise ValueError('MMS1_separated_original_work_residual:'+str((residual,fluct,original_res,constant)))
    info.update(q0_ranks=info['ranks'],ranks=[v+105 for v in info['ranks']],director_ranks=dr,q0_dimension=n0,q1_homogeneous_dimension=nh,original_work_residual=original_res)
    return dict(q=q,Z=Z,poly0=poly0,polycoeff=poly,nc=nc,ni=3*s['ni'],nb=3*(ns-s['ni']),ny=n0+nh,n0=n0,nh=nh,nscalar=ns,info=info,stress=meta,odd_TB=odd,homogeneous=hominfo,fluctuation=fluct,residual=residual,constant_stress=constant,gram_condition=float(np.linalg.cond(gram)),faces=[dict(axis=f['axis'],sign=f['sign'],offset=f['offset'],n=f['n'],Mcoeff=f['Mcoeff']) for f in s['faces'][2:]])
def global_position_kernel(cells,n):
    nx=sum(c['n0'] for c in cells);nh=sum(c['nh'] for c in cells);xo=ho=0;registry={};g0=[];gx=[];gh=[]
    for index,cell in enumerate(cells):
        cell['xoffset']=xo;cell['hoffset']=ho;xo+=cell['n0'];ho+=cell['nh']
        for face,key in zip(cell['faces'],fit.edge_keys(index,n)):
            sl=slice(3*face['offset'],3*(face['offset']+face['n']));a=cell['Z'][sl,:cell['n0']];b=cell['Z'][cell['nc']+sl.start:cell['nc']+sl.stop,:cell['n0']];h=cell['Z'][cell['nc']+sl.start:cell['nc']+sl.stop,cell['n0']:]
            # Basic face spaces use the same normalized polynomial coordinates.
            if key in registry:
                priorcell,pa,pb,ph=registry[key];A=np.zeros((len(a),nx));B=np.zeros_like(A);H=np.zeros((len(a),nh));xs=slice(cell['xoffset'],xo);hs=slice(cell['hoffset'],ho);ps=slice(priorcell['xoffset'],priorcell['xoffset']+priorcell['n0']);qs=slice(priorcell['hoffset'],priorcell['hoffset']+priorcell['nh']);A[:,xs]=a;A[:,ps]-=pa;B[:,xs]=b;B[:,ps]-=pb;H[:,hs]=h;H[:,qs]-=ph;g0.append(A);gx.append(B);gh.append(H)
            else:registry[key]=(cell,a,b,h)
            if key[0]==0 and key[1]==0:
                A=np.zeros((len(a),nx));B=np.zeros_like(A);H=np.zeros((len(a),nh));A[:,cell['xoffset']:xo]=a;B[:,cell['xoffset']:xo]=b;H[:,cell['hoffset']:ho]=h;g0.append(A);gx.append(B);gh.append(H)
    C0=np.vstack(g0);C1=np.vstack(gx);H=np.vstack(gh);U,sv,V=svd(H,full_matrices=False,check_finite=False);r=rank_info(sv)
    if len(set(r))!=1:raise ValueError('MMS1_shared_homogeneous_rank:'+str(r))
    U=U[:,:r[0]];reduced=C1-U@(U.T@C1);C=np.vstack((C0,reduced));Z,info=nullspace(C)
    return Z,dict(info=info,homogeneous_ranks=r,nx=nx,nh=nh,H=H,C1=C1,U=U,singular=sv[:r[0]],V=V[:r[0]])
def recover(cells,globaldata,x):
    h=-(globaldata['V'].T/globaldata['singular'])@(globaldata['U'].T@(globaldata['C1']@x));values=[]
    for cell in cells:
        a=x[cell['xoffset']:cell['xoffset']+cell['n0']];b=h[cell['hoffset']:cell['hoffset']+cell['nh']];values.append(cell['Z']@np.vstack((a,b)))
    return np.array(values)
