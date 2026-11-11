"""Shared zero/first maps with rectangular-H1 nullspace treatment."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,itertools
import numpy as np
from scipy.linalg import solve,svd
from scipy.optimize import brentq
import mq3i1_shared_family as family
import check_mq3_zero_first_map as checks
import q3_thin_expansion as common
ROOT=rd_paths.ROOT
checks.family=family
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_cache():
    family.configure('covariant',True)
    family.register()
    return {}

def digest_bytes(q):return hashlib.sha256(q.tobytes()+b'M-Q3I1-v2').hexdigest()
def make_map(q):
    state=family.choose(q);X0,S,b,F,g=family.coefficients(state,0.);Sc,bc,Fc,gc=family.system(state,1e-25j);Xd=solve(S,bc.imag/1e-25-Sc.imag/1e-25@X0,assume_a='sym',check_finite=False);near=[]
    for h in [0.,-.0001,.0001]:
        br,cr,span=checks.ranks_at(state,h);X,_,_,Fh,gh=family.coefficients(state,h);res=float(np.linalg.norm(Fh@(state['T']@X[:state['T'].shape[1]])-gh)/max(1.,np.linalg.norm(gh)));near.append(dict(h=h,basis_ranks=br,constraint_ranks=cr,span=span,residual=res))
    p,w=family.base.quadrature(11);d0=family.data(q,0.,p);dc=family.data(q,1e-25j,p);R0=family.evaluate(state,0.,X0,p);Rd=family.evaluate(state,1e-25j,X0+1e-25j*Xd,p).imag/1e-25;Xc,_,_,_,_=family.coefficients(state,1e-25j);Ralt=family.evaluate(state,1e-25j,Xc,p).imag/1e-25;weight=w*np.linalg.det(d0['J']);weight/=weight.sum();wc=w*np.linalg.det(dc['J']);wc/=wc.sum();wd=wc.imag/1e-25;zp,zw=np.polynomial.legendre.leggauss(11);shape0=R0.reshape(11,11,11,6,150);shape1=Rd.reshape(11,11,11,6,150);mean=np.einsum('z,abzkn->abkn',zw/2,shape0);linear=2*np.einsum('z,z,abzkn->abkn',zw,zp,shape1)/np.sum(zw*zp*zp);zerror=checks.relative((shape0-mean[:,:,None]).reshape(-1,6,150),R0,weight);lerror=checks.relative((shape1-zp[None,None,:,None,None]*linear[:,:,None]/2).reshape(-1,6,150),Rd,weight);constant=max(float(np.max(abs(R0[:,:,144:]-np.eye(6)))),float(np.max(abs(Rd[:,:,144:]))));cerror=float(np.max(abs(R0[:,:,common.CC])));work0=np.einsum('p,pkn->kn',weight,R0-d0['E']);work1=np.einsum('p,pkn->kn',weight,Rd-dc['E'].imag/1e-25)+np.einsum('p,pkn->kn',wd,R0-d0['E']);work=max(float(np.linalg.norm(work0)),float(np.linalg.norm(work1)))/max(1.,float(np.linalg.norm(np.einsum('p,pkn->kn',weight,d0['E']))));coef=state['T']@X0[:state['T'].shape[1]];coefd=state['T']@Xd[:state['T'].shape[1]];moment=max(float(np.linalg.norm(F@coef-g)/max(1.,np.linalg.norm(g))),float(np.linalg.norm(Fc.imag/1e-25@coef+F@coefd-gc.imag/1e-25)/max(1.,np.linalg.norm(gc.imag/1e-25))));rigid=0.
    for axis in np.eye(3):
        trans=np.zeros((16,3,3));trans[:,0]=axis
        for v in [trans.reshape(144),np.cross(np.broadcast_to(axis,(16,3,3)),q.reshape(16,3,3)).reshape(144)]:
            for B in [R0[:,:,:144],Rd[:,:,:144]]:rigid=max(rigid,float(np.linalg.norm(B@v)/(max(1.,np.linalg.norm(B))*np.linalg.norm(v))))
    valid=dict(rank=all(v['basis_ranks']==state['space_ranks']==[state['space_ranks'][1]]*3 and v['constraint_ranks']==state['constraint_ranks']==[state['constraint_ranks'][1]]*3 for v in near),span=max(v['span'] for v in near)<=1e-10,constraints=max(v['residual'] for v in near)<=1e-10,zero_dependence=zerror<=1e-9,linear_dependence=lerror<=1e-9,constant=constant<=1e-9,zero_c=cerror<=1e-9,work=(not family.WORK or work<=1e-9),moments=moment<=1e-9,rigid=rigid<=1e-9,derivative=checks.relative(Rd-Ralt,Rd,weight)<=1e-9,solve=float(np.linalg.norm(S@X0-b)/max(1.,np.linalg.norm(b)))<=1e-9,derivative_solve=float(np.linalg.norm(S@Xd+Sc.imag/1e-25@X0-bc.imag/1e-25)/max(1.,np.linalg.norm(bc.imag/1e-25)))<=1e-9)
    return dict(state=state,X0=X0,Xd=Xd,validation={key:bool(value) for key,value in valid.items()},existing=False,metrics=dict(near=near,zero_dependence=zerror,linear_dependence=lerror,constant=constant,zero_c=cerror,work=work,moments=moment,rigid=rigid,derivative=checks.relative(Rd-Ralt,Rd,weight),solve=float(np.linalg.norm(S@X0-b)/max(1.,np.linalg.norm(b)))))
def extract(item,plane):
    pts=np.column_stack((plane,np.zeros(len(plane))));state=item['state'];zero=family.evaluate(state,0.,item['X0'],pts);pts[:,2]=1.;first=2*family.evaluate(state,1e-25j,item['X0']+1e-25j*item['Xd'],pts).imag/1e-25;J,_=family.raw(state['q'],0.,pts);return J,zero[:,:,:144],first[:,:,:144]
def assemble(k,n,amp,order,cache,map_paths,mapdir):
    xy,cells,coords,_=common.element.mesh(k,n,amp);free=np.flatnonzero(np.arange(len(xy))%(3*n+1));index={int(v):i for i,v in enumerate(free)};p,w=np.polynomial.legendre.leggauss(order);plane=np.array(list(itertools.product(p,p)));weights=np.array([a*b for a,b in itertools.product(w,w)]);rows=6*len(plane)*len(cells);A0=np.zeros((rows,6*len(free)));C0=np.zeros_like(A0);A1=np.zeros_like(A0);D1=np.zeros((rows,3*len(free)));targets={kind:np.zeros(rows) for kind in ['quadratic','quartic']};zero_c=0.;map_keys=[]
    for cell,(_,_,ids) in enumerate(cells):
        q=coords[ids].reshape(144);key=digest_bytes(q)
        if key not in cache:cache[key]=make_map(q)
        item=cache[key]
        if key not in map_paths:
            path=mapdir/(key+'.npz');np.savez_compressed(path,**{name:item['state'][name] for name in family.KEYS},X0=item['X0'],Xd=item['Xd']);map_paths[key]=dict(artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path),existing=item['existing'],validation=item['validation'],metrics=item.get('metrics',{}))
        if not all(item['validation'].values()):raise ValueError('local_map_validation_failed:'+key)
        map_keys.append(key);J,B0,B1=extract(item,plane);weight=np.sqrt(weights*np.linalg.det(J));points,_=common.element.bilinear(plane,xy[ids][[0,3,15,12]]);H=common.frame_map(k,points[:,0]);L0=H@B0;L1=np.einsum('ij,pjn->pin',common.element.CH,B1);sl=slice(6*len(plane)*cell,6*len(plane)*(cell+1));zero_c=max(zero_c,float(np.linalg.norm(B0[:,:,common.CC])))
        for local,node in enumerate(ids):
            if int(node) not in index:continue
            j=index[int(node)];A0[sl,6*j:6*j+6]=(B0[:,:,9*local:9*local+6]*weight[:,None,None]).reshape(-1,6);C0[sl,6*j:6*j+6]=(L0[:,:,9*local:9*local+6]*weight[:,None,None]).reshape(-1,6);A1[sl,6*j:6*j+6]=(L1[:,:,9*local:9*local+6]*weight[:,None,None]).reshape(-1,6);D1[sl,3*j:3*j+3]=(L1[:,:,9*local+6:9*local+9]*weight[:,None,None]).reshape(-1,3)
        for kind in targets:targets[kind][sl]=(common.target(kind,k,points[:,0])@common.element.CH.T*weight[:,None]).reshape(-1)
    M,area=common.h1tool.h1_matrix(xy,cells,coords,free,order);return dict(xy=xy,cells=cells,coords=coords,free=free,A0=A0,C0=C0,A1=A1,D1=D1,targets=targets,M=M,area=area,zero_c=zero_c,map_keys=map_keys)
def solve_semidefinite_ball(T,b,N,M,radius):
    nf=len(M);nk=N.shape[1];nq=nk+3*nf;P=np.zeros((nf,6,nq));P[:,:3,:nk]=N.reshape(nf,6,nk)[:,3:6]
    for i in range(nf):P[i,3:,nk+3*i:nk+3*i+3]=np.eye(3)
    L=np.linalg.cholesky(M);F=np.einsum('ij,jck->ick',L.T,P).reshape(6*nf,nq);U,s,V=svd(F,full_matrices=True,check_finite=False);del U;rh=[int(np.count_nonzero(s>v*s[0])) for v in [1e-9,1e-10,1e-11]]
    if rh!=[rh[1]]*3:raise ValueError('H1_metric_rank_uncertain:'+str(rh))
    r=rh[1];W=V[:r].T/s[:r];Z=V[r:].T;Ap=T@W;Af=T@Z;free_rank=0;free_s=np.array([]);Uf=np.zeros((len(b),0));Vf=np.zeros((0,Z.shape[1]))
    if Z.shape[1]:
        Uf,free_s,Vf=svd(Af,full_matrices=False,check_finite=False);rr=[int(np.count_nonzero(free_s>v*free_s[0])) for v in [1e-9,1e-10,1e-11]]
        if rr!=[rr[1]]*3:raise ValueError('unpenalized_rank_uncertain:'+str(rr))
        free_rank=rr[1];Uf=Uf[:,:free_rank];Vf=Vf[:free_rank];free_s=free_s[:free_rank]
    A=Ap-Uf@(Uf.T@Ap);e=b-Uf@(Uf.T@b);U,sv,Vh=svd(A,full_matrices=False,check_finite=False);rr=[int(np.count_nonzero(sv>v*sv[0])) for v in [1e-9,1e-10,1e-11]]
    if rr!=[rr[1]]*3:raise ValueError('curvature_rank_uncertain:'+str(rr))
    ra=rr[1];rhs=U[:,:ra].T@e;sv=sv[:ra];Vh=Vh[:ra]
    def z(lam):return Vh.T@(sv*rhs/(sv*sv+lam))
    uncon=z(0.);lam=0.
    if np.linalg.norm(uncon)>radius:
        hi=1.
        while np.linalg.norm(z(hi))>radius:hi*=10
        lam=brentq(lambda value:np.linalg.norm(z(value))/radius-1,0.,hi,xtol=1e-16,rtol=1e-13)
    zp=z(lam);unpen=Vf.T@((Uf.T@(b-Ap@zp))/free_s) if free_rank else np.zeros(Z.shape[1]);q=W@zp+Z@unpen;res=T@q-b;gf=Af.T@res;gp=Ap.T@res+lam*zp;kkt=max(float(np.linalg.norm(gf)/max(1.,np.linalg.norm(Af)*np.linalg.norm(b))),float(np.linalg.norm(gp)/max(1.,np.linalg.norm(Ap)*np.linalg.norm(b)+lam*np.linalg.norm(zp))));return q,dict(error=float(np.linalg.norm(res)/np.linalg.norm(b)),h1=float(np.linalg.norm(F@q)),unconstrained_h1=float(np.linalg.norm(uncon)),multiplier=float(lam),kkt=kkt,H1_ranks=rh,H1_nullity=Z.shape[1],unpenalized_objective_rank=free_rank,unpenalized_nonuniqueness=Z.shape[1]-free_rank,positive_objective_rank=ra,positive_metric_dimension=r)
