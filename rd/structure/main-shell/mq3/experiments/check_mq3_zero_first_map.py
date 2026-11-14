"""Fixed-basis zero/first-thickness KKT map checks, with signed/complex h."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from scipy.linalg import qr,solve,svdvals
import mq3_thickness_family as family
import mq3_joint_moments as original
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def cnorm(A,weight):
    value=np.einsum('ij,pjn->pin',family.base.CH,A);return float(np.linalg.norm(value*np.sqrt(weight)[:,None,None]))
def relative(A,B,weight):return cnorm(A,weight)/max(1.,cnorm(B,weight))
def ranks_at(state,h):
    _,d,w=family.volume(state['q'],h);A=(np.einsum('ij,pjn->pin',family.base.CH,d['W'])*np.sqrt(w)[:,None,None]).reshape(-1,117);norm=np.linalg.norm(A,axis=0);active=norm>1e-14*norm.max();Q,R,p=qr(A[:,active]/norm[active],mode='economic',pivoting=True,check_finite=False);br=family.base.ranks(np.diag(R));Qf,Rf=qr(A@state['T'],mode='economic',check_finite=False);span=float(np.linalg.norm(A-Qf@(Qf.T@A))/np.linalg.norm(A));F,g=family.moments(state['q'],h,d,w);D=F@state['T'];rn=np.linalg.norm(D,axis=1);mask=rn>1e-13*rn.max();_,Rc,_=qr((D[mask]/rn[mask,None]).T,mode='economic',pivoting=True,check_finite=False);cr=family.base.ranks(np.diag(Rc));return br,cr,span
def main():
    dest=ROOT/'rd/structure/main-shell/mq3/records/mq3-zero-first-map-result.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();parent=ROOT/'rd/structure/main-shell/mq3/records/mq3-local-feasibility-result.json';data=json.loads(parent.read_text());assert data['all_checks_pass'];records=[];stop=None;out=ROOT/'build/analysis'/('mq3-zero-first-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();groups={}
    for record in data['records']:groups.setdefault((record['kappa'],record['n'],record['cell']),{})[record['h']]=record
    for source,sha in data['sources'].items():assert digest(ROOT/source)==sha
    points,w=family.base.quadrature(11);zp,zw=np.polynomial.legendre.leggauss(11);eps=1e-25
    try:
        for (k,n,cell),saved_records in groups.items():
            first=saved_records[.01];path=ROOT/first['artifact'];assert digest(path)==first['artifact_sha256'];q=np.load(path)['q'];state=family.choose(q);X0,S0,b0,F0,g0=family.coefficients(state,0.);Sc,bc,Fc,gc=family.system(state,1j*eps);Sd=Sc.imag/eps;bd=bc.imag/eps;Xd=solve(S0,bd-Sd@X0,assume_a='sym',check_finite=False);Xc,_,_,_,_=family.coefficients(state,1j*eps);R0=family.evaluate(state,0.,X0,points);Rd=family.evaluate(state,1j*eps,X0+1j*eps*Xd,points).imag/eps;Rcd=family.evaluate(state,1j*eps,Xc,points).imag/eps;d0=family.data(q,0.,points);dc=family.data(q,1j*eps,points);weight=w*np.linalg.det(d0['J']);weight/=weight.sum();wc=w*np.linalg.det(dc['J']);wc/=wc.sum();wd=wc.imag/eps;Ed=dc['E'].imag/eps
            derivative_difference=relative(Rd-Rcd,Rd,weight);Rshape=R0.reshape(11,11,11,6,150);Dshape=Rd.reshape(11,11,11,6,150);mean=np.einsum('z,abzkn->abkn',zw/2,Rshape);B1=2*np.einsum('z,z,abzkn->abkn',zw,zp,Dshape)/np.sum(zw*zp*zp);zero_dependence=relative((Rshape-mean[:,:,None]).reshape(-1,6,150),R0,weight);linear_dependence=relative((Dshape-zp[None,None,:,None,None]*B1[:,:,None]/2).reshape(-1,6,150),Rd,weight);constant0=float(np.max(abs(R0[:,:,144:]-np.eye(6))));constant1=float(np.max(abs(Rd[:,:,144:])));c_indices=np.array([9*a+b for a in range(16) for b in range(6,9)]);zero_c=float(np.max(abs(R0[:,:,c_indices])));rigid=0.
            for axis in np.eye(3):
                v=np.zeros((16,3,3));v[:,0]=axis
                for mode in [v.reshape(144),np.cross(np.broadcast_to(axis,(16,3,3)),q.reshape(16,3,3)).reshape(144)]:
                    for operator in [R0[:,:,:144],Rd[:,:,:144]]:rigid=max(rigid,float(np.linalg.norm(operator@mode)/(max(1.,np.linalg.norm(operator))*np.linalg.norm(mode))))
            delta=R0-d0['E'];work0=np.einsum('p,pkn->kn',weight,delta);work1=np.einsum('p,pkn->kn',weight,Rd-Ed)+np.einsum('p,pkn->kn',wd,delta);work_error=max(float(np.linalg.norm(work0)),float(np.linalg.norm(work1)))/max(1.,float(np.linalg.norm(np.einsum('p,pkn->kn',weight,d0['E']))));c0=state['T']@X0[:state['T'].shape[1]];cd=state['T']@Xd[:state['T'].shape[1]];moment0=float(np.linalg.norm(F0@c0-g0)/max(1.,np.linalg.norm(g0)));moment1=float(np.linalg.norm(Fc.imag/eps@c0+F0@cd-gc.imag/eps)/max(1.,np.linalg.norm(gc.imag/eps)));linear_residual=float(np.linalg.norm(S0@X0-b0)/np.linalg.norm(b0));derivative_residual=float(np.linalg.norm(S0@Xd+Sd@X0-bd)/max(1.,np.linalg.norm(bd)));s=svdvals(S0,check_finite=False);rcond=float(s[-1]/s[0]);near=[];original_difference=0.;cache={}
            for h in [0.,-.01,.01,-.0001,.0001]:
                br,cr,span=ranks_at(state,h);X,Sh,bh,Fh,gh=family.coefficients(state,h);res=float(np.linalg.norm(Fh@(state['T']@X[:state['T'].shape[1]])-gh)/max(1.,np.linalg.norm(gh)));near.append(dict(h=h,basis_ranks=br,constraint_ranks=cr,span_residual=span,constraint_residual=res));cache[h]=(X,family.evaluate(state,h,X,points))
                if h in saved_records:
                    savedpath=ROOT/saved_records[h]['artifact'];assert digest(savedpath)==saved_records[h]['artifact_sha256'];saved=np.load(savedpath);old={key:saved[key] for key in saved.files};oldR=original.evaluate(old,points);original_difference=max(original_difference,relative(cache[h][1]-oldR,oldR,weight))
            fd=[]
            for h in [.01,.005,.0025]:
                fields=[]
                for sign in [1,-1]:
                    if sign*h not in cache:
                        X,_,_,_,_=family.coefficients(state,sign*h);cache[sign*h]=(X,family.evaluate(state,sign*h,X,points))
                    fields.append(cache[sign*h][1])
                fd.append(relative((fields[0]-fields[1])/(2*h)-Rd,Rd,weight))
            convergence=(max(fd)<=1e-8) or (fd[-1]<fd[0] and all(2<=fd[i]/fd[i+1]<=6 for i in range(2)));rankcheck=all(v['basis_ranks']==state['space_ranks'] and v['constraint_ranks']==state['constraint_ranks'] for v in near) and state['space_ranks']==[state['space_ranks'][1]]*3 and state['constraint_ranks']==[state['constraint_ranks'][1]]*3;checks=dict(rank=rankcheck,span=max(v['span_residual'] for v in near)<=1e-10,full_constraints=max(v['constraint_residual'] for v in near)<=1e-10,condition=rcond>=1e-10,solve=max(linear_residual,derivative_residual)<=1e-10,derivative=derivative_difference<=1e-8,saved_finite_h=original_difference<=1e-8,zero_dependence=zero_dependence<=1e-9,linear_dependence=linear_dependence<=1e-9,constant=max(constant0,constant1,zero_c)<=1e-9,rigid=rigid<=1e-9,work=work_error<=1e-9,moments=max(moment0,moment1)<=1e-9,finite_difference=convergence,resources=time.monotonic()-start<=300 and memory.memory_bytes()<=2*1024**3);artifact=out/f'k{k:g}-n{n}-c{cell}.npz';np.savez_compressed(artifact,q=q,T=state['T'],rows=state['rows'],U=state['U'],X0=X0,Xd=Xd,S0=S0,Sd=Sd,b0=b0,bd=bd,B0=mean[:,:,:144],B1=B1[:,:,:144]);records.append(dict(kappa=k,n=n,cell=cell,basis_ranks=state['space_ranks'],constraint_ranks=state['constraint_ranks'],near=near,KKT_rcond=rcond,solve_residual=linear_residual,derivative_residual=derivative_residual,derivative_difference=derivative_difference,saved_difference=original_difference,zero_dependence=zero_dependence,linear_dependence=linear_dependence,constant0=constant0,constant1=constant1,zero_c=zero_c,rigid=rigid,work_error=work_error,moment0=moment0,moment1=moment1,finite_difference_errors=fd,checks={key:bool(value) for key,value in checks.items()},artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact)));print(k,n,cell,'checks',checks,'thickness',zero_dependence,linear_dependence,'fd',fd,flush=True)
            if not all(checks.values()):stop='required_zero_first_map_gate_failed';break
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(family.__file__),Path(family.base.__file__),Path(family.base.geometry.__file__),Path(original.__file__),parent,ROOT/'rd/structure/main-shell/mq3/mq3-zero-thickness-linearization.md',ROOT/'rd/structure/main-shell/mq3/records/mq3-zero-first-map-protocol.md']
    dest.write_text(json.dumps(dict(records=records,all_checks_pass=len(records)==10 and stop is None,stop_reason=stop,physical_adoption=False,shared_thin_kernel_evaluated=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
