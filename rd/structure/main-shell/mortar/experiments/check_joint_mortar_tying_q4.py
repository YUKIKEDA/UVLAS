"""Pre-registered necessary gates and adoption disposition for mortar EL2."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time,itertools
import numpy as np
from scipy.linalg import svd,qr
import joint_mortar_tying_q4 as f
import check_q9_plate_convergence as memory
m=f.m;ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))

def generic():
    p,_=m.native_plane();x=(p[:,0]+1)/2;y=p[:,1]/2
    return np.vstack((np.column_stack((x,y,.2*x*x+.08*x**4*y**4)),np.column_stack((-.4*x-.32*x**3*y**4,-.32*x**4*y**3,np.ones(25))),np.zeros((25,3)))).reshape(225)

def local(q,name,h,out):
    plane,wp=m.quadrature(9,2);z,wz=np.polynomial.legendre.leggauss(5)
    p=np.array([[s,t,c] for (s,t),c in itertools.product(plane,z)]);w=np.array([a*b for a,b in itertools.product(wp,wz)])
    _,J,_=m.raw(q,p,h*p[:,2]/2);B=f.linearized(q,h,p)
    mu=w*np.linalg.det(J)*h/2;W=(np.einsum('ij,pjn->pin',m.CH,B)*np.sqrt(mu)[:,None,None]).reshape(-1,225)
    singular=svd(W/np.linalg.norm(W,axis=0),compute_uv=False,check_finite=False)
    affine=np.column_stack([np.einsum('ij,nj->ni',b,q.reshape(75,3)).reshape(-1) for b in m.MB]);rigid=[]
    for axis in np.eye(3):
        t=np.zeros((75,3));t[:25]=axis;rigid.extend([t.reshape(-1),np.cross(axis,q.reshape(75,3)).reshape(-1)])
    rigid=np.array(rigid).T;inside=np.einsum('p,pki->ik',mu,B);outside=m.boundary(q,h,independent=True)
    rng=np.random.default_rng(20261104);current=q+.001*rng.normal(size=225);v1=rng.normal(size=225);v2=rng.normal(size=225)
    U,force,tangent=f.raw_nonlinear(q,current,h,v1,v2)
    Uc,_,_=f.raw_nonlinear(q,current+1e-25j*v1,h,v1,v2);_,fc,_=f.raw_nonlinear(q,current+1e-25j*v2,h,v1,v2)
    axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis);A=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    Q=np.eye(3)+np.sin(.71)*A+(1-np.cos(.71))*A@A
    rotated=current.reshape(75,3)@Q.T;rotated[:25]+=np.array([.2,-.1,.3])
    Ur,fr,_=f.raw_nonlinear(q,rotated.reshape(-1),h,v1,v2)
    _,g,_=m.tables(p,h*p[:,2]/2);grad=g@np.linalg.inv(J);F=np.eye(3)+np.einsum('na,pnb->pab',(current-q).reshape(75,3),grad)
    metrics=dict(affine=rel(B@affine-np.eye(6),np.broadcast_to(np.eye(6),(len(B),6,6))),
        rigid=rel(B@rigid,B),constant_stress_work=rel(inside-outside,outside),
        energy_force=float(abs(Uc.imag/1e-25-force@v1)/max(1e-30,abs(force@v1))),
        force_tangent=float(abs(v1@(fc.imag/1e-25)-tangent)/max(1e-30,abs(tangent))),rotation_energy=rel(Ur-U,U),
        rotation_force=rel(fr-(force.reshape(75,3)@Q.T).reshape(-1),force),min_J=float(np.min(np.linalg.det(J))),min_F=float(np.min(np.linalg.det(F))))
    if name=='flat':
        u=np.zeros((75,3));u[:25,0]=q.reshape(75,3)[:25,0]**3;u[:25,1]=q.reshape(75,3)[:25,1]**3
        N,_=m.shape(p,0);xy=N@q.reshape(75,3)[:25];expected=np.zeros((len(p),6));expected[:,0]=3*xy[:,0]**2;expected[:,1]=3*xy[:,1]**2
        metrics['nonuniform_membrane_patch']=rel(B@u.reshape(-1)-expected,expected)
    checks=dict(rank=f.ranks(singular)==[219]*3,orientation=min(metrics['min_J'],metrics['min_F'])>0,
        identities=all(v<=1e-9 for k,v in metrics.items() if k not in ['min_J','min_F']))
    path=out/f'local-{name}-{h:g}.npz';np.savez_compressed(path,q=q,current=current,v1=v1,v2=v2,energy=U,force=force,tangent=tangent)
    return dict(geometry=name,h=h,ranks=f.ranks(singular),metrics=metrics,checks={k:bool(v) for k,v in checks.items()},artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))

def joint(mesh,k,amp,degree,out,n=2):
    a,Zrd,Zc,meta=f.shared(mesh,k,degree,9);b,_,_,metab=f.shared(mesh,k,degree,11,True)
    reduced=a['A0']@Zrd;scale=1/np.linalg.norm(reduced,axis=0)
    _,s,Vh=svd(reduced*scale,full_matrices=False,check_finite=False);_,s11,Vh11=svd((b['A0']@Zrd)*scale,full_matrices=False,check_finite=False)
    rank=f.ranks(s)[1];V=Vh[rank:].T;Vo=Vh11[f.ranks(s11)[1]:].T
    controls=np.column_stack([Zrd.T@m.displacement(mesh,kind,k)[a['free_rd']].reshape(-1)/scale for kind in ['quadratic','cubic']])
    A,_=qr(controls,mode='economic');Q,_=qr(V.T@A,mode='full');Vn=np.column_stack((A,V@Q[:,2:]));N=Zrd@(scale[:,None]*Vn)
    G11,_=f.constraints(mesh,0.,degree,11,7,True);G9,_=f.constraints(mesh,0.,degree)
    metrics=dict(anchor=float(np.linalg.norm(Vn@Vn.T-V@V.T)),orthogonal=float(np.linalg.norm(Vn.T@Vn-np.eye(Vn.shape[1]))),
        projector=float(np.linalg.norm(Vn@Vn.T-Vo@Vo.T)),kernel=rel(a['A0']@N,a['A0'])/float(np.linalg.norm(N)),
        rd_constraint=rel(np.kron(G11[:,a['free_rd']],np.eye(3))@N,G11)/max(1.,np.linalg.norm(N)),
        c_constraint=rel(np.kron(G11[:,a['free_c']],np.eye(3))@Zc,G11)/max(1.,np.linalg.norm(Zc)))
    checks=dict(rank=f.ranks(s)==f.ranks(s11)==[rank]*3,anchor=metrics['anchor']<=1e-7 and metrics['orthogonal']<=1e-10,
        projector=metrics['projector']<=1e-7,kernel=metrics['kernel']<=1e-10,constraints=max(metrics['rd_constraint'],metrics['c_constraint'])<=1e-10)
    solutions={};stored=dict(ref=mesh['ref'],cells=mesh['cells'],free_rd=a['free_rd'],free_c=a['free_c'],kernel=N,c_basis=Zc,
        scale=scale,Zrd=Zrd,singular=s,singular11=s11,xy0=mesh['xy'][0],xy1=mesh['xy'][1],xy2=mesh['xy'][2])
    stop=None
    if all(checks.values()):
        try:
            T=np.column_stack((a['A1']@N,a['D1']@Zc));nd=3*len(mesh['free'][1]);nc=3*len(mesh['free'][2]);nr=3*len(mesh['free'][0]);nk=N.shape[1]
            P=np.zeros((nd+nc,T.shape[1]));P[:nd,:nk]=N[nr:];P[nd:,nk:]=Zc
            F=np.linalg.cholesky(np.kron(a['M'],np.eye(3))).T@P
            for kind in ['quadratic','cubic','quartic','exponential']:
                refH=m.reference_h1(kind,k,40);q,fit=m.solve_ball(T,a['targets'][kind],F,2*refH)
                u=np.zeros_like(mesh['ref']);u[a['free_rd']]=(N@q[:nk]).reshape(-1,3);u[a['free_c']]=(Zc@q[nk:]).reshape(-1,3)
                e=float(np.linalg.norm(b['A1']@u[a['free_rd']].reshape(-1)+b['D1']@u[a['free_c']].reshape(-1)-b['targets'][kind])/np.linalg.norm(b['targets'][kind]))
                dc=u[np.r_[mesh['free'][1],mesh['free'][2]]];h1=float(np.sqrt(np.sum(dc*(b['M']@dc))))
                checks[kind]=fit['kkt']<=1e-8 and abs(e-fit['error'])<=1e-8 and abs(h1/fit['h1']-1)<=1e-8 and h1<=2*refH*(1+1e-8) and (kind in ['quartic','exponential'] or e<=1e-8)
                solutions[kind]=dict(error=e,h1=h1,fit=fit);stored['u_'+kind]=u
        except Exception as exc:stop=repr(exc);checks['solve']=False
    path=out/f'joint-p{degree}-k{k:g}-a{amp:g}-n{n}.npz';np.savez_compressed(path,**stored)
    return dict(kappa=k,amplitude=amp,n=n,degree=degree,nullity=N.shape[1],ranks=f.ranks(s),ranks11=f.ranks(s11),constraint_metadata=meta,
        metrics=metrics,solutions=solutions,checks={k:bool(v) for k,v in checks.items()},error=stop,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))

def interface_checks(mesh,k,amp,degree):
    rng=np.random.default_rng(701);scans=[f.rigid_scan(mesh,h,degree) for h in [.01,.0001]]
    h=.01;G,meta=f.constraints(mesh,h,degree);Gi,_=f.constraints(mesh,h,degree,11,7,True)
    u=rng.normal(size=mesh['ref'].shape);lam=rng.normal(size=(len(G),3));reaction=G.T@lam
    work=float(np.sum(lam*(Gi@u)));direct=float(np.sum(reaction*u));Q,_=np.linalg.qr(rng.normal(size=(3,3)))
    qr=mesh['ref']+u;qrrot=qr@Q.T+np.array([.3,-.1,.2]);eq=G@qr;rot=Gi@qrrot
    metrics=dict(integral=rel(G-Gi,G),dual=abs(direct-work)/max(1e-30,abs(direct)),rotation=rel(rot-eq@Q.T,eq),
        reaction_rotation=rel(Gi.T@(lam@Q.T)-reaction@Q.T,reaction))
    patches=[f.patch(mesh,h,degree,np.eye(3)+.03*b) for b in m.MB]
    checks=dict(interface=all(v<=1e-9 for v in metrics.values()),patch=all(r['residual']<=1e-9 and r['gap']<=1e-9 for r in patches),rigid=all(r['passed'] for r in scans))
    return dict(kappa=k,amplitude=amp,degree=degree,scans=scans,metrics=metrics,patches=patches,checks={k:bool(v) for k,v in checks.items()})

def main():
    dest=ROOT/'rd/structure/main-shell/mortar/records/joint-mortar-tying-q4-result.json';assert not dest.exists();start=time.monotonic()
    out=ROOT/'build/analysis'/('joint-mortar-tying-q4-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();locals=[];interfaces=[];records=[];stop=None
    try:
        for name,q in [('flat',m.mesh(0.,1,0.)['ref'][m.mesh(0.,1,0.)['cells'][0]].reshape(-1)),('curved-distorted',m.mesh(.4,2,.2)['ref'][m.mesh(.4,2,.2)['cells'][0]].reshape(-1)),('quartic-reference',generic())]:
            for h in [.01,.001,.0001]:locals.append(local(q,name,h,out))
            print('local',name,[r['checks'] for r in locals if r['geometry']==name],flush=True)
        for degree,k,amp in itertools.product([1],[0.,.4],[0.,.2]):
            tick=time.monotonic();r=interface_checks(f.broken_mesh(k,2,amp),k,amp,degree)
            r['checks']['resources']=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3
            interfaces.append(r);print('interface',degree,k,amp,[(s['spurious_free_modes'],s['spurious_supported_modes']) for s in r['scans']],r['checks'],flush=True)
        for degree in [1]:
            selected=[r for r in interfaces if r['degree']==degree]
            if not all(all(r['checks'].values()) for r in selected):continue
            for k,amp in itertools.product([0.,.4],[0.,.2]):
                tick=time.monotonic();r=joint(f.broken_mesh(k,2,amp),k,amp,degree,out)
                r['checks']['resources']=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3
                records.append(r);print('joint',degree,k,amp,r['nullity'],{k:v['error'] for k,v in r['solutions'].items()},r['checks'],r['error'],flush=True)
        baseline=json.loads((ROOT/'rd/structure/main-shell/q4-projection/records/tying-projection-v2-result.json').read_text())['records'];dispositions=[]
        for degree in [1]:
            selected=[r for r in records if r['degree']==degree];algebra=all(all(r['checks'].values()) for r in locals) and all(all(r['checks'].values()) for r in interfaces if r['degree']==degree) and len(selected)==4 and all(all(r['checks'].values()) for r in selected)
            improved=algebra and all(r['solutions']['quartic']['error']<.99*next(b for b in baseline if b['method']=='raw' and b['kappa']==r['kappa'] and b['amplitude']==r['amplitude'])['solutions']['quartic']['error'] for r in selected if r['kappa'])
            dispositions.append(dict(degree=degree,algebra_pass=bool(algebra),both_curved_improved=bool(improved),main_shell_adopted=False,
                status='refinement_required' if improved else 'rejected_necessary_gate'))
            if improved:
                tick=time.monotonic();r=joint(f.broken_mesh(.4,4,.2),.4,.2,degree,out,4)
                r['checks']['resources']=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3;records.append(r)
                print('refined',degree,r['nullity'],r['checks'],flush=True)
    except Exception as exc:stop=repr(exc);dispositions=[]
    sources=[Path(__file__),Path(f.__file__),Path(m.__file__),ROOT/'rd/structure/main-shell/mortar/joint-mortar-tying-q4-formulation.md',ROOT/'rd/structure/main-shell/mortar/records/joint-mortar-tying-q4-protocol.md',ROOT/'rd/structure/main-shell/mortar/records/joint-mortar-tying-q4-preparation.json',ROOT/'rd/structure/main-shell/q4-projection/records/tying-projection-v2-result.json']
    dest.write_text(json.dumps(dict(local_records=locals,interfaces=interfaces,records=records,dispositions=dispositions,stop_reason=stop,
        physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
