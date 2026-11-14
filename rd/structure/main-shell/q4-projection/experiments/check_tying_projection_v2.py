import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
from scipy.linalg import svdvals,svd,qr
import tying_mid_projection_q4_v2 as m
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))
def relative(a,b):return float(np.linalg.norm(a)/max(1e-30,np.linalg.norm(b)))
def ranks(s):return [int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]]
def generic():
    p,_=m.native_plane();x=(p[:,0]+1)/2;y=p[:,1]/2;f=.2*x*x+.08*x**4*y**4;fx=.4*x+.32*x**3*y**4;fy=.32*x**4*y**3
    return np.vstack((np.column_stack((x,y,f)),np.column_stack((-fx,-fy,np.ones(len(x)))),np.zeros((25,3)))).reshape(225)

def local(q,name,out):
    records=[];rng=np.random.default_rng(20261103);A=np.column_stack([np.einsum('ij,nj->ni',mb,q.reshape(75,3)).reshape(225) for mb in m.MB])
    rigid=[]
    for axis in np.eye(3):
        u=np.zeros((75,3));u[:25]=axis;rigid.extend([u.reshape(225),np.cross(axis,q.reshape(75,3)).reshape(225)])
    rigid=np.array(rigid).T
    _,Jc,_=m.raw(q,np.zeros((1,3)),0.);L=float(np.sqrt(np.linalg.norm(Jc[0,:,0])*np.linalg.norm(Jc[0,:,1])))
    scales=np.r_[np.full(75,L),np.ones(75),np.full(75,1/L)];v1=rng.normal(size=225)*scales;v2=rng.normal(size=225)*scales;current=q+.002*v1
    p0,w0,J0,B00,R00,*_=m.local_system(q,0.)
    pc,wc,Jc,Bc,Rc,*_=m.local_system(q,1e-25j)
    plane,_=m.quadrature(9,2);B1=m.expansion(q,plane)[2]
    expected=(B1[:,None]*p0.reshape(81,5,3)[:,:,2,None,None]/2).reshape(-1,6,225)
    derivative=rel(Rc.imag/1e-25-expected,expected)
    angle=.37;Qrigid=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
    for h in [.01,.001,.0001]:
        p,w,J,B,R,pn,wn,Jn,Bn=m.local_system(q,h,True)
        mu=w*np.linalg.det(J);mu/=mu.sum();volume=(w*np.linalg.det(J)).sum()*h/2
        W=(np.einsum('ij,pjn->pin',m.CH,R)*np.sqrt(mu)[:,None,None]).reshape(-1,225);scale=1/np.linalg.norm(W,axis=0)
        sc=svdvals(W*scale,check_finite=False);Wn=(np.einsum('ij,pjn->pin',m.CH,Bn)*np.sqrt(wn*np.linalg.det(Jn))[:,None,None]).reshape(-1,225)
        sn=svdvals(Wn/np.linalg.norm(Wn,axis=0),check_finite=False)
        boundary=m.boundary(q,h);force_c=volume*np.einsum('p,pki->ik',mu,R);force_n=h/2*np.einsum('p,pki->ik',wn*np.linalg.det(Jn),Bn)
        metrics=dict(constant=rel(R@A-np.eye(6),np.eye(6)),rigid=rel(R@rigid,R),
            conservative_work=relative(force_c-boundary,boundary),native_work=relative(force_n-boundary,boundary),
            min_det=float(min(np.min(np.linalg.det(J)),np.min(np.linalg.det(Jn)))),thin_derivative=derivative)
        energy,force,tangent=m.nonlinear(q,current,h,v1,v2,True)
        ec,_,_=m.nonlinear(q,current+1e-25j*v1,h,v1,v2,True);_,fc,_=m.nonlinear(q,current+1e-25j*v2,h,v1,v2,True)
        metrics['energy_force']=float(abs(ec.imag/1e-25-force@v1)/max(1.,abs(force@v1)))
        metrics['force_tangent']=float(abs(v1@(fc.imag/1e-25)-tangent)/max(1.,abs(tangent)))
        _,g,_=m.tables(p,h*p[:,2]/2);G=g@np.linalg.inv(J)
        Fcur=np.eye(3)+np.einsum('na,pnb->pab',(current-q).reshape(75,3),G)
        metrics['min_F']=float(np.min(np.linalg.det(Fcur)))
        qrot=np.einsum('ij,nj->ni',Qrigid,current.reshape(75,3));qrot[:25]+=L*np.array([.2,-.3,.1])
        ur,fr,tr=m.nonlinear(q,qrot.reshape(225),h,v1,v2)
        expected_force=np.einsum('ij,nj->ni',Qrigid,force.reshape(75,3)).reshape(225)
        metrics['rotation_energy']=rel(ur-energy,energy);metrics['rotation_force']=rel(fr-expected_force,expected_force)
        if name=='flat':
            D=np.zeros((225,50))
            for i in range(25):D[3*i:3*i+2,2*i:2*i+2]=np.eye(2)
            metrics['flat_membrane_identity']=rel(R@D-B@D,B@D)
        tau=np.column_stack([np.sin((i+1)*p[:,0])+.2*p[:,1]**2+.1*p[:,2] for i in range(6)])
        dual=m.dual_force(q,h,B,B.reshape(81,5,6,225)[:,2],tau)
        direct=np.einsum('p,pki,pk->i',mu,R,tau)
        metrics['general_dual_work']=rel(direct-dual,direct)
        twice=m.project_fields(q,h,R)[2]
        metrics['idempotence']=rel(twice-R,R)
        checks=dict(native_rank=ranks(sn)==[219]*3,conservative_rank=ranks(sc)==[219]*3,orientation=min(metrics['min_det'],metrics['min_F'])>0,
            conservative_fields=all(v<=1e-9 for k,v in metrics.items() if k not in ['min_det','min_F','native_work']),native_work=metrics['native_work']<=1e-9)
        path=out/f'local-{name}-h{h:g}.npz';np.savez_compressed(path,q=q,current=current,v1=v1,v2=v2,energy=energy,force=force,tangent=tangent)
        records.append(dict(geometry=name,h=h,metrics=metrics,native_ranks=ranks(sn),conservative_ranks=ranks(sc),checks={k:bool(v) for k,v in checks.items()},
            artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path)))
    return records

def shared(mesh,k,amp,method,out,n=2):
    a=m.assemble(mesh,k,9,method);b=m.assemble(mesh,k,11,method);scale=1/np.linalg.norm(a['A0'],axis=0)
    _,s9,V9=svd(a['A0']*scale,full_matrices=False,check_finite=False);_,s11,V11=svd(b['A0']*scale,full_matrices=False,check_finite=False);r=ranks(s9)[1]
    V=V9[r:].T;Vo=V11[ranks(s11)[1]:].T
    controls=np.column_stack([m.displacement(mesh,kind,k)[a['free_rd']].reshape(-1)/scale for kind in ['quadratic','cubic']])
    Z,_=qr(controls,mode='economic',check_finite=False);Q,_=qr(V.T@Z,mode='full',check_finite=False);Vn=np.column_stack((Z,V@Q[:,2:]));N=scale[:,None]*Vn
    metrics=dict(anchor=float(np.linalg.norm(Vn@Vn.T-V@V.T)),orthogonal=float(np.linalg.norm(Vn.T@Vn-np.eye(Vn.shape[1]))),
        projector=float(np.linalg.norm(Vn@Vn.T-Vo@Vo.T)),kernel=float(np.linalg.norm((a['A0']*scale)@Vn)/(np.linalg.norm(a['A0']*scale)*np.linalg.norm(Vn))))
    checks=dict(rank=ranks(s9)==ranks(s11)==[r]*3,anchor=metrics['anchor']<=1e-7 and metrics['orthogonal']<=1e-10,
        projector=metrics['projector']<=1e-7,kernel=metrics['kernel']<=1e-10)
    solutions={};stored=dict(ref=mesh['ref'],cells=mesh['cells'],free_rd=a['free_rd'],free_c=a['free_c'],kernel=N,scale=scale,
        xy0=mesh['xy'][0],xy1=mesh['xy'][1],xy2=mesh['xy'][2],full_s9=s9,full_s11=s11)
    error=None
    if all(checks.values()):
        try:
            T=np.column_stack((a['A1']@N,a['D1']));F=m.metric_factor(a,N,mesh)
            for kind in ['quadratic','cubic','quartic','exponential']:
                rh=m.reference_h1(kind,k,40);q,fit=m.solve_ball(T,a['targets'][kind],F,2*rh);u=np.zeros_like(mesh['ref'])
                u[a['free_rd']]=(N@q[:N.shape[1]]).reshape(-1,3);u[a['free_c']]=q[N.shape[1]:].reshape(-1,3)
                rd=u[a['free_rd']].reshape(-1);cc=u[a['free_c']].reshape(-1)
                e=float(np.linalg.norm(b['A1']@rd+b['D1']@cc-b['targets'][kind])/np.linalg.norm(b['targets'][kind]))
                dc=u[np.r_[mesh['free'][1],mesh['free'][2]]];h1=float(np.sqrt(np.sum(dc*(b['M']@dc))))
                checks[kind]=fit['kkt']<=1e-8 and abs(e-fit['error'])<=1e-8 and abs(h1/fit['h1']-1)<=1e-8 and h1<=2*rh*(1+1e-8) and (kind in ['quartic','exponential'] or e<=1e-8)
                solutions[kind]=dict(error=e,h1=h1,fit=fit);stored['u_'+kind]=u
        except Exception as exc:error=repr(exc);checks['first_order_solve']=False
    path=out/f'{method}-k{k:g}-a{amp:g}-n{n}.npz';np.savez_compressed(path,**stored)
    return dict(method=method,kappa=k,amplitude=amp,n=n,nullity=N.shape[1],ranks9=ranks(s9),ranks11=ranks(s11),metrics=metrics,solutions=solutions,
        checks={k:bool(v) for k,v in checks.items()},error=error,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))

def main():
    dest=ROOT/'rd/structure/main-shell/q4-projection/records/tying-projection-v2-result.json';assert not dest.exists();start=time.monotonic();out=ROOT/'build/analysis'/('tying-projection-v2-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    locals=[];records=[];stop=None
    try:
        for name,q in [('flat',m.mesh(0.,1,0.)['ref'][m.mesh(0.,1,0.)['cells'][0]].reshape(225)),
            ('curved-distorted',m.mesh(.4,2,.2)['ref'][m.mesh(.4,2,.2)['cells'][0]].reshape(225)),('quartic-reference',generic())]:
            locals.extend(local(q,name,out));print('local',name,[(r['h'],r['metrics']['native_work'],r['metrics']['conservative_work']) for r in locals if r['geometry']==name],flush=True)
        for k,amp in itertools.product([0.,.4],[0.,.2]):
            mesh=m.mesh(k,2,amp)
            for method in ['raw','conservative']:
                tick=time.monotonic();r=shared(mesh,k,amp,method,out)
                r['checks']['resources']=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3
                records.append(r);print('shared',method,k,amp,r['nullity'],{k:v['error'] for k,v in r['solutions'].items()},r['checks'],r['error'],flush=True)
                if not r['checks']['resources']:raise RuntimeError('shared_resource_gate_failed')
        # One prescribed refinement, only if the conservative operator has all
        # necessary algebraic gates and has improved both curved comparisons.
        selected=[r for r in records if r['method']=='conservative'];raws=[r for r in records if r['method']=='raw']
        local_ok=all(all(v for k,v in r['checks'].items() if k!='native_work') for r in locals)
        improved=all(s['solutions'].get('quartic',{}).get('error',1)>=0 and s['solutions'].get('quartic',{}).get('error',1)<.99*next(r for r in raws if (r['kappa'],r['amplitude'])==(s['kappa'],s['amplitude']))['solutions']['quartic']['error'] for s in selected if s['kappa'])
        if local_ok and all(all(r['checks'].values()) for r in selected) and improved:
            tick=time.monotonic();r=shared(m.mesh(.4,4,.2),.4,.2,'conservative',out,4)
            r['checks']['resources']=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3
            records.append(r);print('refined',r['nullity'],{k:v['error'] for k,v in r['solutions'].items()},r['checks'],flush=True)
    except Exception as exc:stop=repr(exc)
    conservative_local=bool(locals) and all(all(v for k,v in r['checks'].items() if k!='native_work') for r in locals)
    native_work=bool(locals) and all(r['checks']['native_work'] for r in locals)
    sources=[Path(__file__),Path(m.__file__),ROOT/"rd/structure/main-shell/q4-projection/experiments/covariant_projection_basis.py",ROOT/'rd/structure/main-shell/q4-projection/tying-projection-v2-formulation.md',ROOT/'rd/structure/main-shell/q4-projection/records/tying-projection-v2-protocol.md',ROOT/'rd/structure/main-shell/q4-projection/records/tying-projection-v2-preparation.json']
    dest.write_text(json.dumps(dict(local_records=locals,records=records,conservative_local_pass=conservative_local,native_work_pass=native_work,
        all_shared_numerical_checks=len(records)>=8 and all(all(r['checks'].values()) for r in records) and stop is None,stop_reason=stop,physical_adoption=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
