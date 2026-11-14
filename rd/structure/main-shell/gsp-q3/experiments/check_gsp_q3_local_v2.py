"""Fixed 30-case gate for geometry-based stress projection."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import qr,solve,svdvals,eigvalsh
import gsp_q3_v2 as candidate
import q3_thin_expansion as geometry
import check_q9_plate_convergence as memory
from check_mq3i1_local_v2 import patch_vectors,exact_patches,reparameterize,mapped_points

ROOT=rd_paths.ROOT
KEYS=['q','frame_scale','length','center','T','active','column_norms','coeff','denominator']


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def relative(A,B):return float(np.linalg.norm(A)/max(1.,np.linalg.norm(B)))


def affine_fields(state,h,points):
    N,_,_=candidate.kinematics(state['q'],h,points)
    q=state['q'].reshape(16,3,3);z=h*points[:,2]/2
    X=N@q[:,0]+z[:,None]*(N@q[:,1])+z[:,None]**2/2*(N@q[:,2])
    phi=np.column_stack((np.ones(len(points)),(X-state['center'])/state['length']))
    return (phi[:,None,:,None]*np.eye(6)[None,:,None,:]).reshape(len(points),6,24)


def basis_gate(state):
    reports=[]
    for order in [9,11]:
        p,w=candidate.BASE.quadrature(order);weight,_=candidate.measure(state,0.,p,w)
        raw=candidate.stress_generator(state,p)
        A=(np.einsum('ij,pjn->pin',candidate.CIH,raw)*np.sqrt(weight)[:,None,None]).reshape(-1,263)
        active=state['active'];scaled=A[:,active]/state['column_norms'][active]
        rank=candidate.ranks(svdvals(scaled,check_finite=False))
        S=candidate.stress_field(state,p)
        AT=(np.einsum('ij,pjn->pin',candidate.CIH,S)*np.sqrt(weight)[:,None,None]).reshape(-1,S.shape[2]);Q,_=qr(AT,mode='economic',check_finite=False)
        reports.append(dict(order=order,ranks=rank,span=relative(A-Q@(Q.T@A),A),
            gram_identity=relative(AT.T@AT-np.eye(AT.shape[1]),np.eye(AT.shape[1]))))
    passed=all(v['ranks']==state['stress_ranks']==[state['T'].shape[1]]*3 and max(v['span'],v['gram_identity'])<=1e-10 for v in reports)
    return reports,passed


def differential(state):
    p,w=candidate.BASE.quadrature(9);_,E0=candidate.linear_raw(state['q'],0.,p)
    P0=candidate.project(state,0.,E0)
    _,Ec=candidate.linear_raw(state['q'],1e-25j,p)
    Pc=candidate.project(state,1e-25j,Ec)
    Hdot=Pc['H'].imag/1e-25;ldot=Pc['rhs'].imag/1e-25
    ad=solve(P0['H'],ldot-Hdot@P0['a'],assume_a='sym',check_finite=False)
    rd=np.einsum('ij,pjn->pin',candidate.CI,P0['S']@ad)
    alt=Pc['R'].imag/1e-25
    R0=P0['R'];shape0=R0.reshape(9,9,9,6,150);shape1=rd.reshape(9,9,9,6,150)
    z,zw=np.polynomial.legendre.leggauss(9)
    mean=np.einsum('z,abzkn->abkn',zw/2,shape0)
    slope=np.einsum('z,z,abzkn->abkn',zw,z,shape1)/np.sum(zw*z*z)
    metrics=dict(derivative=relative(rd-alt,rd),zero_dependence=relative(shape0-mean[:,:,None],shape0),
        first_linearity=relative(shape1-z[None,None,:,None,None]*slope[:,:,None],shape1),
        zero_c=float(np.max(abs(R0[:,:,geometry.CC]))),constant_derivative=float(np.max(abs(rd[:,:,144:]))))
    for h in [-.0001,.0001]:
        pp,ww,V,S,H=candidate.system(state,h)
        vals=svdvals(H,check_finite=False)
        metrics['near_rcond_violation_'+str(h)]=float(vals[-1]/vals[0]<1e-10)
    return P0['a'],ad,metrics


def main():
    dest=ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-local-v2-result.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();out=ROOT/'build/analysis'/('gsp-q3-local-v2-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    records=[];states=[];stop=None;failure=None
    try:
        for k in [0.,.4]:
            for n,amp in [(1,0.),(2,.2)]:
                xy,cells,coords,_=geometry.element.mesh(k,n,amp)
                for cell,(_,_,ids) in enumerate(cells):
                    q=coords[ids].reshape(144);state,prep=candidate.prepare(q)
                    gates,passed=basis_gate(state)
                    path=out/f'k{k:g}-n{n}-c{cell}-state.npz'
                    a0,ad,diff=differential(state)
                    np.savez_compressed(path,**{key:state[key] for key in KEYS},a0=a0,ad=ad,svd_singular=state['svd_singular'],qr_diagonal=state['qr_diagonal'])
                    states.append(dict(kappa=k,n=n,amplitude=amp,cell=cell,preparation=prep,basis_gates=gates,
                        differential=diff,passed=bool(passed and max(diff.values())<=1e-9),
                        artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path)))
                    if not states[-1]['passed']:raise RuntimeError('stress_basis_or_thickness_gate_failed')
                    alternate={}
                    for mode in ['swap','s','t']:
                        qt,P=reparameterize(q,mode);other,report=candidate.prepare(qt)
                        alternate[mode]=(other,P)
                    for h in [.01,.001,.0001]:
                        tick=time.monotonic();p,w=candidate.BASE.quadrature(9)
                        J,E=candidate.linear_raw(q,h,p);P=candidate.project(state,h,E)
                        R=P['R'];weight=P['weight'];H=P['H'];Hsv=svdvals(H,check_finite=False)
                        metrics=dict(rcond=float(Hsv[-1]/Hsv[0]),solve=relative(H@P['a']-P['rhs'],P['rhs']),
                            constant=float(np.max(abs(R[:,:,144:]-np.eye(6)))),
                            work=relative(np.einsum('p,pkn->kn',weight,R-E),np.einsum('p,pkn->kn',weight,E)))
                        _,legacy=candidate.reference.raw(q,h,p)
                        metrics['raw_agreement']=relative(E-legacy,legacy)
                        fields=affine_fields(state,h,p);FP=candidate.project(state,h,fields)
                        metrics['affine_strain']=relative(FP['R']-fields,fields)
                        metrics['affine_stress_adjoint']=relative(candidate.adjoint(P,fields)-fields,fields)
                        twice=candidate.project(state,h,R)
                        metrics['idempotence']=relative(twice['R']-R,R)
                        metrics['energy_adjoint']=relative(candidate.adjoint(P,P['sigma'])-P['sigma'],P['sigma'])
                        tau=np.concatenate((np.sin(p.sum(axis=1))[:,None,None]*np.eye(6),
                            (p[:,0]**3+p[:,1]**3+p[:,2]**3)[:,None,None]*np.eye(6)),axis=2)
                        lhs=np.einsum('p,pkj,pki->ji',weight,tau,R)
                        rhs=np.einsum('p,pkj,pki->ji',weight,candidate.adjoint(P,tau),E)
                        metrics['adjoint_pairing']=relative(lhs-rhs,lhs)
                        rawmat=(np.einsum('ij,pjn->pin',candidate.CH,E[:,:,:144])*np.sqrt(weight)[:,None,None]).reshape(-1,144)
                        modmat=(np.einsum('ij,pjn->pin',candidate.CH,R[:,:,:144])*np.sqrt(weight)[:,None,None]).reshape(-1,144)
                        scale=1/np.linalg.norm(rawmat,axis=0)
                        branks=candidate.ranks(svdvals(modmat*scale,check_finite=False))
                        gap=(modmat*scale).T@(modmat*scale)-(rawmat*scale).T@(rawmat*scale)
                        metrics['energy_upper_violation']=max(0.,float(eigvalsh(gap,check_finite=False)[-1]))
                        rigid=0.
                        for axis in np.eye(3):
                            v=np.zeros((16,3,3));v[:,0]=axis
                            for mode in [v.reshape(144),np.cross(np.broadcast_to(axis,(16,3,3)),q.reshape(16,3,3)).reshape(144)]:
                                rigid=max(rigid,float(np.linalg.norm(modmat@mode)/(np.linalg.norm(modmat)*np.linalg.norm(mode))))
                        metrics['rigid']=rigid
                        u=patch_vectors(q);exact=exact_patches(q,h,p)
                        metrics['raw_patch']=relative(E[:,:,:144]@u-exact,exact)
                        metrics['patch']=relative(R[:,:,:144]@u-exact,exact)
                        er=np.einsum('ij,pjc->pic',candidate.CH,exact);em=np.einsum('ij,pjc->pic',candidate.CH,R[:,:,:144]@u)
                        e0=np.einsum('p,pkc,pkc->c',weight,er,er);e1=np.einsum('p,pkc,pkc->c',weight,em,em)
                        metrics['patch_energy']=float(np.max(abs(e1/e0-1)))
                        axis_metrics={}
                        for mode,(other,transform) in alternate.items():
                            _,Ea=candidate.linear_raw(other['q'],h,p);Pa=candidate.project(other,h,Ea)
                            oldpoints=mapped_points(p,mode);_,Eold=candidate.linear_raw(q,h,oldpoints)
                            # Same original coefficient map, evaluated at the mapped points.
                            Sold=candidate.stress_field(state,oldpoints)
                            Rold=np.einsum('ij,pjn->pin',candidate.CI,Sold@P['a'])
                            mapped=Pa['R'].copy();mapped[:,:,:144]=Pa['R'][:,:,:144]@transform
                            axis_metrics[mode]=relative(mapped-Rold,Rold)
                        metrics['axes']=max(axis_metrics.values())
                        rng=np.random.default_rng(20261103);O,_=np.linalg.qr(rng.standard_normal((len(H),len(H))))
                        other=dict(state,T=state['T']@O,coeff=state['coeff']@O);OP=candidate.project(other,h,E)
                        metrics['basis_invariance']=relative(OP['R']-R,R)
                        # Nonlinear energy, raw-B force and full material+geometric tangent.
                        amplitude=np.array([state['length'],1.,1/state['length']])[:,None]
                        v1=(rng.standard_normal((16,3,3))*amplitude[None]).reshape(144)
                        v2=(rng.standard_normal((16,3,3))*amplitude[None]).reshape(144)
                        current=q+.002*v1
                        energy,force,material,geometric=candidate.variations(state,current,h,v1,v2)
                        _,_,Fcur,_,_=candidate.green(state,current,h,p)
                        metrics['current_orientation_violation']=float(np.min(np.linalg.det(Fcur))<=0.)
                        ec,_,_,_,_,_=candidate.energy_force(state,current+1e-25j*v1,h)
                        _,fc,_,_,_,_=candidate.energy_force(state,current+1e-25j*v2,h)
                        metrics['energy_force_derivative']=float(abs(ec.imag/1e-25-force@v1)/max(1.,abs(force@v1)))
                        metrics['full_tangent_derivative']=float(abs(v1@(fc.imag/1e-25)-material-geometric)/max(1.,abs(material+geometric)))
                        axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis)
                        W=np.array([[0.,-axis[2],axis[1]],[axis[2],0.,-axis[0]],[-axis[1],axis[0],0.]])
                        Qrigid=np.eye(3)+np.sin(.37)*W+(1-np.cos(.37))*(W@W)
                        qrigid=np.einsum('ij,npj->npi',Qrigid,q.reshape(16,3,3));qrigid[:,0]+=state['length']*np.array([.3,-.2,.1])
                        ur,fr,_,_,_,_=candidate.energy_force(state,qrigid.reshape(144),h)
                        metrics['finite_rigid']=max(float(abs(ur)),float(np.linalg.norm(fr)))
                        # Separate evaluation points, keeping the 9^3 operator fixed.
                        p11,w11=candidate.BASE.quadrature(11);weight11,volume11=candidate.measure(state,h,p11,w11)
                        S11=candidate.stress_field(state,p11);R11=np.einsum('ij,pjn->pin',candidate.CI,S11@P['a'])
                        _,E11=candidate.linear_raw(q,h,p11)
                        metrics['work11']=relative(np.einsum('p,pkn->kn',weight11,R11-E11),np.einsum('p,pkn->kn',weight11,E11))
                        mod11=(np.einsum('ij,pjn->pin',candidate.CH,R11[:,:,:144])*np.sqrt(weight11)[:,None,None]).reshape(-1,144)
                        branks11=candidate.ranks(svdvals(mod11*scale,check_finite=False))
                        checks=dict(condition=metrics['rcond']>=1e-10,solve=metrics['solve']<=1e-10,
                            local_rank=branks==branks11==[138]*3,
                            fields=all(value<=1e-9 for name,value in metrics.items() if name not in ['rcond','solve']))
                        artifact=out/f'k{k:g}-n{n}-c{cell}-h{h:g}.npz'
                        np.savez_compressed(artifact,a=P['a'],H=H,rhs=P['rhs'],B_singular=svdvals(modmat*scale,check_finite=False),
                            scale=scale,patch_vectors=u,raw_patch_energy=e0,modified_patch_energy=e1,
                            nonlinear_current=current,direction1=v1,direction2=v2,material_tangent=material,geometric_tangent=geometric)
                        checks['resources']=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3
                        records.append(dict(kappa=k,n=n,amplitude=amp,cell=cell,h=h,B_ranks=branks,B_ranks11=branks11,
                            metrics=metrics,axis_metrics=axis_metrics,checks={key:bool(value) for key,value in checks.items()},
                            state_artifact=states[-1]['artifact'],state_sha256=states[-1]['artifact_sha256'],
                            artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact),seconds=time.monotonic()-tick))
                        print(k,n,cell,h,branks,metrics,checks,flush=True)
                        if not all(checks.values()):raise RuntimeError('required_GSP_local_gate_failed')
    except Exception as exc:
        stop=repr(exc)
        if isinstance(exc,candidate.PreparationFailure):
            path=out/'stress-preparation-failure.npz';np.savez_compressed(path,**exc.witness)
            failure=dict(report=exc.report,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))
    sources=[ROOT/"rd/structure/main-shell/gsp-q3/experiments/gsp_q3_basis_v2.py",ROOT/"rd/structure/main-shell/gsp-q3/experiments/gsp_q3.py",Path(__file__),Path(candidate.__file__),Path(candidate.reference.__file__),Path(geometry.__file__),
        ROOT/'rd/structure/main-shell/gsp-q3/gsp-q3-formulation.md',ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-local-v2-protocol.md',ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-local-v2-preparation.json']
    dest.write_text(json.dumps(dict(records=records,states=states,preparation_failure=failure,
        all_checks_pass=len(records)==30 and stop is None,stop_reason=stop,physical_adoption=False,shared_bending_evaluated=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')


if __name__=='__main__':main()
