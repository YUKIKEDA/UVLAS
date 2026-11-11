"""Stop-on-failure M-Q3I1 local solvability, conservation and membrane checks."""
from pathlib import Path
import hashlib
import itertools
import json
import time
import numpy as np
from scipy.linalg import qr, solve, svdvals
import mq3i1_moments as candidate
import q3_thin_expansion as geometry
import check_q9_plate_convergence as memory

ROOT = Path(__file__).resolve().parents[2]


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def relative(A, B):
    return float(np.linalg.norm(A)/max(1., np.linalg.norm(B)))


def patch_vectors(q):
    q = q.reshape(16,3,3)
    u = np.zeros((16,3,3,4))
    for col, (component, variable) in enumerate([(0,0), (1,1), (0,1), (1,0)]):
        r, d, c = q[:,0,variable], q[:,1,variable], q[:,2,variable]
        u[:,0,component,col] = r*r/2
        u[:,1,component,col] = r*d
        u[:,2,component,col] = d*d+r*c
    return u.reshape(144,4)


def exact_patches(q, h, points):
    N, _ = candidate.family.base.geometry.shape(points)
    z = h*points[:,2]/2
    values = q.reshape(16,3,3)
    xyz = N@values[:,0]+z[:,None]*(N@values[:,1])+z[:,None]**2/2*(N@values[:,2])
    result = np.zeros((len(points),6,4))
    for col,(component,variable) in enumerate([(0,0),(1,1),(0,1),(1,0)]):
        grad = np.zeros((len(points),3,3))
        grad[:,component,variable] = xyz[:,variable]
        E = (grad+grad.transpose(0,2,1))/2
        result[:,:,col] = np.einsum('kij,pij->pk', candidate.family.base.MB, E)
    return result


def reparameterize(q, mode):
    nodes = np.array(list(itertools.product(candidate.family.base.geometry.NODES,
        candidate.family.base.geometry.NODES)))[:,::-1]
    old = nodes[:,::-1] if mode=='swap' else nodes*np.array([-1.,1.] if mode=='s' else [1.,-1.])
    order = np.array([np.argmin(np.linalg.norm(nodes-v,axis=1)) for v in old])
    P = np.zeros((144,144))
    for new, prev in enumerate(order):
        for power in range(3):
            sign = -1. if power==1 else 1.
            for axis in range(3):
                P[9*new+3*power+axis,9*prev+3*power+axis] = sign
    return P@q, P


def mapped_points(points, mode):
    out = points.copy()
    if mode=='swap':
        out[:,:2] = points[:,[1,0]]
    else:
        out[:,0 if mode=='s' else 1] *= -1
    out[:,2] *= -1
    return out


def save_state(path, state):
    np.savez_compressed(path, **{key:state[key] for key in ['q','T','rows','U','Qs','hard_U',
        'length','geometry_projection','geometry_vectors','A','b']})


def differential(state):
    X0, S0, b0, F0, g0 = candidate.coefficients(state, 0.)
    Sc, bc, _, _ = candidate.system(state, 1e-25j)
    Xd = solve(S0, bc.imag/1e-25-Sc.imag/1e-25@X0, assume_a='sym', check_finite=False)
    Xc, _, _, _, _ = candidate.coefficients(state, 1e-25j)
    p, w = candidate.family.base.quadrature(9)
    R0 = candidate.evaluate(state, 0., X0, p)
    Rd = candidate.evaluate(state, 1e-25j, X0+1e-25j*Xd, p).imag/1e-25
    alt = candidate.evaluate(state, 1e-25j, Xc, p).imag/1e-25
    c0 = state['T']@X0[:state['T'].shape[1]]
    cd = state['T']@Xd[:state['T'].shape[1]]
    metrics = dict(derivative=relative(Rd-alt, Rd),
        solve=relative(S0@X0-b0,b0),
        derivative_solve=relative(S0@Xd+Sc.imag/1e-25@X0-bc.imag/1e-25,bc.imag/1e-25),
        zero_moments=relative(F0@c0-g0,g0))
    metrics['zero_c'] = float(np.max(np.abs(R0[:,:,geometry.CC])))
    metrics['constant'] = max(float(np.max(abs(R0[:,:,144:]-np.eye(6)))),float(np.max(abs(Rd[:,:,144:]))))
    shape0, shape1 = R0.reshape(9,9,9,6,150), Rd.reshape(9,9,9,6,150)
    zp, zw = np.polynomial.legendre.leggauss(9)
    mean = np.einsum('z,abzkn->abkn',zw/2,shape0)
    slope = np.einsum('z,z,abzkn->abkn',zw,zp,shape1)/np.sum(zw*zp*zp)
    metrics['zero_dependence'] = relative(shape0-mean[:,:,None],shape0)
    metrics['first_linearity'] = relative(shape1-zp[None,None,:,None,None]*slope[:,:,None],shape1)
    for h in [-.0001,.0001]:
        X, S, b, F, g = candidate.coefficients(state,h)
        coef = state['T']@X[:state['T'].shape[1]]
        metrics['near_'+str(h)] = max(relative(S@X-b,b), relative(F@coef-g,g))
        C,_,_,_,_,_=candidate.gram(state,h)
        C=state['hard_U']@C
        Q,_=qr(C.T,mode='full',check_finite=False)
        ranksC=candidate.ranks(svdvals(C,check_finite=False))
        ranksA=candidate.ranks(svdvals(state['A']@Q[:,len(C):],check_finite=False))
        metrics['near_rank_'+str(h)] = float(ranksC!=[len(C)]*3 or ranksA!=[C.shape[1]-len(C)]*3)
    return X0, Xd, metrics


def main():
    dest = ROOT/'docs/verification/mq3i1-local-result.json'
    if dest.exists():
        raise FileExistsError(dest)
    start = time.monotonic()
    out = ROOT/'build/analysis'/('mq3i1-local-'+time.strftime('%Y%m%d-%H%M%S'))
    out.mkdir()
    records, states, stop, failure = [], [], None, None
    try:
        for k in [0.,.4]:
            for n,amp in [(1,0.),(2,.2)]:
                xy,cells,coords,_ = geometry.element.mesh(k,n,amp)
                for cell,(_,_,ids) in enumerate(cells):
                    q = coords[ids].reshape(144)
                    state = candidate.choose(q)
                    state_path = out/f'k{k:g}-n{n}-c{cell}-state.npz'
                    save_state(state_path,state)
                    _,_,diff = differential(state)
                    state_record = dict(kappa=k,n=n,amplitude=amp,cell=cell,
                        preparation=state['preparation'], differential=diff,
                        artifact=state_path.relative_to(ROOT).as_posix(),artifact_sha256=digest(state_path))
                    states.append(state_record)
                    if not all(state['preparation']['checks'].values()) or max(diff.values())>1e-9:
                        raise RuntimeError('preparation_or_differential_failed')
                    alternate = {}
                    for mode in ['swap','s','t']:
                        qt,P = reparameterize(q,mode)
                        alternate[mode] = (candidate.choose(qt),P)
                    for h in [.01,.001,.0001]:
                        tick = time.monotonic()
                        X,S,b,F,g = candidate.coefficients(state,h)
                        coef = state['T']@X[:state['T'].shape[1]]
                        singular = svdvals(S,check_finite=False)
                        points,w = candidate.family.base.quadrature(9)
                        R = candidate.evaluate(state,h,X,points)
                        raw = candidate.family.data(q,h,points)
                        measure = w*np.linalg.det(raw['J'])
                        weight = measure/measure.sum()
                        sr = (np.einsum('ij,pjn->pin',candidate.family.base.CH,raw['E'][:,:,:144])*np.sqrt(weight)[:,None,None]).reshape(-1,144)
                        sm = (np.einsum('ij,pjn->pin',candidate.family.base.CH,R[:,:,:144])*np.sqrt(weight)[:,None,None]).reshape(-1,144)
                        scale = 1/np.linalg.norm(sr,axis=0)
                        bs = svdvals(sm*scale,check_finite=False)
                        branks = candidate.ranks(bs)
                        metrics = dict(rcond=float(singular[-1]/singular[0]),
                            solve=relative(S@X-b,b), all_moments=relative(F@coef-g,g),
                            constant=float(np.max(abs(R[:,:,144:]-np.eye(6)))),
                            work=relative(np.einsum('p,pkn->kn',weight,R-raw['E']),np.einsum('p,pkn->kn',weight,raw['E'])))
                        C,_,_,_,_,_=candidate.gram(state,h)
                        C=state['hard_U']@C
                        QC,_=qr(C.T,mode='full',check_finite=False)
                        ranksC=candidate.ranks(svdvals(C,check_finite=False))
                        ranksA=candidate.ranks(svdvals(state['A']@QC[:,len(C):],check_finite=False))
                        target_raw=(np.einsum('ij,pjn->pin',candidate.family.base.CH,raw['W'])*np.sqrt(weight)[:,None,None]).reshape(-1,117)
                        tn=np.linalg.norm(target_raw,axis=0);active=tn>1e-14*tn.max()
                        target_ranks=candidate.ranks(svdvals(target_raw[:,active]/tn[active],check_finite=False))
                        rigid = 0.
                        for axis in np.eye(3):
                            trans = np.zeros((16,3,3));trans[:,0]=axis
                            rot = np.cross(np.broadcast_to(axis,(16,3,3)),q.reshape(16,3,3))
                            for v in [trans.reshape(144),rot.reshape(144)]:
                                rigid = max(rigid,float(np.linalg.norm(sm@v)/(np.linalg.norm(sm)*np.linalg.norm(v))))
                        metrics['rigid'] = rigid
                        covR = candidate.family.base.covariant(raw['J'],R)
                        covE = candidate.family.base.covariant(raw['J'],raw['E'])
                        mid = points.copy();mid[:,2]=0.
                        Jmid,_ = candidate.family.raw(q,h,mid)
                        rmid = candidate.family.base.covariant(Jmid,candidate.evaluate(state,h,X,mid))
                        emid = candidate.family.base.covariant(Jmid,candidate.family.raw(q,h,mid)[1])
                        metrics['zz'] = relative(covR[:,2,2]-covE[:,2,2],covE[:,2,2])
                        metrics['thickness_difference'] = relative((covR[:,:2,:2]-rmid[:,:2,:2])-(covE[:,:2,:2]-emid[:,:2,:2]),covE[:,:2,:2])
                        u = patch_vectors(q)
                        exact = exact_patches(q,h,points)
                        metrics['raw_patch'] = relative(raw['E'][:,:,:144]@u-exact,exact)
                        metrics['patch'] = relative(R[:,:,:144]@u-exact,exact)
                        Eraw = np.einsum('ij,pjk->pik',candidate.family.base.CH,exact)
                        Emod = np.einsum('ij,pjk->pik',candidate.family.base.CH,R[:,:,:144]@u)
                        e0 = np.einsum('p,pkc,pkc->c',weight,Eraw,Eraw)
                        e1 = np.einsum('p,pkc,pkc->c',weight,Emod,Emod)
                        metrics['patch_energy'] = float(np.max(abs(e1/e0-1)))
                        target = (raw['W']@state['T']*np.sqrt(weight)[:,None,None]).reshape(-1,len(state['T'][0]))
                        Qtarget,_ = np.linalg.qr(target,mode='reduced')
                        fp = (raw['F'][:,:,:144]@u*np.sqrt(weight)[:,None,None]).reshape(-1,4)
                        metrics['patch_target_span'] = relative(fp-Qtarget@(Qtarget.T@fp),fp)
                        axis_metrics = {}
                        for mode,(other,P) in alternate.items():
                            Xa,_,_,_,_ = candidate.coefficients(other,h)
                            Ra = candidate.evaluate(other,h,Xa,points)
                            Rold = candidate.evaluate(state,h,X,mapped_points(points,mode))
                            mapped = Ra.copy();mapped[:,:,:144]=Ra[:,:,:144]@P
                            axis_metrics[mode] = relative(mapped-Rold,Rold)
                        metrics['axes'] = max(axis_metrics.values())
                        rng = np.random.default_rng(31415)
                        O,_ = np.linalg.qr(rng.standard_normal((len(state['A']),len(state['A']))))
                        temp = dict(state,A=O@state['A'],b=O@state['b'])
                        Xo,_,_,_,_ = candidate.coefficients(temp,h)
                        metrics['test_basis'] = relative(Xo-X,X)
                        checks = dict(condition=metrics['rcond']>=1e-10,solve=metrics['solve']<=1e-10,
                            moments=metrics['all_moments']<=1e-10,local_rank=branks==[138]*3,
                            coefficient_ranks=ranksC==[len(C)]*3 and ranksA==[C.shape[1]-len(C)]*3 and target_ranks==state['preparation']['space_ranks'],
                            fields=all(value<=1e-9 for key,value in metrics.items() if key not in ['rcond','solve','all_moments']))
                        p11,w11 = candidate.family.base.quadrature(11)
                        R11 = candidate.evaluate(state,h,X,p11)
                        d11 = candidate.family.data(q,h,p11)
                        weight11 = w11*np.linalg.det(d11['J']);weight11/=weight11.sum()
                        work11 = np.einsum('p,pkn->kn',weight11,R11-d11['E'])
                        metrics['work11'] = relative(work11,np.einsum('p,pkn->kn',weight11,d11['E']))
                        sm11=(np.einsum('ij,pjn->pin',candidate.family.base.CH,R11[:,:,:144])*np.sqrt(weight11)[:,None,None]).reshape(-1,144)
                        branks11=candidate.ranks(svdvals(sm11*scale,check_finite=False))
                        checks['separate_points'] = metrics['work11']<=1e-9 and branks11==[138]*3
                        artifact = out/f'k{k:g}-n{n}-c{cell}-h{h:g}.npz'
                        np.savez_compressed(artifact,X=X,S=S,rhs=b,B_singular=bs,patch_vectors=u,
                            patch_energies_raw=e0,patch_energies_modified=e1,orthogonal_test_basis=O)
                        checks['resources'] = time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3
                        records.append(dict(kappa=k,n=n,amplitude=amp,cell=cell,h=h,B_ranks=branks,B_ranks11=branks11,
                            C_ranks=ranksC,A_kernel_ranks=ranksA,target_ranks=target_ranks,
                            metrics=metrics,axis_metrics=axis_metrics,checks={key:bool(v) for key,v in checks.items()},
                            state_artifact=state_record['artifact'],state_sha256=state_record['artifact_sha256'],
                            artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact),seconds=time.monotonic()-tick))
                        print(k,n,cell,h,'rank',branks,'metrics',metrics,'checks',checks,flush=True)
                        if not all(checks.values()):
                            raise RuntimeError('required_local_gate_failed')
    except Exception as exc:
        stop = repr(exc)
        if isinstance(exc,candidate.PreparationFailure):
            path=out/'preparation-failure.npz'
            np.savez_compressed(path,**exc.witness)
            failure=dict(report=exc.report,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))
    sources = [Path(__file__),Path(candidate.__file__),Path(candidate.family.__file__),
        Path(candidate.interior.__file__),Path(geometry.__file__),
        ROOT/'docs/mq3-interior-moment-alternative.md',ROOT/'docs/verification/mq3i1-local-preparation.json']
    dest.write_text(json.dumps(dict(records=records,states=states,all_checks_pass=len(records)==30 and stop is None,
        stop_reason=stop,preparation_failure=failure,physical_adoption=False,shared_bending_evaluated=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},
        seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')


if __name__=='__main__':
    main()
