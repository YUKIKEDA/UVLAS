"""Frozen local and four-configuration mixed-order kinematic screen."""
from pathlib import Path
import hashlib,json,time,itertools
import numpy as np
from scipy.linalg import svdvals,svd,qr
import mixed_order_150 as m
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm(a)/max(1.,np.linalg.norm(b)))
def ranks(s):return [int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]]

def local_gate(mesh,k,cell,out):
    ids=mesh['cells'][cell];q=mesh['ref'][ids].reshape(150);p,w=m.quadrature(11,2);J,B0,B1=m.expansion(q,p)
    _,Ji,Bi=m.raw(q,p,0.,True);_,_,Bc=m.raw(q,p,1e-25j,True)
    metrics=dict(zero_path=rel(B0-Bi,B0),first_derivative=rel(B1-Bc.imag/1e-25,B1),zero_c=float(np.max(abs(B0[:,:,123:]))))
    Nr,_=m.shape(p,0);x=(Nr@q.reshape(50,3)[:25])[:,0]
    for kind in ['quadratic','cubic']:
        u=m.displacement(mesh,kind,k)[ids].reshape(150);target=m.target(kind,k,x)
        metrics[kind+'_zero']=rel(B0@u,np.zeros((len(p),6)));metrics[kind+'_first']=rel(B1@u-target,target)
    A=np.column_stack([np.einsum('ij,nj->ni',mb,q.reshape(50,3)).reshape(150) for mb in m.MB])
    rigid=[]
    for axis in np.eye(3):
        u=np.zeros((50,3));u[:25]=axis;rigid.append(u.reshape(150))
        rigid.append(np.cross(axis,q.reshape(50,3)).reshape(150))
    rigid=np.array(rigid).T
    finite=[];saved=dict(q=q)
    for index,h in enumerate([.01,.001,.0001]):
        p9,w9=m.quadrature(9);_,J9,B9=m.raw(q,p9,h*p9[:,2]/2);weight=w9*np.linalg.det(J9)
        E=(np.einsum('ij,pjn->pin',m.CH,B9)*np.sqrt(weight)[:,None,None]).reshape(-1,150);scale=1/np.linalg.norm(E,axis=0)
        s9=svdvals(E*scale,check_finite=False);p11,w11=m.quadrature(11);_,J11,B11=m.raw(q,p11,h*p11[:,2]/2)
        E11=(np.einsum('ij,pjn->pin',m.CH,B11)*np.sqrt(w11*np.linalg.det(J11))[:,None,None]).reshape(-1,150)
        s11=svdvals(E11*scale,check_finite=False);surface=m.boundary(q,h)
        volume=h/2*np.einsum('p,pki->ik',weight,B9)
        fm=dict(min_det=float(min(np.min(np.linalg.det(J9)),np.min(np.linalg.det(J11)))),constant=rel(B9@A-np.eye(6),np.eye(6)),
            rigid=rel(B9@rigid,B9),work=rel(volume-surface,surface))
        checks=dict(rank=ranks(s9)==ranks(s11)==[144]*3,orientation=fm['min_det']>0,fields=max(fm['constant'],fm['rigid'],fm['work'])<=1e-9)
        finite.append(dict(h=h,ranks9=ranks(s9),ranks11=ranks(s11),metrics=fm,checks={a:bool(b) for a,b in checks.items()}))
        saved[f'scale{index}']=scale;saved[f's9_{index}']=s9;saved[f's11_{index}']=s11
    path=out/f'k{k:g}-local-{cell}.npz';np.savez_compressed(path,**saved)
    passed=max(metrics.values())<=1e-9 and all(all(f['checks'].values()) for f in finite)
    return dict(kappa=k,cell=cell,metrics=metrics,finite=finite,passed=bool(passed),artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))

def main():
    dest=ROOT/'docs/verification/mixed-order-150-result.json';assert not dest.exists();start=time.monotonic()
    out=ROOT/'build/analysis'/('mixed-order-150-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();records=[];locals=[];stop=None
    baseline=json.loads((ROOT/'docs/verification/mq3i1-shared-result.json').read_text())
    try:
        for k,amp in itertools.product([0.,.4],[0.,.2]):
            tick=time.monotonic();mesh=m.mesh(k,2,amp);folder=out/f'k{k:g}-a{amp:g}';folder.mkdir()
            for cell in range(4):
                result=local_gate(mesh,k,cell,folder);result['amplitude']=amp;locals.append(result)
                print('local',k,amp,cell,result['passed'],flush=True)
                if not result['passed']:raise RuntimeError('mixed_order_local_gate_failed')
            a=m.assemble(mesh,k,9);b=m.assemble(mesh,k,11);scale=1/np.linalg.norm(a['A0'],axis=0)
            _,s9,V9=svd(a['A0']*scale,full_matrices=False,check_finite=False);_,s11,V11=svd(b['A0']*scale,full_matrices=False,check_finite=False)
            r=ranks(s9)[1];V=V9[r:].T;Vo=V11[ranks(s11)[1]:].T
            controls=np.column_stack([m.displacement(mesh,kind,k)[a['free_rd']].reshape(-1)/scale for kind in ['quadratic','cubic']])
            Z,_=qr(controls,mode='economic',check_finite=False);rotation,_=qr(V.T@Z,mode='full',check_finite=False)
            Vn=np.column_stack((Z,V@rotation[:,2:]));shift=float(np.linalg.norm(Vn@Vn.T-V@V.T));orth=float(np.linalg.norm(Vn.T@Vn-np.eye(Vn.shape[1])));N=scale[:,None]*Vn
            projection=float(np.linalg.norm(Vn@Vn.T-Vo@Vo.T));kernel=float(np.linalg.norm((a['A0']*scale)@Vn)/(np.linalg.norm(a['A0']*scale)*np.linalg.norm(Vn)))
            checks=dict(rank=ranks(s9)==ranks(s11)==[r]*3,anchor=shift<=1e-7 and orth<=1e-10,projector=projection<=1e-7,kernel=kernel<=1e-10,zero_c=max(a['zero_c'],b['zero_c'])<=1e-9)
            solutions={};F=m.metric_factor(a,N,mesh);T=np.column_stack((a['A1']@N,a['D1']))
            stored=dict(ref=mesh['ref'],cells=mesh['cells'],free_rd=a['free_rd'],free_c=a['free_c'],kernel=N,scale=scale,full_s9=s9,full_s11=s11,
                xy0=mesh['xy'][0],xy1=mesh['xy'][1],xy2=mesh['xy'][2],field0=mesh['field_ids'][0],field1=mesh['field_ids'][1],field2=mesh['field_ids'][2])
            stage_error=None
            if all(checks.values()):
                try:
                    for kind in ['quadratic','cubic','quartic']:
                        rh=m.reference_h1(kind,k,40);rh32=m.reference_h1(kind,k,32)
                        fitq,fit=m.solve_ball(T,a['targets'][kind],F,2*rh);u=np.zeros_like(mesh['ref']);u[a['free_rd']]=(N@fitq[:N.shape[1]]).reshape(-1,3);u[a['free_c']]=fitq[N.shape[1]:].reshape(-1,3)
                        rd=u[a['free_rd']].reshape(-1);cc=u[a['free_c']].reshape(-1)
                        error=float(np.linalg.norm(b['A1']@rd+b['D1']@cc-b['targets'][kind])/np.linalg.norm(b['targets'][kind]))
                        dc=u[np.r_[mesh['free'][1],mesh['free'][2]]];h1=float(np.sqrt(np.sum(dc*(b['M']@dc))))
                        zero=float(np.linalg.norm(b['A0']@rd)/(np.linalg.norm(b['A0'])*np.linalg.norm(rd)))
                        sc=dict(kkt=fit['kkt']<=1e-8,sphere=h1<=2*rh*(1+1e-8),quadrature=abs(error-fit['error'])<=1e-8,H1=abs(h1/fit['h1']-1)<=1e-8,
                            zero=zero<=1e-9,reference=abs(rh32/rh-1)<=1e-12,positive_fit=kind=='quartic' or error<=1e-8)
                        checks[kind]=all(sc.values());solutions[kind]=dict(fit=fit,observed_error=error,observed_h1=h1,reference_h1=rh,zero_constraint=zero,checks={a:bool(b) for a,b in sc.items()});stored['u_'+kind]=u
                except Exception as exc:stage_error=repr(exc);checks['first_order_solve']=False
            artifact=folder/'shared.npz';np.savez_compressed(artifact,**stored)
            old=next(x for x in baseline['records'] if (x['kappa'],x['amplitude'])==(k,amp))['solutions']['quartic']['observed_error']
            error=solutions.get('quartic',{}).get('observed_error');reduction=1-error/old if error is not None else None
            checks['resources']=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3
            records.append(dict(kappa=k,amplitude=amp,n=2,rd_columns=a['A0'].shape[1],full_ranks9=ranks(s9),full_ranks11=ranks(s11),nullity=N.shape[1],
                projector_difference=projection,anchor_shift=shift,anchor_orthogonality=orth,kernel_residual=kernel,solutions=solutions,stage_error=stage_error,
                baseline_quartic_error=old,relative_reduction=reduction,progress_gate=bool(reduction>=.01) if k and reduction is not None else None,
                checks={a:bool(b) for a,b in checks.items()},seconds=time.monotonic()-tick,artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact)))
            print('shared',k,amp,'nullity',N.shape[1],{kind:value['observed_error'] for kind,value in solutions.items()},checks,stage_error,flush=True)
            if not all(checks.values()):raise RuntimeError('mixed_order_shared_numerical_gate_failed')
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(m.__file__),ROOT/'docs/mixed-order-150-screening-contract.md',ROOT/'docs/verification/mixed-order-150-protocol.md',
        ROOT/'docs/verification/mixed-order-150-preparation.json',ROOT/'docs/verification/mq3i1-shared-result.json']
    dest.write_text(json.dumps(dict(records=records,local_records=locals,all_numerical_checks=len(records)==4 and stop is None,
        both_curved_improved=len(records)==4 and all(r['progress_gate'] for r in records if r['kappa']),stop_reason=stop,physical_adoption=False,
        finite_thickness_accuracy_evaluated=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
