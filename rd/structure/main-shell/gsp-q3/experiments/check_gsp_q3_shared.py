"""Four fixed GSP shared thin-space comparisons, never a finite-h bound."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
from scipy.linalg import qr
import gsp_q3_shared_space as method
import check_q3_thin_limit as spectra
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def anchor(a,V,scale,k):
    value=method.common.field('quadratic',k,a['xy'])[a['free'],:6].reshape(-1)/scale
    z=value/np.linalg.norm(value);Q,_=qr((V.T@z)[:,None],mode='full',check_finite=False)
    new=np.column_stack((z,V@Q[:,1:]));return new,float(np.linalg.norm(new@new.T-V@V.T)),float(np.linalg.norm(new.T@new-np.eye(new.shape[1])))

def main():
    dest=ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-shared-result.json';assert not dest.exists()
    start=time.monotonic();cache=method.load_cache();maps={};out=ROOT/'build/analysis'/('gsp-q3-shared-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    records=[];stop=None;baseline=json.loads((ROOT/'rd/structure/main-shell/mq3/records/mq3i1-shared-result.json').read_text())
    try:
        for k,amp in itertools.product([0.,.4],[0.,.2]):
            tick=time.monotonic();a=method.assemble(k,2,amp,9,cache,maps);b=method.assemble(k,2,amp,11,cache,maps)
            scale=1/np.linalg.norm(a['A0'],axis=0)
            s9,r9,V9=spectra.rank_info(a['A0']*scale,True);s11,r11,V11=spectra.rank_info(b['A0']*scale,True);r=r9[1]
            Nv,shift9,orth9=anchor(a,V9[r:].T,scale,k);No,shift11,orth11=anchor(b,V11[r11[1]:].T,scale,k);N=scale[:,None]*Nv
            projection=float(np.linalg.norm(Nv@Nv.T-No@No.T));kernel=float(np.linalg.norm((a['A0']*scale)@Nv)/(np.linalg.norm(a['A0']*scale)*np.linalg.norm(Nv)))
            checks=dict(rank=r9==r11==[r]*3,kernel=kernel<=1e-10,projector=projection<=1e-7,
                anchor=max(shift9,shift11)<=1e-7 and max(orth9,orth11)<=1e-10,zero_c=max(a['zero_c'],b['zero_c'])<=1e-9)
            control=method.common.field('quadratic',k,a['xy']);rd=control[a['free'],:6].reshape(-1);cc=control[a['free'],6:].reshape(-1)
            control0=float(np.linalg.norm(a['A0']@rd)/max(1.,np.linalg.norm(rd)));control1=float(np.linalg.norm(a['A1']@rd+a['D1']@cc-a['targets']['quadratic'])/np.linalg.norm(a['targets']['quadratic']))
            checks['positive_control']=max(control0,control1)<=1e-9
            stored=dict(xy=a['xy'],coords=a['coords'],ids=np.array([c[2] for c in a['cells']]),free=a['free'],kernel=N,
                common_scale=scale,full_s9=s9,full_s11=s11,control_u=control)
            solutions={};stage_error=None
            if all(checks.values()):
                try:
                    T=np.column_stack((a['A1']@N,a['D1']))
                    for kind in ['quadratic','quartic']:
                        rh=method.common.reference_h1(kind,k,40);rh32=method.common.reference_h1(kind,k,32)
                        q,fit=method.solve_semidefinite_ball(T,a['targets'][kind],N,a['M'],2*rh);u=method.common.full_field(a,N,q)
                        rd=u[a['free'],:6].reshape(-1);cc=u[a['free'],6:].reshape(-1)
                        error=float(np.linalg.norm(b['A1']@rd+b['D1']@cc-b['targets'][kind])/np.linalg.norm(b['targets'][kind]))
                        zero=float(np.linalg.norm(b['A0']@rd)/(np.linalg.norm(b['A0'])*np.linalg.norm(rd)));h1=method.common.h1tool.h1(u,b['free'],b['M'])
                        sc=dict(kkt=fit['kkt']<=1e-8,sphere=h1<=2*rh*(1+1e-8),quadrature=abs(error-fit['error'])<=1e-8,
                            H1=abs(h1/fit['h1']-1)<=1e-8,zero=zero<=1e-9,reference=abs(rh32/rh-1)<=1e-12,positive_fit=kind!='quadratic' or error<=1e-8)
                        checks[kind]=all(sc.values());solutions[kind]=dict(fit=fit,observed_error=error,observed_h1=h1,reference_h1=rh,
                            h1_bound=2*rh,zero_constraint=zero,checks={name:bool(value) for name,value in sc.items()});stored['u_'+kind]=u
                except Exception as exc:stage_error=repr(exc);checks['first_order_solve']=False
            path=out/f'k{k:g}-a{amp:g}-n2.npz';np.savez_compressed(path,**stored)
            old=next(v for v in baseline['records'] if (v['kappa'],v['amplitude'],v['n'])==(k,amp,2))
            olderror=old['solutions']['quartic']['observed_error'];newerror=solutions.get('quartic',{}).get('observed_error')
            reduction=1-newerror/olderror if newerror is not None else None
            checks['resources']=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3
            records.append(dict(kappa=k,amplitude=amp,n=2,free_nodes=len(a['free']),rd_columns=a['A0'].shape[1],full_ranks9=r9,full_ranks11=r11,
                nullity=N.shape[1],anchor_shifts=[shift9,shift11],anchor_orthogonality=[orth9,orth11],kernel_residual=kernel,
                projector_difference=projection,control0=control0,control1=control1,solutions=solutions,stage_error=stage_error,
                baseline_quartic_error=olderror,relative_reduction=reduction,progress_gate_passed=bool(reduction>=.01) if k and reduction is not None else None,
                checks={name:bool(value) for name,value in checks.items()},seconds=time.monotonic()-tick,map_keys=a['map_keys'],
                artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path)))
            print('result',k,amp,'nullity',N.shape[1],'errors',{kind:s['observed_error'] for kind,s in solutions.items()},'checks',checks,flush=True)
            if not all(checks.values()):stop='required_shared_space_gate_failed';break
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(method.__file__),Path(method.candidate.__file__),ROOT/'rd/structure/main-shell/gsp-q3/experiments/gsp_q3_basis_v2.py',
        Path(method.common.__file__),Path(method.common.h1tool.__file__),Path(spectra.__file__),ROOT/'rd/structure/main-shell/mq3/experiments/mq3i1_shared_space.py',
        ROOT/'rd/structure/main-shell/gsp-q3/gsp-q3-shared-and-membrane-contract.md',ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-shared-protocol.md',ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-shared-preparation.json',
        ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-additional-local-result.json',ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-additional-local-audit.json',
        ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-local-v2-result.json',ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-local-v2-audit.json',ROOT/'rd/structure/main-shell/mq3/records/mq3i1-shared-result.json']
    dest.write_text(json.dumps(dict(records=records,maps=maps,all_checks_pass=len(records)==4 and stop is None,
        both_curved_improved=len(records)==4 and all(v['progress_gate_passed'] for v in records if v['kappa']),
        stop_reason=stop,physical_adoption=False,finite_thickness_lower_bound=False,finite_thickness_accuracy_evaluated=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
