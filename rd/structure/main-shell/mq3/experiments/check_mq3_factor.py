"""M-Q3 shared zero-kernel and H1-bounded first-strain comparison."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
import mq3_factor_space as method
import check_q3_thin_limit as spectra
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    target=ROOT/'rd/structure/main-shell/mq3/records/mq3-factor-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();cache={};paths={};out=ROOT/'build/analysis'/('mq3-factor-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();mapdir=out/'maps';mapdir.mkdir();records=[];stop=None
    try:
        for metric,work,k,amp,n in [(metric,work,k,amp,n) for metric,work in itertools.product(["physical","covariant"],[True,False]) for k,amp,n in [(.4,0.,1),(.4,0.,2),(.4,.2,1),(.4,.2,2),(0.,.2,2)]]:
            method.family.configure(metric,work)
            tick=time.monotonic();print('assemble',k,amp,n,flush=True);a=method.assemble(k,n,amp,9,cache,paths,mapdir);b=method.assemble(k,n,amp,11,cache,paths,mapdir);scale=1/np.linalg.norm(a['A0'],axis=0);s9,r9,V9=spectra.rank_info(a['A0']*scale,True);s11,r11,V11=spectra.rank_info(b['A0']*scale,True);pr9,ps9=spectra.subsets(a['C0'],scale);pr11,ps11=spectra.subsets(b['C0'],scale);r=r9[1];Nv=V9[r:].T;N=scale[:,None]*Nv;No=V11[r11[1]:].T;projection=float(np.linalg.norm(Nv@Nv.T-No@No.T));kernel=float(np.linalg.norm((a['A0']*scale)@Nv)/(np.linalg.norm(a['A0']*scale)*np.linalg.norm(Nv)));stable=r9==r11==[r]*3 and all(pr9[key]==pr11[key]==[pr9[key][1]]*3 for key in pr9);checks=dict(rank=stable,kernel=kernel<=1e-10,projector=projection<=1e-7,zero_c=max(a['zero_c'],b['zero_c'])<=1e-9);control=method.common.field('quadratic',k,a['xy']);rd=control[a['free'],:6].reshape(-1);c=control[a['free'],6:].reshape(-1);control0=float(np.linalg.norm(a['A0']@rd)/max(1.,np.linalg.norm(rd)));control1=float(np.linalg.norm(a['A1']@rd+a['D1']@c-a['targets']['quadratic'])/np.linalg.norm(a['targets']['quadratic']));checks['positive_control']=max(control0,control1)<=1e-9;solutions={};stored=dict(xy=a['xy'],coords=a['coords'],ids=np.array([v[2] for v in a['cells']]),free=a['free'],kernel=N,common_scale=scale,full_s9=s9,full_s11=s11,**{'s9_'+key:value for key,value in ps9.items()},**{'s11_'+key:value for key,value in ps11.items()});stage_error=None
            if all(checks.values()):
                try:
                    T=np.column_stack((a['A1']@N,a['D1']))
                    for kind in ['quadratic','quartic']:
                        rh=method.common.reference_h1(kind,k,40);rh32=method.common.reference_h1(kind,k,32);q,fit=method.solve_semidefinite_ball(T,a['targets'][kind],N,a['M'],2*rh);u=method.common.full_field(a,N,q);rd=u[a['free'],:6].reshape(-1);c=u[a['free'],6:].reshape(-1);error=float(np.linalg.norm(b['A1']@rd+b['D1']@c-b['targets'][kind])/np.linalg.norm(b['targets'][kind]));zero=float(np.linalg.norm(b['A0']@rd)/(np.linalg.norm(b['A0'])*np.linalg.norm(rd)));h1=method.common.h1tool.h1(u,b['free'],b['M']);sc=dict(kkt=fit['kkt']<=1e-8,sphere=h1<=2*rh*(1+1e-8),quadrature=abs(error-fit['error'])<=1e-8,H1=abs(h1/fit['h1']-1)<=1e-8,zero=zero<=1e-9,reference=abs(rh32/rh-1)<=1e-12,positive_fit=kind!='quadratic' or error<=1e-8);checks[kind]=all(sc.values());solutions[kind]=dict(fit=fit,observed_error=error,observed_h1=h1,reference_h1=rh,h1_bound=2*rh,zero_constraint=zero,checks={key:bool(v) for key,v in sc.items()});stored['u_'+kind]=u
                except Exception as exc:stage_error=repr(exc);checks['first_order_solve']=False
            artifact=out/f'{metric}-w{int(work)}-k{k:g}-a{amp:g}-n{n}.npz';np.savez_compressed(artifact,**stored);
            if metric=='physical' and work and solutions:
                old=json.loads((ROOT/'rd/structure/main-shell/mq3/records/mq3-shared-thin-space-result.json').read_text());ref=next(v for v in old['records'] if (v['kappa'],v['amplitude'],v['n'])==(k,amp,n));checks['regression']=N.shape[1]==ref['nullity'] and all(abs(solutions[kind]['observed_error']-ref['solutions'][kind]['observed_error'])<=1e-8 for kind in solutions)
            checks['resources']=time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3;records.append(dict(metric=metric,work=work,kappa=k,amplitude=amp,n=n,free_nodes=len(a['free']),rd_columns=a['A0'].shape[1],full_ranks9=r9,full_ranks11=r11,partial_ranks9=pr9,partial_ranks11=pr11,nullity=N.shape[1] if stable else None,kernel_residual=kernel,projector_difference=projection,control0=control0,control1=control1,solutions=solutions,stage_error=stage_error,checks={key:bool(v) for key,v in checks.items()},seconds=time.monotonic()-tick,map_keys=a['map_keys'],artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact)));print('result',k,amp,n,'ranks',r9,r11,'null',records[-1]['nullity'],'errors',{key:value['observed_error'] for key,value in solutions.items()},'checks',checks,flush=True);del a,b,V9,V11,N,Nv,No,stored
            if not all(checks.values()):stop='required_shared_space_gate_failed';break
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(method.__file__),Path(method.family.__file__),Path(method.checks.__file__),Path(method.common.__file__),Path(spectra.__file__),ROOT/'rd/structure/main-shell/mq3/records/mq3-zero-first-map-result.json',ROOT/'rd/structure/main-shell/mq3/mq3-constraint-factor-diagnostic.md',ROOT/'rd/structure/main-shell/mq3/records/mq3-factor-preparation.json']
    target.write_text(json.dumps(dict(records=records,maps=paths,all_checks_pass=len(records)==20 and stop is None,stop_reason=stop,physical_adoption=False,finite_thickness_lower_bound=False,finite_thickness_accuracy_evaluated=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
