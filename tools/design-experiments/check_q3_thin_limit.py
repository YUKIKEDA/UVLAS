"""Finite shared-kernel diagnostic, not a finite-thickness lower-bound proof."""
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
from scipy.linalg import svd,svdvals
import q3_thin_expansion as method
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rank_info(A,vectors=False):
    if vectors:
        U,s,V=svd(A,full_matrices=False,check_finite=False);del U
    else:s=svdvals(A,check_finite=False);V=None
    ranks=[int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]];return s,ranks,V
def subsets(C,scale):
    C=C.reshape(-1,6,C.shape[1]);out={};spectra={}
    for label,rows in [('membrane',[0,1,5]),('shear',[3,4]),('normal',[2])]:
        spectrum,ranks,_=rank_info(C[:,rows].reshape(-1,C.shape[-1])*scale);out[label]=ranks;spectra[label]=spectrum
    return out,spectra
def main():
    dest=ROOT/'docs/verification/q3-thin-limit-result.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();out=ROOT/'build/analysis'/('q3-thin-limit-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();records=[];stop=None
    try:
        for k,amp,n in itertools.product([0.,.4],[0.,.2],[1,2,4]):
            tick=time.monotonic();a=method.assemble(k,n,amp,9);b=method.assemble(k,n,amp,11);scale=1/np.linalg.norm(a['A0'],axis=0);s9,r9,V9=rank_info(a['A0']*scale,True);s11,r11,V11=rank_info(b['A0']*scale,True);partial9,ps9=subsets(a['C0'],scale);partial11,ps11=subsets(b['C0'],scale);rank=r9[1];Nv=V9[rank:].T;N=scale[:,None]*Nv;Vother=V11[r11[1]:].T;projector=float(np.linalg.norm(Nv@Nv.T-Vother@Vother.T));kernel=float(np.linalg.norm((a['A0']*scale)@Nv)/(np.linalg.norm(a['A0']*scale)*np.linalg.norm(Nv)));rankcheck=r9==r11==[rank]*3 and all(partial9[key]==partial11[key]==[partial9[key][1]]*3 for key in partial9);solutions={};stored=dict(xy=a['xy'],coords=a['coords'],ids=np.array([v[2] for v in a['cells']]),free=a['free'],kernel=N,common_scale=scale,full_s9=s9,full_s11=s11,**{'s9_'+key:value for key,value in ps9.items()},**{'s11_'+key:value for key,value in ps11.items()});checks=dict(rank=rankcheck,kernel=kernel<=1e-10,projector=projector<=1e-7,zero_c=max(a['zero_c'],b['zero_c'])<=1e-12)
            control=method.field('quadratic',k,a['xy']);rd=control[a['free'],:6].reshape(-1);cc=control[a['free'],6:].reshape(-1);c0=float(np.linalg.norm(a['A0']@rd)/max(np.linalg.norm(rd),1e-300));c1=float(np.linalg.norm(a['A1']@rd+a['D1']@cc-a['targets']['quadratic'])/np.linalg.norm(a['targets']['quadratic']));checks['positive_control']=max(c0,c1)<=1e-9
            if all(checks.values()):
                T=np.column_stack((a['A1']@N,a['D1']));H=method.induced_metric(N,a['M'])
                for kind in ['quadratic','quartic']:
                    rh32=method.reference_h1(kind,k,32);rh=method.reference_h1(kind,k,40);q,fit=method.solve_ball(T,a['targets'][kind],H,2*rh);u=method.full_field(a,N,q);rdu=u[a['free'],:6].reshape(-1);cdu=u[a['free'],6:].reshape(-1);res=b['A1']@rdu+b['D1']@cdu-b['targets'][kind];observed=float(np.linalg.norm(res)/np.linalg.norm(b['targets'][kind]));zero=float(np.linalg.norm(b['A0']@rdu)/(np.linalg.norm(b['A0'])*np.linalg.norm(rdu)));h1=method.h1tool.h1(u,b['free'],b['M']);fitchecks=dict(kkt=fit['kkt']<=1e-8,sphere=h1<=2*rh*(1+1e-8),quadrature=abs(observed-fit['error'])<=1e-8,h1=abs(h1/fit['h1']-1)<=1e-8,zero_constraint=zero<=1e-9,reference_h1=abs(rh32/rh-1)<=1e-12,positive_fit=kind!='quadratic' or observed<=1e-8);checks[kind]=all(fitchecks.values());solutions[kind]=dict(fit=fit,observed_error=observed,observed_h1=h1,reference_h1=rh,h1_bound=2*rh,zero_constraint=zero,checks={key:bool(value) for key,value in fitchecks.items()});stored['u_'+kind]=u
            artifact=out/f'k{k:g}-a{amp:g}-n{n}.npz';np.savez_compressed(artifact,**stored);checks['resources']=time.monotonic()-tick<=120 and time.monotonic()-start<=600 and memory.memory_bytes()<=2*1024**3;records.append(dict(kappa=k,amplitude=amp,n=n,free_nodes=len(a['free']),rd_columns=a['A0'].shape[1],full_ranks9=r9,full_ranks11=r11,partial_ranks9=partial9,partial_ranks11=partial11,nullity=N.shape[1],kernel_residual=kernel,projector_difference=projector,control0=c0,control1=c1,solutions=solutions,checks={key:bool(value) for key,value in checks.items()},seconds=time.monotonic()-tick,artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact)));print(k,amp,n,'null',N.shape[1],'errors',{key:value['observed_error'] for key,value in solutions.items()},'checks',checks,flush=True)
            del a,b,V9,V11,T,H,N,Nv,stored
            if not all(checks.values()):stop='required_gate_failed';break
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(method.__file__),Path(method.element.__file__),Path(method.h1tool.__file__),ROOT/'docs/q3-thin-limit-compatibility-contract.md',ROOT/'docs/verification/q3-thin-limit-experiment.md']
    dest.write_text(json.dumps(dict(records=records,all_checks_pass=len(records)==12 and stop is None,stop_reason=stop,physical_adoption=False,finite_thickness_lower_bound=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
