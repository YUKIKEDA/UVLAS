"""Five fixed full-position necessary comparisons for MMS-1."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import svd,lstsq
import mixed_mean_stress_thin_v2 as t
import check_moment_enriched_thin_global_v2 as shared
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def case(k,n,amp,out):
    tick=time.monotonic();g=t.m.mesh(k,n,amp);cells=[]
    for index,ids in enumerate(g['cells']):
        cell=t.local(g['ref'][ids].reshape(-1),k);cells.append(cell);print('LOCAL',k,n,amp,index,'stress',cell['stress']['stress_rank'],'q01',cell['ny'],'constraint',cell['info']['ranks'],flush=True)
        if time.monotonic()-tick>900 or memory.memory_bytes()>2*1024**3:raise RuntimeError('MMS1_local_resource_stop')
    Z,interface,_,_=shared.assemble_constraints(cells,n);A=[];b=[]
    for cell in cells:
        E,target=t.sample(cell,k,13);A.append(E@Z[cell['yoffset']:cell['yoffset']+cell['ny']]);b.append(target)
    A=np.vstack(A);b=np.vstack(b);U,sv,V=svd(A,full_matrices=False,check_finite=False);ranks=t.rank_info(sv)
    if len(set(ranks))!=1:raise ValueError('MMS1_global_position_rank_ambiguous:'+str(ranks))
    rank=ranks[0];value=(V[:rank].T/sv[:rank])@(U[:,:rank].T@b);res=A@value-b;error=np.linalg.norm(res,axis=0)/np.linalg.norm(b,axis=0);kkt=np.linalg.norm(A.T@res,axis=0)/np.maximum(1e-30,np.linalg.norm(A)*np.linalg.norm(b,axis=0))
    latent=Z@value;err=np.zeros(3);ref=np.zeros(3);maxwork=0.;values=[]
    for cell in cells:
        y=latent[cell['yoffset']:cell['yoffset']+cell['ny']];E,target=t.sample(cell,k,15);rr=E@y-target;err+=np.sum(rr*rr,axis=0);ref+=np.sum(target*target,axis=0);values.append(cell['Z']@y);maxwork=max(maxwork,cell['residual'])
    ind=np.sqrt(err/ref);checks=dict(rank=True,work=maxwork<=1e-9,fluctuation=max(c['fluctuation'] for c in cells)<=1e-9,constant_stress=max(c['constant_stress'] for c in cells)<=1e-9,kkt=max(kkt)<=1e-8,quadrature=max(abs(ind-error))<=1e-8,resources=time.monotonic()-tick<=900 and memory.memory_bytes()<=2*1024**3)
    path=out/f'k{k:g}-n{n}-a{amp:g}.npz';np.savez_compressed(path,values=np.array(values),latent=latent,reference=g['ref'],cells=g['cells'])
    return dict(kappa=k,n=n,amplitude=amp,global_q0_ranks=ranks,q01_dimension=interface['nullity'],local=[dict(stress=c['stress'],q01=c['info'],odd=c['odd_TB'],gram_condition=c['gram_condition'],fluctuation=c['fluctuation'],constant_stress=c['constant_stress']) for c in cells],error=error.tolist(),independent=ind.tolist(),kkt=kkt.tolist(),checks={k:bool(v) for k,v in checks.items()},seconds=time.monotonic()-tick,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
def main():
    dest=ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-thin-v2-result.json';assert not dest.exists();start=time.monotonic();out=ROOT/'build/analysis'/('mixed-mean-stress-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();records=[];stops=[]
    for k,n,a in [(0.,2,0.),(0.,2,.2),(.4,2,0.),(.4,2,.2),(.4,4,.2)]:
        try:
            r=case(k,n,a,out);records.append(r);print('GLOBAL',k,n,a,r['global_q0_ranks'],r['error'],r['checks'],flush=True)
            if not all(r['checks'].values()):raise RuntimeError('MMS1_numeric_gate')
            if time.monotonic()-start>1800:raise RuntimeError('MMS1_series_resource_stop')
        except Exception as exc:stops.append(dict(kappa=k,n=n,amplitude=a,reason=repr(exc)));print('STOP',repr(exc),flush=True);break
        (out/'checkpoint.json').write_text(json.dumps(dict(records=records,stops=stops),indent=2),encoding='utf-8')
    complete=len(records)==5 and not stops;improvement=False
    if complete:
        a,b=records[-2:];improvement=b['global_q0_ranks'][0]>a['global_q0_ranks'][0] and all(b['error'][j]<a['error'][j] for j in [1,2])
    ps=[Path(__file__),ROOT/'rd/structure/main-shell/mixed-mean-stress/experiments/mixed_mean_stress_thin_v2.py',ROOT/'rd/structure/main-shell/mixed-mean-stress/experiments/mixed_mean_stress_space_v2.py',ROOT/'rd/structure/main-shell/mixed-mean-stress/mixed-mean-stress-formulation.md',ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-protocol.md']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_numeric_checks=bool(complete),necessary_approximation_gate=bool(complete and improvement),main_shell_adopted=False,finite_h_checks_executed=False,geometry_conditional=True,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in ps}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
