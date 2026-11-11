"""Energy-minimizing MMS-1 leading q0/q1/q2, five frozen configurations."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.linalg import svd,lstsq
import mixed_mean_stress_bending as t
import check_moment_enriched_thin_global_v2 as shared
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def case(k,n,amp,out):
    start=time.monotonic();grid=t.m.mesh(k,n,amp);cells=[]
    for index,ids in enumerate(grid['cells']):
        cell=t.local(grid['ref'][ids].reshape(-1),k);cells.append(cell);print('LOCAL',k,n,amp,index,'q01',cell['ny'],'q2',cell['q2_internal_ranks'],flush=True)
        if time.monotonic()-start>900 or memory.memory_bytes()>2*1024**3:raise RuntimeError('MMS1_bending_local_resources')
    shared.t=t;Zy,interface,free,count=shared.assemble_constraints(cells,n);ny=Zy.shape[1];lookup={int(v):i for i,v in enumerate(free)};rows=sum(len(c['R']) for c in cells);A=np.zeros((rows,ny+3*len(free)));b=np.zeros((rows,3));offset=0;reference=np.zeros(3)
    for cell in cells:
        R=cell['R'];yc=cell['ny'];nb=cell['nb'];sl=slice(offset,offset+len(R));offset+=len(R);ys=slice(cell['yoffset'],cell['yoffset']+yc);A[sl,:ny]=R[:,:yc]@Zy[ys];tail=R[:,yc:yc+nb]@cell['trace_transform']
        for j,node in enumerate(cell['trace_ids']):
            if int(node) in lookup:A[sl,ny+3*lookup[int(node)]:ny+3*lookup[int(node)]+3]=tail[:,3*j:3*j+3]
        b[sl]=R[:,yc+nb:yc+nb+3];reference+=cell['reference']
    Ay=A[:,:ny];Ab=A[:,ny:];U,sb,Vb=svd(Ab,full_matrices=False,check_finite=False);rb=t.rank_info(sb)
    if len(set(rb))!=1:raise ValueError('MMS1_bending_global_q2_rank:'+str(rb))
    Qb=U[:,:rb[0]];Ar=Ay-Qb@(Qb.T@Ay);br=b-Qb@(Qb.T@b);y,_,rank,sv=lstsq(Ar,br,cond=1e-10,lapack_driver='gelsd',check_finite=False);ranks=t.rank_info(sv)
    if ranks!=[ny]*3:
        path=out/f'k{k:g}-n{n}-a{amp:g}-stopped-objective.npz';np.savez_compressed(path,A=Ar,b=br,singular=sv,q01_map=Zy);(path.with_suffix('.json')).write_text(json.dumps(dict(ranks=ranks,columns=ny,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path)),indent=2),encoding='utf-8');raise ValueError('MMS1_leading_objective_not_full_rank:'+str(ranks)+'/'+str(ny))
    xb=(Vb[:rb[0]].T/sb[:rb[0]])@(Qb.T@(b-Ay@y));x=np.vstack((y,xb));residual=A@x-b;fit_error=np.sqrt(np.sum(residual*residual,axis=0)/reference);kkt=np.linalg.norm(A.T@residual,axis=0)/np.maximum(1e-30,np.linalg.norm(A)*np.linalg.norm(b,axis=0));del A,b,Ay,Ab,Ar,br,U,Qb,Vb
    Y=Zy@x[:ny];trace=np.zeros((count,3,3));trace[free]=x[ny:].reshape(-1,3,3);values=[]
    for cell in cells:
        localy=Y[cell['yoffset']:cell['yoffset']+cell['ny']];q2b=cell['trace_transform']@trace[cell['trace_ids']].reshape(-1,3);values.append(t.recover(cell,localy,q2b,k))
    primary=t.observe(cells,values,k);ind=t.observe(cells,values,k,True);geometry=max(t.geometry_difference(c['q'],k) for c in cells)
    checks=dict(interface=len(set(interface['ranks']))==1,local=max(c['residual'] for c in cells)<=1e-9,constraints=max(ind['constraints'].values())<=1e-9,objective_full_rank=ranks==[ny]*3,kkt=max(kkt)<=1e-8,projection=max(abs(a-b) for a,b in zip(primary['error'],fit_error))<=1e-8,independent_energy=max(abs(a-b) for a,b in zip(primary['error'],ind['error']))<=1e-8,independent_position=max(abs(a-b) for a,b in zip(primary['position_error'],ind['position_error']))<=1e-8,geometry=geometry<=1e-12,stress_certificate=all(c['stress']['exact_space_certificate']['exact_lower_rank']==c['stress']['exact_space_certificate']['upper_rank'] for c in cells),resources=time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3)
    path=out/f'k{k:g}-n{n}-a{amp:g}.npz';np.savez_compressed(path,values=np.array(values),ref=grid['ref'],cells=grid['cells'],trace_q2=trace)
    return dict(kappa=k,n=n,amplitude=amp,primary=primary,independent=ind,fit_error=fit_error.tolist(),kkt=kkt.tolist(),q01_columns=ny,objective_ranks=ranks,global_q2_ranks=rb,interface=interface,local_q2_ranks=[c['q2_internal_ranks'] for c in cells],geometry_difference=geometry,checks={k:bool(v) for k,v in checks.items()},seconds=time.monotonic()-start,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))
def main():
    dest=ROOT/'docs/verification/mixed-mean-stress-bending-result.json';assert not dest.exists();assert json.loads((ROOT/'docs/verification/mixed-mean-stress-bending-coefficient-result.json').read_text())['all_checks_pass'];start=time.monotonic();records=[];stops=[];out=ROOT/'build/analysis'/('mixed-mean-stress-bending-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    for k,n,amp in [(0.,2,0.),(0.,2,.2),(.4,2,0.),(.4,2,.2),(.4,4,.2)]:
        try:
            r=case(k,n,amp,out);records.append(r);print('BENDING',k,n,amp,r['primary'],r['checks'],r['seconds'],flush=True)
            if not all(r['checks'].values()):raise RuntimeError('MMS1_bending_numerical_gate')
            if time.monotonic()-start>1800:raise RuntimeError('MMS1_bending_series_resources')
        except Exception as exc:stops.append(dict(kappa=k,n=n,amplitude=amp,reason=repr(exc)));print('STOP',repr(exc),flush=True);break
        (out/'checkpoint.json').write_text(json.dumps(dict(records=records,stops=stops),indent=2),encoding='utf-8')
    complete=len(records)==5 and not stops;accuracy=False
    if complete:
        a,b=records[-2:];accuracy=all(r['primary']['error'][0]<=1e-8 and r['primary']['position_error'][0]<=1e-8 for r in records) and all(b['primary'][field][j]<a['primary'][field][j] and b['primary'][field][j]<=.01 for field in ['error','position_error'] for j in [1,2])
    ps=[Path(__file__),Path(t.__file__),ROOT/'tools/design-experiments/mixed_mean_stress_thin_v6.py',ROOT/'docs/verification/mixed-mean-stress-bending-protocol.md',ROOT/'docs/verification/mixed-mean-stress-bending-coefficient-result.json']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_numerical_checks=bool(complete),necessary_accuracy_gate=bool(complete and accuracy),main_shell_adopted=False,finite_h_global_executed=False,geometry_conditional=True,kinds=t.KINDS,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),sources={p.relative_to(ROOT).as_posix():digest(p) for p in ps}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
