"""Global leading bending minimization with all-position interface conditions."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time,itertools
import numpy as np
from scipy.linalg import lstsq
import moment_enriched_thin as t
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def assemble_constraints(cells,n):
    ny=sum(c['ny'] for c in cells);registry={};offset=0;trace_count=0;fixed=[];equations=[];connections=[]
    for cell_index,cell in enumerate(cells):
        cell['yoffset']=offset;offset+=cell['ny'];nodes=[];transforms=[]
        for face,key in zip(cell['faces'],t.fit.edge_keys(cell_index,n)):
            Mc=face['Mcoeff'];nd=3*face['n'];sl=slice(3*face['offset'],3*(face['offset']+face['n']))
            if key not in registry:
                canonical=Mc.copy();nodeids=np.arange(trace_count,trace_count+face['n']);trace_count+=face['n'];registry[key]=dict(M=canonical,nodes=nodeids)
                if key[0]==0 and key[1]==0:fixed.extend(nodeids)
            r=registry[key];transform=Mc.T@r['M'];Tc=np.kron(transform,np.eye(3));nodes.extend(r['nodes']);transforms.append(Tc)
            if np.linalg.norm(transform.T@transform-np.eye(face['n']))>1e-10:raise ValueError('thin_face_space')
            local=np.vstack((Tc.T@cell['Z'][sl],Tc.T@cell['Z'][cell['nc']+sl.start:cell['nc']+sl.stop]))
            if 'seen' in r:
                prior,block=r['seen'];G=np.zeros((2*nd,ny));G[:,cell['yoffset']:cell['yoffset']+cell['ny']]=local;G[:,prior:prior+block.shape[1]]-=block;equations.append(G)
            else:r['seen']=(cell['yoffset'],local)
            if key[0]==0 and key[1]==0:
                G=np.zeros((2*nd,ny));G[:,cell['yoffset']:cell['yoffset']+cell['ny']]=local;equations.append(G)
        cell['trace_ids']=np.array(nodes);cell['trace_transform']=np.zeros((cell['nb'],cell['nb']));p=0
        for block in transforms:cell['trace_transform'][p:p+len(block),p:p+len(block)]=block;p+=len(block)
    G=np.vstack(equations);Z,info=t.nullspace(G);del G
    free=np.setdiff1d(np.arange(trace_count),fixed);free=free[free%3==2];return Z,info,free,trace_count

def observe(cells,values,kappa,independent=False):
    err=np.zeros(2);ref=np.zeros(2);max_constraint=0.;vertical=np.zeros(2);elastic=np.zeros(2)
    for cell,q in zip(cells,values):
        op=t.operators(cell['q'],15 if independent else 13,9 if independent else 7,independent);s=op['state'];roots,b=t.leading_rows(op,kappa)
        residual=roots[0]@q[0]+roots[1]@q[1]+roots[2]@q[2]-b;nr=len(op['p'])*6
        elastic+=np.sum(residual[:nr]**2,axis=0);vertical+=np.sum(residual[nr:]**2,axis=0);ref+=np.sum(b*b,axis=0)
        C=t.constraints(op);v=np.vstack((q[0],q[1]));max_constraint=max(max_constraint,float(np.linalg.norm(C@v)/max(1.,np.linalg.norm(C)*np.linalg.norm(v))))
    return dict(error=np.sqrt((elastic+vertical)/ref).tolist(),elastic_error=np.sqrt(elastic/ref).tolist(),vertical_error=np.sqrt(vertical/ref).tolist(),constraint_residual=max_constraint)

def case(k,n,amp,out):
    start=time.monotonic();grid=t.m.mesh(k,n,amp);cells=[]
    for index,ids in enumerate(grid['cells']):
        cells.append(t.local(grid['ref'][ids].reshape(-1),k))
        print('local',k,n,amp,index,'ny',cells[-1]['ny'],'q0dim',cells[-1]['info']['certificate']['exact_membrane_nullity'],flush=True)
        if time.monotonic()-start>900 or memory.memory_bytes()>2*1024**3:raise RuntimeError('thin_local_resource_limit')
    Zy,info,free,count=assemble_constraints(cells,n);ny=Zy.shape[1];lookup={int(node):i for i,node in enumerate(free)};rows=sum(len(c['R']) for c in cells);A=np.zeros((rows,ny+3*len(free)));b=np.zeros((rows,2));offset=0;reference=np.zeros(2)
    for cell in cells:
        R=cell['R'];yc=cell['ny'];nb=cell['nb'];sl=slice(offset,offset+len(R));offset+=len(R);ys=slice(cell['yoffset'],cell['yoffset']+yc)
        A[sl,:ny]=R[:,:yc]@Zy[ys];tail=R[:,yc:yc+nb]@cell['trace_transform']
        for j,node in enumerate(cell['trace_ids']):
            if int(node) in lookup:A[sl,ny+3*lookup[int(node)]:ny+3*lookup[int(node)]+3]=tail[:,3*j:3*j+3]
        b[sl]=R[:,yc+nb:yc+nb+2];reference+=cell['reference']
    from scipy.linalg import svd
    Ay=A[:,:ny];Ab=A[:,ny:];U,sb,Vb=svd(Ab,full_matrices=False,check_finite=False);rb=t.rank_info(sb)
    if len(set(rb))!=1:raise ValueError('global_q2_rank_ambiguous:'+str(rb))
    Qb=U[:,:rb[0]];Ar=Ay-Qb@(Qb.T@Ay);br=b-Qb@(Qb.T@b)
    y,_,rank,sv=lstsq(Ar,br,cond=1e-10,lapack_driver='gelsd',check_finite=False);ranks=t.rank_info(sv)
    if len(set(ranks))!=1:raise ValueError('leading_objective_rank_ambiguous:'+str(ranks))
    xb=(Vb[:rb[0]].T/sb[:rb[0]])@(Qb.T@(b-Ay@y));x=np.vstack((y,xb));del Ar,br,U,Qb,Vb,Ay,Ab
    residual=A@x-b;fit_error=np.sqrt(np.sum(residual*residual,axis=0)/reference);kkt=np.linalg.norm(A.T@residual,axis=0)/np.maximum(1e-30,np.linalg.norm(A)*np.linalg.norm(b,axis=0));del A,b
    Y=Zy@x[:ny];del Zy;trace=np.zeros((count,3,2));trace[free]=x[ny:].reshape(-1,3,2);values=[]
    for cell in cells:
        y=Y[cell['yoffset']:cell['yoffset']+cell['ny']];q2b=cell['trace_transform']@trace[cell['trace_ids']].reshape(-1,2);values.append(t.recover(cell,y,q2b,k))
    primary=observe(cells,values,k);ind=observe(cells,values,k,True)
    checks=dict(interface_rank=len(set(info['ranks']))==1,local_constraint=max(c['info']['residual'] for c in cells)<=1e-9,independent_constraint=ind['constraint_residual']<=1e-9,
        objective_rank=len(set(ranks))==1,kkt=max(kkt)<=1e-8,projection=max(abs(a-b) for a,b in zip(primary['error'],fit_error))<=1e-8,
        independent=max(abs(a-b) for a,b in zip(primary['error'],ind['error']))<=1e-8,exact_certificates=all(c['info']['certificate']['all_rational_substitutions_zero'] for c in cells),
        geometry_match=max(c['info']['certificate']['max_rounded_reference_difference'] for c in cells)<=1e-12,resources=time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3)
    path=out/f'k{k:g}-n{n}-a{amp:g}.npz';np.savez_compressed(path,values=np.array(values),ref=grid['ref'],cells=grid['cells'],trace_q2=trace)
    return dict(kappa=k,n=n,amplitude=amp,local_dimensions=[c['info']['nullity'] for c in cells],local_isometry_certificates=[c['info']['certificate'] for c in cells],interface=info,
        objective_ranks=ranks,global_q2_ranks=rb,q01_reduced_columns=ny,objective_columns=len(x),fit_error=fit_error.tolist(),primary=primary,independent=ind,kkt=kkt.tolist(),checks={k:bool(v) for k,v in checks.items()},
        geometry_conditional=True,not_a_finite_h_lower_bound=True,seconds=time.monotonic()-start,artifact=path.relative_to(ROOT).as_posix(),artifact_sha256=digest(path))

def main():
    dest=ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-thin-global-result.json';assert not dest.exists();records=[];stops=[];start=time.monotonic();out=ROOT/'build/analysis'/('moment-enriched-thin-global-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    assert json.loads((ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-thin-local-result.json').read_text())['all_checks_pass']
    for k,n,amp in [(0.,2,0.),(0.,2,.2),(.4,2,0.),(.4,2,.2),(.4,4,.2)]:
        try:
            r=case(k,n,amp,out);records.append(r);print('GLOBAL',k,n,amp,r['primary'],r['checks'],r['seconds'],flush=True)
            if not all(r['checks'].values()):stops.append(dict(kappa=k,n=n,amplitude=amp,reason='leading_gate_failed'));break
        except Exception as exc:stops.append(dict(kappa=k,n=n,amplitude=amp,reason=repr(exc)));print('STOP',repr(exc),flush=True);break
        (out/'checkpoint.json').write_text(json.dumps(dict(records=records,stops=stops),indent=2),encoding='utf-8')
        if time.monotonic()-start>1800:stops.append(dict(reason='analysis_series_resource_limit'));break
    paths=[Path(__file__),ROOT/'rd/structure/main-shell/moment-enriched/experiments/moment_enriched_thin.py',ROOT/'rd/structure/main-shell/moment-enriched/experiments/moment_enriched_exact_isometry.py',ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-thin-protocol.md',ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-thin-preparation.json']
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_analysis_checks=len(records)==5 and not stops,main_shell_adopted=False,geometry_conditional=True,physical_acceptance_thresholds_unchanged=True,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
