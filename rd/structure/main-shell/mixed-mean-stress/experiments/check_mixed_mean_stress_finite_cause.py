"""R1: original-point energies and assembled component virtual work; no solve."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import gc,hashlib,json,time
import numpy as np
import mixed_mean_stress_finite as fit
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT
LABELS=['thickness_variation','projected_mean','top_bottom','vertical']
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def evaluate_case(record):
    tick=time.monotonic();path=ROOT/record['artifact'];assert digest(path)==record['artifact_sha256']
    with np.load(path) as data:values=data['local'].copy();ref=data['ref'].copy();cells=data['cells'].copy();trace=data['trace'].copy()
    h=record['h'];k=record['kappa'];n=record['n'];energy=np.zeros((4,3));refnorm=np.zeros(3);identity=np.zeros(3)
    registry={};trace_ref=[];fixed=[];rows=[];work=np.zeros((2,4,3));gwork=np.zeros_like(work);cell_times=[]
    for index,(ids,u) in enumerate(zip(cells,values)):
        start=time.monotonic();q=ref[ids].reshape(-1);s=fit.f.prepare(q,h,k)
        nodes,scalarT=fit.topology.register(s,fit.topology.edge_keys(index,n),registry,trace_ref,fixed);T=np.kron(scalarT,np.eye(3));ni=3*s['ni'];nc=u.shape[0]
        star,_=fit.reference(q,s['p'],h,k);e=s['B']@u-star;stress,high,low=fit.f.project(s,e)
        Sh=np.einsum('ij,pjt->pit',fit.m.C,high,optimize=True);Sm=stress-Sh
        energy[0]+=np.einsum('p,pkt,pkt->t',s['weight'],high,Sh);energy[1]+=np.sum(low*low,axis=0)
        identity+=np.einsum('p,pkt,pkt->t',s['weight'],e,stress)
        sr=fit.strain_rows(s,star);refnorm+=np.sum(sr*sr,axis=0)
        gradients=np.zeros((4,nc,3));gradients[0]=np.einsum('p,pkn,pkt->nt',s['weight'],s['B'],Sh,optimize=True)
        gradients[1]=np.einsum('p,pkn,pkt->nt',s['weight'],s['B'],Sm,optimize=True)
        direct=np.einsum('p,pkn,pkt->nt',s['weight'],s['B'],stress,optimize=True)
        face_info=[]
        for face in s['faces']:
            group=2 if face['axis']==2 else 3;jump=face['jump']@u.reshape(s['nscalar'],-1);jump=jump.reshape(len(jump),3,3);ww=s['tau']*face['measure']
            term=np.einsum('p,pat,pat->t',ww,jump,jump);energy[group]+=term;identity+=term
            g=np.einsum('p,pn,pat->nat',ww,face['jump'],jump,optimize=True).reshape(nc,3);gradients[group]+=g;direct+=g
            face_info.append((group,face,jump,ww))
        W=s['root']/h**1.5;norm2=np.sum(W*W,axis=0);canonical_norm=np.sum((W[:,ni:]@T)**2,axis=0).reshape(-1,3)
        for direction in range(2):
            phase=.137 if direction==0 else .271
            interior=np.sin(phase*(1+index*ni+np.arange(ni))[:,None]+.37*np.arange(3)[None,:])
            bgid=100000+3*nodes[:,None]+np.arange(3)[None,:]
            boundary=np.sin(phase*(1+bgid).reshape(-1,1)+.37*np.arange(3)[None,:])
            boundary[np.repeat(np.isin(nodes,fixed),3)]=0
            d=h**1.5*np.vstack((interior,T@boundary));de=s['B']@d;_,dh,dl=fit.f.project(s,de)
            work[direction,0]+=np.einsum('p,pkt,pkt->t',s['weight'],dh,Sh)
            work[direction,1]+=np.sum(dl*low,axis=0)
            for group,face,jump,ww in face_info:
                dj=(face['jump']@d.reshape(s['nscalar'],-1)).reshape(len(jump),3,3)
                work[direction,group]+=np.einsum('p,pat,pat->t',ww,dj,jump)
            gwork[direction]+=np.einsum('cnt,nt->ct',gradients,d)
        rows.append(dict(nodes=nodes,T=T,ni=ni,gradients=gradients/h**3,direct=direct/h**3,norm2=norm2,canonical_norm=canonical_norm))
        cell_times.append(time.monotonic()-start)
        del s,W;gc.collect()
        if index%4==3:print('CELL',record['amplitude'],h,index+1,'seconds',time.monotonic()-tick,flush=True)
        if time.monotonic()-tick>300 or memory.memory_bytes()>2*1024**3:raise RuntimeError('R1_resources')
    free=np.setdiff1d(np.arange(len(trace_ref)),fixed);trace_g=np.zeros((4,len(trace_ref),3,3));direct_g=np.zeros((len(trace_ref),3,3));trace_norm=np.zeros((len(trace_ref),3));inside=[];inside_direct=[]
    for row in rows:
        ni=row['ni'];g=row['gradients'];inside.append(g[:,:ni]/np.sqrt(row['norm2'][:ni])[None,:,None]);inside_direct.append(row['direct'][:ni]/np.sqrt(row['norm2'][:ni,None]))
        for group in range(4):trace_g[group,row['nodes']]+=(row['T'].T@g[group,ni:]).reshape(-1,3,3)
        direct_g[row['nodes']]+=(row['T'].T@row['direct'][ni:]).reshape(-1,3,3);trace_norm[row['nodes']]+=row['canonical_norm']
    assembled=np.concatenate(inside+[trace_g[:,free].reshape(4,-1,3)/np.sqrt(trace_norm[free].reshape(1,-1,1))],axis=1)
    assembled_direct=np.vstack(inside_direct+[direct_g[free].reshape(-1,3)/np.sqrt(trace_norm[free].reshape(-1,1))])
    count=assembled.shape[1];denominator=np.sqrt(count)*np.sqrt(refnorm/h**3)
    total=assembled.sum(axis=0);gram=np.einsum('cit,dit->tcd',assembled,assembled);norm=np.sqrt(np.maximum(0,np.einsum('tcc->tc',gram)))
    cosine=gram/np.maximum(1e-300,norm[:,:,None]*norm[:,None,:]);kkt=np.linalg.norm(total,axis=0)/denominator
    energy_error=np.sqrt(energy.sum(axis=0)/refnorm);identity_error=np.abs(identity-energy.sum(axis=0))/np.maximum(refnorm,1e-300)
    work_scale=np.maximum(1e-30,np.sum(np.abs(work),axis=1)+np.sum(np.abs(gwork),axis=1))
    work_error=np.max(np.abs(work-gwork)/work_scale[:,None,:])
    gradient_error=np.max(np.linalg.norm(total-assembled_direct,axis=0)/denominator)
    parent=np.array(record['fit']['error']);error_difference=np.max(np.abs(energy_error-parent))
    checks=dict(reconstruction=error_difference<=1e-8,constitutive_identity=max(identity_error)<=1e-8,virtual_work=work_error<=1e-8,gradient_sum=gradient_error<=1e-8,original_stationarity=max(kkt)<=1e-8,resources=time.monotonic()-tick<=300 and memory.memory_bytes()<=2*1024**3)
    return dict(amplitude=record['amplitude'],h=h,n=n,energy=energy.tolist(),energy_over_h3=(energy/h**3).tolist(),reference_norm=refnorm.tolist(),reference_over_h3=(refnorm/h**3).tolist(),normalized_energy_squared=(energy/refnorm).tolist(),energy_fractions=(energy/np.maximum(energy.sum(axis=0),1e-300)).tolist(),error=energy_error.tolist(),component_gradient_indicators=(norm/denominator[:,None]).tolist(),component_gradient_cosines=cosine.tolist(),original_kkt=kkt.tolist(),parent_kkt=record['fit']['direct_kkt'],gradient_sum_error=float(gradient_error),virtual_work_error=float(work_error),constitutive_identity_error=identity_error.tolist(),parent_error_difference=float(error_difference),virtual_work=work.tolist(),gradient_virtual_work=gwork.tolist(),position_error_from_parent=record['observed']['recovered_position_error'],checks={key:bool(value) for key,value in checks.items()},cell_seconds=cell_times,seconds=time.monotonic()-tick,artifact=record['artifact'],artifact_sha256=record['artifact_sha256'])
def main():
    source=ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-curved-result.json';parent=json.loads(source.read_text());dest=ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-cause-result.json';assert not dest.exists()
    selected=[r for r in parent['records'] if r['n']==4];assert len(selected)==3
    records=[];stops=[];start=time.monotonic()
    for record in selected:
        try:
            r=evaluate_case(record);records.append(r);print('CAUSE',r['amplitude'],r['h'],r['energy_fractions'],r['component_gradient_indicators'],r['checks'],flush=True)
            if not all(r['checks'].values()):raise RuntimeError('R1_numerical_gate')
            if time.monotonic()-start>1200:raise RuntimeError('R1_series_resources')
        except Exception as exc:stops.append(dict(amplitude=record['amplitude'],h=record['h'],reason=repr(exc)));print('STOP',repr(exc),flush=True);break
    paths=[Path(__file__),Path(fit.__file__),ROOT/'rd/structure/main-shell/mixed-mean-stress/experiments/mixed_mean_stress_element_v2.py',ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-cause-protocol.md',source]
    dest.write_text(json.dumps(dict(records=records,stops=stops,all_checks_pass=len(records)==3 and not stops,labels=LABELS,kinds=fit.KINDS,new_solve_performed=False,new_candidate_tested=False,main_shell_adopted=False,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),sources={p.relative_to(ROOT).as_posix():digest(p) for p in paths}),indent=2),encoding='utf-8')
if __name__=='__main__':main()
