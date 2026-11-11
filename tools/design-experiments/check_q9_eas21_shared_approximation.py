"""Selected global strain best approximation with EAS energy condensation."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import q9_eas21 as method
import q9_developable_bending_approximation as base
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def assemble(n,h):
    xy,cells,coords,u,free=base.mesh(.4,n,.2);lookup={int(v):i for i,v in enumerate(free)};nr=len(free)*3;np_=len(free)*6
    route=dict(W=np.zeros((82*len(cells),nr)),D=np.zeros((82*len(cells),np_)),e=np.zeros(82*len(cells)))
    points,w,N,g=base.table(9);NS,gs,ws=base.surface_table(5);M=np.zeros((len(xy),len(xy)));area=0.;total=0.;states=[];interp=0.;compression=0.
    for cell,(_,_,ids,_) in enumerate(cells):
        q=coords[ids].reshape(81);state=method.prepare(q,h);states.append(state);J,_,B,_=method.evaluate(state,points);weight=w*np.linalg.det(J)
        star=base.target(.4,(N@xy[ids])[:,0],h*points[:,2]/2);e=(star@method.CH.T*np.sqrt(weight)[:,None]).reshape(-1);total+=e@e
        A=(np.einsum('ij,pjk->pik',method.CH,B)*np.sqrt(weight)[:,None,None]).reshape(-1,81);R,log=base.compress(np.column_stack((A,e)));compression=max(compression,log['error_upper']);sl=slice(82*cell,82*(cell+1));route['e'][sl]=R[:,81]
        for local,node in enumerate(ids):
            if int(node) not in lookup:continue
            j=lookup[int(node)];route['W'][sl,3*j:3*j+3]=R[:,9*local:9*local+3];route['D'][sl,6*j:6*j+6]=R[:,9*local+3:9*local+9]/np.repeat([1.,base.core.LENGTH],3)[None,:]
        residual=A@u[ids].reshape(81)-e;interp+=residual@residual
        surfaceJ=np.einsum('ia,pib->pab',coords[ids,0,:],gs);metric=surfaceJ.transpose(0,2,1)@surfaceJ;aw=ws*np.sqrt(np.linalg.det(metric));area+=aw.sum();localH=np.einsum('p,pi,pj->ij',aw,NS,NS)+base.core.LENGTH**2*np.einsum('p,pia,pab,pjb->ij',aw,gs,np.linalg.inv(metric),gs);M[np.ix_(ids,ids)]+=localH
    M=M[np.ix_(free,free)]/area;pref=u[free,3:9].copy();pref[:,3:]*=base.core.LENGTH
    return dict(routes={'eas':route},M=M,xy=xy,coords=coords,free=free,cells=cells,reference=u,reference_p=pref,reference_norm=float(total),interpolation=float(np.sqrt(interp/total)),compression= compression),states
def observe(a,states,h,u):
    points,w,N,g=base.table(11);total=error=energy=upper2=0.
    for state,(_,_,ids,_) in zip(states,a['cells']):
        J,_,B,_=method.evaluate(state,points);weight=w*np.linalg.det(J);star=base.target(.4,(N@a['xy'][ids])[:,0],h*points[:,2]/2);e=(star@method.CH.T*np.sqrt(weight)[:,None]).reshape(-1)
        A=(np.einsum('ij,pjk->pik',method.CH,B)*np.sqrt(weight)[:,None,None]).reshape(-1,81);value,ve=base.rounding.dot_rows(A,u[ids].reshape(81));delta=value-e;_,high=base.rounding.norms(delta,ve+base.rounding.gamma(3)*(abs(value)+abs(e)));total+=e@e;error+=delta@delta;energy+=value@value;upper2+=high*high
    margin=base.rounding.gamma(len(points)*6*len(a['cells'])+10)
    return dict(error=float(np.sqrt(error/total)),error_upper=float(np.sqrt(upper2*(1+margin)/(total*(1-margin)))),energy_ratio=float(energy/total),reference_norm=float(total))
def main():
    dest=ROOT/'docs/verification/q9-eas21-shared-approximation-result.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();out=ROOT/'build/analysis'/('q9-eas21-shared-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();records=[];stop=None
    try:
        for n,h in [(4,.01),(4,.0001),(8,.001),(8,.0001)]:
            a,states=assemble(n,h);R,qr,solutions,condition=base.solve_route(a,'eas',[2.]);fit=solutions[0];r=fit.pop('r');p=fit.pop('p');u=base.full_field(a,r,p);obs=observe(a,states,h,u);position=base.position_error(.4,a,u,11);artifact=out/f'n{n}-h{h:g}.npz';np.savez_compressed(artifact,u=u,r=r,p=p)
            checks=dict(kkt=max(fit['kkt_r'],fit['kkt_p'])<=1e-8,quadrature=abs(fit['full_error']-obs['error'])<=1e-8,sphere=base.core.pnorm(p.reshape(-1,6),a['M'])<=2*fit['reference_h1']*(1+1e-8),resources=memory.memory_bytes()<4*1024**3 and time.monotonic()-start<600)
            records.append(dict(n=n,h=h,fit=fit,observation=obs,interpolation=a['interpolation'],position_error=position,local_compression_upper=a['compression'],qr=qr,position_condition=condition['position_condition'],checks=checks,artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact)));print(n,h,obs['error'],position,checks,flush=True)
            if not all(checks.values()):stop='numerical_gate';break
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(method.__file__),Path(base.__file__),ROOT/'docs/verification/q9-eas21-shared-approximation-experiment.md',ROOT/'docs/verification/q9-eas21-feasibility-result.json']
    gate=None
    if len(records)==4 and stop is None:
        amplification=records[1]['observation']['error']/records[0]['observation']['error'];gate=dict(absolute_1pct=all(r['observation']['error_upper']<=.01 for r in records[2:]),thin_amplification=amplification,thin_uniform=amplification<=2.)
    dest.write_text(json.dumps(dict(candidate='EAS21-linear-homogeneous',physical_adoption=False,records=records,gate=gate,stop_reason=stop,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
