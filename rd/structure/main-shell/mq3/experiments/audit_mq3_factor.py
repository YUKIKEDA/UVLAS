"""General-LU complex-thickness and tensor-energy evaluation of saved solves."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,itertools,time
import numpy as np
from scipy.linalg import solve
import mq3_factor_family as family
ROOT=rd_paths.ROOT
LAM=.3/(1.3*.4);MU=1/2.6
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def energy(E):return LAM*np.trace(E,axis1=1,axis2=2)**2+2*MU*np.einsum('pij,pij->p',E,E)
def main():
    dest=ROOT/'rd/structure/main-shell/mq3/records/mq3-factor-audit.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();parent=ROOT/'rd/structure/main-shell/mq3/records/mq3-factor-result.json';data=json.loads(parent.read_text());records=[];cache={};stop=None;p,w=np.polynomial.legendre.leggauss(11);plane=np.array(list(itertools.product(p,p)));points=np.column_stack((plane,np.ones(len(plane))));weights=np.array([a*b for a,b in itertools.product(w,w)]);N,g=family.base.geometry.shape(plane)
    for source,sha in data['sources'].items():assert digest(ROOT/source)==sha
    try:
        for record in data['records']:
            family.configure(record['metric'],record['work'])
            if not all(record['checks'].values()):continue
            path=ROOT/record['artifact'];assert digest(path)==record['artifact_sha256'];saved=np.load(path);free=saved['free'];lookup={int(v):i for i,v in enumerate(free)};totals={kind:dict(error=0.,reference=0.,zero=0.,h1=0.) for kind in ['quadratic','quartic']};kernel_square=matrix_square=area=0.
            for ids,key in zip(saved['ids'],record['map_keys']):
                if key not in cache:
                    mp=ROOT/data['maps'][key]['artifact'];assert digest(mp)==data['maps'][key]['artifact_sha256'];m=np.load(mp);state={name:m[name] for name in ['q','T','rows','U']};S,b,_,_=family.system(state,0.);X0=solve(S,b,assume_a='gen',check_finite=False);Sc,bc,_,_=family.system(state,1e-25j);Xc=solve(Sc,bc,assume_a='gen',check_finite=False);B0=family.evaluate(state,0.,X0,points)[:,:,:144];B1=2*family.evaluate(state,1e-25j,Xc,points).imag[:,:,:144]/1e-25;J,_=family.raw(state['q'],0.,points);cache[key]=(B0,B1,J)
                B0,B1,J=cache[key];weight=weights*np.linalg.det(J);K=np.zeros((144,saved['kernel'].shape[1]))
                for local,node in enumerate(ids):
                    if int(node) in lookup:K[9*local:9*local+6]=saved['kernel'][6*lookup[int(node)]:6*lookup[int(node)]+6]
                kernel_square+=float(np.sum(weight[:,None,None]*(B0@K)**2));matrix_square+=float(np.sum(weight[:,None,None]*B0[:,:,np.array([9*a+b for a in range(16) for b in range(6)])]**2));x=(N@saved['coords'][ids,0])[:,0];k=record['kappa'];avec=np.column_stack((np.ones_like(x),np.zeros_like(x),k*x));D=1+k*k*x*x;Jsurf=np.einsum('na,pnb->pab',saved['coords'][ids,0],g);G=Jsurf.transpose(0,2,1)@Jsurf;aw=weights*np.sqrt(np.linalg.det(G));area+=aw.sum()
                for kind in totals:
                    u=saved['u_'+kind][ids];v=np.ones_like(x) if kind=='quadratic' else x*x;target=-v[:,None,None]*avec[:,:,None]*avec[:,None,:]/D[:,None,None];E=np.einsum('kij,pk->pij',family.base.MB,B1@u.reshape(144));totals[kind]['error']+=weight@energy(E-target);totals[kind]['reference']+=weight@energy(target);totals[kind]['zero']+=float(np.sum(weight[:,None]*(B0@u.reshape(144))**2));val=N@u[:,3:9];du=np.einsum('na,pnb->pab',u[:,3:9],g);totals[kind]['h1']+=aw@(np.sum(val*val,axis=1)+np.einsum('pai,pij,paj->p',du,np.linalg.inv(G),du))
            kernel=float(np.sqrt(kernel_square/matrix_square)/np.linalg.norm(saved['kernel']));observations={};checks=dict(kernel=kernel<=1e-9)
            for kind,value in totals.items():
                error=float(np.sqrt(value['error']/value['reference']));h1=float(np.sqrt(value['h1']/area));rdnorm=np.linalg.norm(saved['u_'+kind][free,:6]);zero=float(np.sqrt(value['zero']/matrix_square)/rdnorm);difference=abs(error-record['solutions'][kind]['observed_error']);hdiff=abs(h1/record['solutions'][kind]['observed_h1']-1);good=difference<=1e-8 and hdiff<=1e-8 and zero<=1e-9;checks[kind]=good;observations[kind]=dict(error=error,error_difference=difference,h1_difference=hdiff,zero_residual=zero)
            records.append(dict(metric=record['metric'],work=record['work'],kappa=record['kappa'],amplitude=record['amplitude'],n=record['n'],kernel_residual=kernel,observations=observations,checks={key:bool(value) for key,value in checks.items()}));print(record['kappa'],record['amplitude'],record['n'],checks,flush=True)
            if not all(checks.values()):stop='required_audit_gate_failed';break
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(family.__file__),Path(family.base.geometry.__file__),parent,ROOT/'rd/structure/main-shell/mq3/records/mq3-factor-audit-protocol.md'];complete=sum(all(r['checks'].values()) for r in data['records']);dest.write_text(json.dumps(dict(records=records,validated_primary_configurations=complete,all_checks_pass=len(records)==complete and stop is None,primary_all_configurations_complete=data['all_checks_pass'],stop_reason=stop,independent_global_resolve=False,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
if __name__=='__main__':main()
