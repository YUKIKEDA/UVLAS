"""Independent Green-field and full 81-coordinate energy differentiation audit."""
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
import q9_edge_moment_fields as field
import audit_q9_curved_strain_screen as independent
import check_q9_curved_load as geometry
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def relative(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-30))
def main():
    target=D/'q9-edge-moment-field-audit.json'
    if target.exists():raise FileExistsError(target)
    parent=D/'q9-edge-moment-field-result.json';record=json.loads(parent.read_text(encoding='utf-8'))
    assert record['all_criteria_pass'];assert record['source_sha256']==digest(ROOT/'tools/design-experiments/check_q9_edge_moment_fields.py')
    assert record['protocol_sha256']==digest(D/'q9-edge-moment-field-experiment.md')
    assert record['contract_sha256']==digest(ROOT/'docs/q9-edge-moment-field-interface.md')
    for name,value in record['dependencies'].items():assert digest(ROOT/'tools/design-experiments'/name)==value
    start=time.monotonic();cases=[];stop=None
    M=independent.E.copy();M[3:]*=np.sqrt(2.)
    A=np.arange(36,dtype=float).reshape(6,6)/40;C=A.T@A+np.diag([1.,2.,3.,4.,5.,6.])
    try:
        for name,h in itertools.product(['RECT','SKEW','TRAP','SADDLE','NORMAL-C'],[.1,.01]):
            q=geometry.geometry(name);ref=field.ReferenceReduction(q,h)
            points=[ref.points,ref.ties]+ref.edge_points
            samples=[independent.data(q,h,p) for p in points]
            maps=[v[2] for v in samples];U=[]
            for (axis,_),(_,J,_) in zip(field.EDGES,samples[2:]):
                Ae=J[:,:,[axis,2]];L=np.linalg.cholesky(Ae.transpose(0,2,1)@Ae)
                U.append(np.linalg.solve(L,Ae.transpose(0,2,1)).transpose(0,2,1))
            def inputs(trial):
                fields=[]
                for g in maps:
                    F=np.einsum('ni,pnj->pij',trial.reshape(27,3),g)
                    fields.append((F.transpose(0,2,1)@F-np.eye(3))/2)
                v,t=[np.einsum('aij,pij->pa',M,e) for e in fields[:2]]
                e=np.array([np.stack((np.einsum('pi,pij,pj->p',u[:,:,0],x,u[:,:,0]),np.einsum('pi,pij,pj->p',u[:,:,1],x,u[:,:,1]),np.sqrt(2.)*np.einsum('pi,pij,pj->p',u[:,:,0],x,u[:,:,1])),axis=1) for u,x in zip(U,fields[2:])])
                return v,t,e
            F0=np.array([[1.04,.03,-.02],[.01,.97,.02],[0.,-.01,1.02]])
            natural=(q.reshape(27,3)@F0.T).reshape(81);initial=inputs(natural)
            trial=natural.copy();qc=trial.reshape(9,3,3);qc[:,0,2]+=.02*np.sin(np.arange(9));qc[:,1,0]+=.015*np.cos(np.arange(9));qc[:,2,1]+=.01*np.sin(2*np.arange(9))
            primary=field.green_sources(ref,trial);other=inputs(trial)
            difference=max(relative(a,b) for a,b in zip(primary,other))
            def energy(x):
                elastic=ref.kernel(tuple(a-b for a,b in zip(inputs(x),initial)))
                return np.einsum('p,pi,ij,pj->',ref.measure,elastic,C,elastic)/2
            _,force=field.energy_force(ref,trial,initial,C)
            grad=[];K=[];Kcs=[]
            for basis in np.eye(81):
                perturbed=trial.astype(complex)+1e-25j*basis
                grad.append(float(energy(perturbed).imag/1e-25))
                K.append(field.tangent_action(ref,trial,initial,C,basis)[0])
                Kcs.append(field.energy_force(ref,perturbed,initial,C)[1].imag/1e-25)
            K=np.array(K).T;Kcs=np.array(Kcs).T
            min_det=min(float(np.min(np.linalg.det(np.einsum('ni,pnj->pij',trial.reshape(27,3),g)))) for g in maps)
            force_error=relative(np.array(grad),force);tangent_error=relative(K,Kcs);symmetry=relative(K,K.T)
            checks=dict(field=difference<=1e-10,force=force_error<=1e-9,tangent=tangent_error<=1e-9,symmetry=symmetry<=1e-10,orientation=min_det>0)
            cases.append(dict(geometry=name,h=h,field=difference,force=force_error,tangent=tangent_error,symmetry=symmetry,minimum_sample_detF=min_det,checks={k:bool(v) for k,v in checks.items()}))
            print(json.dumps(cases[-1]),flush=True)
            if time.monotonic()-start>300 or memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
    except Exception as exc:stop=repr(exc)
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),protocol_sha256=digest(D/'q9-edge-moment-field-audit-protocol.md'),dependencies={Path(m.__file__).name:digest(m.__file__) for m in [field,independent,geometry,memory]},cases=cases,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=bool(stop is None and len(cases)==10 and all(all(c['checks'].values()) for c in cases)),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['cases','dependencies']}))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
