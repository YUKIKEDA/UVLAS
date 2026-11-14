"""Predeclared field/dual and nonlinear differentiation gates for G-EM1."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
import q9_edge_moment_fields as field
import q9_edge_moment as old
import check_q9_curved_load as geometry
import check_q9_independent_load as material
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;DEST=ROOT/'rd/structure/main-shell/q9-edge-trace/records/q9-edge-moment-field-result.json'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def relative(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-30))
def dot(x,y):return sum(np.sum(a*b) for a,b in zip(x,y))
def source_linear(ref):
    result=[]
    for p in [ref.points,ref.ties]:result.append(old.mandel(old.source(ref.q,ref.width,ref.offset,p)[1]))
    edges=[]
    for p,U in zip(ref.edge_points,ref.edge_U):
        E=old.source(ref.q,ref.width,ref.offset,p)[1]
        edges.append(np.einsum('aij,pki,pklb,plj->pab',field.M2,U,E,U,optimize=True))
    return result+[np.array(edges)]

def main():
    if DEST.exists():raise FileExistsError(DEST)
    rng=np.random.default_rng(20261102);start=time.monotonic();cases=[];stop=None
    D=np.diag([1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2)]);iso=D@material.elastic@D
    A=rng.normal(size=(6,6));aniso=A.T@A+np.diag([1.,2.,3.,4.,5.,6.])
    try:
        for name,h,rot,chart in itertools.product(['RECT','SKEW','TRAP','SADDLE','NORMAL-C'],[.1,.01],[False,True],[False,True]):
            q=geometry.geometry(name)
            if rot:q=geometry.transform(q,geometry.ROTATION)
            width=h;offset=0.
            if chart:
                a=1.3;b=.1*h;K=np.array([[1,b,b*b/2],[0,a,a*b],[0,0,a*a]])
                q=np.einsum('ab,nbc->nac',K,q.reshape(9,3,3)).reshape(81);width=h/a;offset=-b/a
            ref=field.ReferenceReduction(q,width,offset,5)
            shapes=[(125,6),(36,6),(4,25,3)];x=tuple(rng.normal(size=s) for s in shapes);z=tuple(rng.normal(size=s) for s in shapes)
            snapshot=[a.copy() for a in x];cache={k:v.copy() for k,v in vars(ref).items() if isinstance(v,np.ndarray)}
            y=rng.normal(size=(125,6));rx=ref.apply(x);adj=ref.adjoint(y)
            lhs=np.sum(y*rx);rhs=dot(x,adj);dual=abs(lhs-rhs)/max(np.linalg.norm(y)*np.linalg.norm(rx),1e-30)
            separate=[]
            for k in range(3):
                single=tuple(v if i==k else np.zeros_like(v) for i,v in enumerate(x));r=ref.apply(single)
                separate.append(abs(np.sum(y*r)-dot(single,adj))/max(np.linalg.norm(y)*np.linalg.norm(r),1e-30))
            linear=relative(ref.apply(tuple(1.7*a-.4*b for a,b in zip(x,z))),1.7*rx-.4*ref.apply(z))
            mean=relative(ref.av@rx,ref.av@x[0])
            inputs_unchanged=all(np.array_equal(a,b) for a,b in zip(x,snapshot))
            rejects=0
            for bad in [(np.zeros((1,6)),x[1],x[2]),(x[0]*np.nan,x[1],x[2]),(x[0].astype(complex),x[1],x[2])]:
                try:ref.apply(bad)
                except ValueError:rejects+=1
            sources=source_linear(ref);new=np.stack([ref.apply(tuple(a[...,i] for a in sources)) for i in range(87)],axis=-1)
            expected=old.evaluate(old.prepare(q,width,offset,5),ref.points)[2];agreement=relative(new,expected)
            # Compatible affine natural state, followed by a finite spatial rotation.
            F0=np.array([[1.04,.03,-.02],[.01,.97,.02],[0.,-.01,1.02]])
            natural=(q.reshape(27,3)@F0.T).reshape(81);initial=field.green_sources(ref,natural)
            moved=geometry.transform(natural,geometry.ROTATION)
            natural_errors=[];derivatives=[];negative=[]
            trial=natural.copy();qc=trial.reshape(9,3,3)
            qc[:,0,2]+=.02*np.sin(np.arange(9));qc[:,1,0]+=.015*np.cos(np.arange(9));qc[:,2,1]+=.01*np.sin(2*np.arange(9))
            for C in [iso,aniso]:
                strain=ref.kernel(tuple(a-b for a,b in zip(field.green_sources(ref,moved),initial)))
                _,nf=field.energy_force(ref,moved,initial,C)
                natural_errors.append(max(float(np.linalg.norm(strain)/np.sqrt(strain.size)),float(np.linalg.norm(nf)/(np.linalg.norm(C)*ref.measure.sum()))))
                W,f=field.energy_force(ref,trial,initial,C);directions=rng.normal(size=(3,81));directions/=np.linalg.norm(directions,axis=1)[:,None]
                H=[];energy_errors=[];tangent_errors=[]
                for direction in directions:
                    hd,material_only=field.tangent_action(ref,trial,initial,C,direction)
                    complexW,complexf=field.energy_force(ref,trial.astype(complex)+1e-25j*direction,initial,C)
                    energy_errors.append(abs(float(complexW.imag/1e-25)-f@direction)/max(np.linalg.norm(f),1e-30))
                    tangent_errors.append(relative(hd,complexf.imag/1e-25));H.append(hd)
                    negative.append(relative(material_only,complexf.imag/1e-25))
                H=np.array(H).T;projected=directions@H
                derivatives.append(dict(energy=max(energy_errors),tangent=max(tangent_errors),symmetry=relative(projected,projected.T),energy_value=float(W)))
            cache_unchanged=all(np.array_equal(v,vars(ref)[k]) and not vars(ref)[k].flags.writeable for k,v in cache.items())
            checks=dict(linearity=linear<=1e-12,duality=max([dual]+separate)<=1e-12,mean=mean<=1e-12,old_agreement=agreement<=1e-10,input_immutable=inputs_unchanged,cache_immutable=cache_unchanged,rejection=rejects==3,natural=max(natural_errors)<=1e-10,energy=max(d['energy'] for d in derivatives)<=1e-9,tangent=max(d['tangent'] for d in derivatives)<=1e-9,symmetry=max(d['symmetry'] for d in derivatives)<=1e-10)
            record=dict(geometry=name,h=h,rotated=rot,chart=chart,linearity=linear,duality=dual,input_group_duality=separate,mean=mean,old_agreement=agreement,natural=max(natural_errors),derivatives=derivatives,omit_geometric_defect=max(negative),checks={k:bool(v) for k,v in checks.items()})
            cases.append(record)
            print(json.dumps(dict(geometry=name,h=h,rot=rot,chart=chart,pass_all=all(checks.values()),tangent=max(d['tangent'] for d in derivatives))),flush=True)
            if time.monotonic()-start>300 or memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
    except Exception as exc:stop=repr(exc)
    result=dict(date='2026-11-02',source_sha256=digest(__file__),dependencies={Path(m.__file__).name:digest(m.__file__) for m in [field,old,geometry,material]},protocol_sha256=digest(ROOT/'rd/structure/main-shell/q9-edge-trace/records/q9-edge-moment-field-experiment.md'),contract_sha256=digest(ROOT/'rd/structure/main-shell/q9-edge-trace/q9-edge-moment-field-interface.md'),cases=cases,stop_reason=stop,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),negative_control_pass=bool(cases and max(c['omit_geometric_defect'] for c in cases)>=1e-6),product_shell_adopted=False)
    result['all_criteria_pass']=bool(stop is None and len(cases)==40 and result['negative_control_pass'] and all(all(c['checks'].values()) for c in cases))
    DEST.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['cases','dependencies']}))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
