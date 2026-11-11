"""Independent engineering maps and pointwise adjoint for shared-edge screen."""
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from numpy.polynomial.legendre import leggauss
import audit_q9_curved_strain_screen as independent
import check_q9_curved_load as geometry
import check_q9_curved_strain_boundary as definition
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def maps(q,h,points,scheme):
    _,J,g=independent.data(q,h,points);S,ties,components=independent.interpolation(points);_,Js,gs=independent.data(q,h,ties);_,Jc,_=independent.data(q,h,[[0,0,0]]);A=J if scheme=='G-COV' else np.broadcast_to(Jc,J.shape);At=Js if scheme=='G-COV' else np.broadcast_to(Jc,Js.shape);inv=np.linalg.inv(A)
    V=independent.strain_matrix(A);Vs=independent.strain_matrix(At);samples=np.array([Vs[i,4 if c==0 else 3] for i,c in enumerate(components)]);replacement=np.einsum('pai,ik->pak',S,samples);V[:,4]=replacement[:,0];V[:,3]=replacement[:,1];O=np.einsum('pij,pjk->pik',independent.strain_matrix(inv),V)
    sigma=(np.einsum('pi,pj->pij',J[:,:,1],J[:,:,2])+np.einsum('pi,pj->pij',J[:,:,2],J[:,:,1]))/2
    dual=np.einsum('pai,pij,pbj->pab',inv,sigma,inv);retained=dual.copy();retained[:,0,2]=retained[:,2,0]=0;retained[:,1,2]=retained[:,2,1]=0;physical=np.einsum('pia,pab,pjb->pij',A,retained,A)
    operator=np.einsum('pib,pab->pia',g,physical).reshape(len(points),81);moments=np.array([(np.outer(a[:,c],a[:,2])+np.outer(a[:,2],a[:,c]))/2 for a,c in zip(At,components)])
    sample_operator=np.einsum('pib,pab->pia',gs,moments).reshape(36,81);coefficients=np.column_stack([2*dual[:,c,2]*S[:,c,i] for i,c in enumerate(components)]);operator+=coefficients@sample_operator
    return O,operator
def error(O):
    w=np.array([1,1,1,.5,.5,.5]);d=O-np.eye(6)[None];return float(np.max(np.sqrt(np.sum(d*d*w[None,:,None],axis=1)/w[None,:])))
def main():
    out=F/'q9-curved-strain-boundary-audit.json'
    if out.exists():raise FileExistsError(out)
    start=time.monotonic();path=F/'q9-curved-strain-boundary-result.json';data=json.loads(path.read_text());assert data['diagnostic_checks_pass'] and data['source_sha256']==digest(definition.__file__);assert data['protocol_sha256']==digest(F/'q9-curved-strain-boundary-experiment.md')
    for name,h in data['dependency_sha256'].items():assert digest(ROOT/'tools/design-experiments'/name)==h
    p,_=leggauss(9);points=[]
    for axis,sign in itertools.product(range(2),[-1,1]):
        for t,z in itertools.product(p,[-1,0,1]):points.append([sign,t,z] if axis==0 else [t,sign,z])
    records=[]
    for case in data['cases']:
        name=case['geometry'];h=case['h'];scheme=case['scheme'];q=geometry.geometry(name);O,_=maps(q,h,points,scheme);value=error(O)
        coords,ids=definition.two_elements(name);left=np.array([[1,t,z] for t,z in itertools.product([-.8,-.1,.6],[-.6,0,.7])]);right=left.copy();right[:,0]=-1
        _,a=maps(coords[ids[0]].reshape(81),h,left,scheme);_,b=maps(coords[ids[1]].reshape(81),h,right,scheme);A=np.zeros((9,135));B=A.copy();A[:,(9*ids[0][:,None]+np.arange(9)).reshape(-1)]=a;B[:,(9*ids[1][:,None]+np.arange(9)).reshape(-1)]=b
        trace=float(np.linalg.norm(A-B)/max(np.linalg.norm(A),np.linalg.norm(B)));dv=abs(value-case['boundary_constant_value_error']);dt=abs(trace-case['shared_trace_error']);assert max(dv,dt)<=1e-10
        assert (value<=1e-10)==case['value_gate'] and (trace<=1e-10)==case['trace_gate'];records.append(dict(geometry=name,h=h,scheme=scheme,value_difference=dv,trace_difference=dt))
        assert time.monotonic()-start<=120 and memory.memory_bytes()<=512*1024**2
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(path),dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [independent,geometry,definition,memory]},cases=records,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=len(records)==20,operators_accepted=False)
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_criteria_pass=result['all_criteria_pass'],maximum_difference=max(max(c['value_difference'],c['trace_difference']) for c in records),seconds=result['seconds'])))
if __name__=='__main__':main()
