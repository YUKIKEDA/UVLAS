"""Independent Mandel reconstruction, least-squares edge fit and 11-point audit."""
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from scipy.linalg import lstsq
from numpy.polynomial.legendre import leggauss
import audit_q9_curved_strain_screen as independent
import check_q9_curved_strain_boundary as pair_geometry
import check_q9_independent_load as material
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification';MB=[]
for k in range(6):
    E=np.zeros((3,3))
    if k<3:E[k,k]=1
    else:
        i,j=[(1,2),(0,2),(0,1)][k-3];E[i,j]=E[j,i]=1/np.sqrt(2.)
    MB.append(E)
MB=np.array(MB);C=np.diag([1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2)])@material.elastic@np.diag([1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2)])
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def quad(n):
    p,w=leggauss(n);return np.array(list(itertools.product(p,p,p))),np.array([a*b*c for a,b,c in itertools.product(w,w,w)])
def fields(q,h,p):
    _,J,g=independent.data(q,h,p);B=np.einsum('kab,pib->pkia',MB,g).reshape(len(p),6,81)
    return J,np.concatenate((B,np.broadcast_to(np.eye(6),(len(p),6,6))),axis=2)
def transform(J):
    covariance=np.einsum('pia,kij,pjb->pkab',J,MB,J,optimize=True)
    return np.einsum('cab,pkab->pck',MB,covariance,optimize=True)
def defect(state,points):
    J,_=fields(state['q'],state['h'],points);H=transform(J);S,_,_=independent.interpolation(points)
    return H[:,[4,3],:]-np.einsum('pai,ik->pak',S,state['samples'][:,81:])
def edge(state,axis,sign):
    p,w=leggauss(5);points=np.array([[t,sign,z] if axis==0 else [sign,t,z] for t,z in itertools.product(p,p)]);J,E=fields(state['q'],state['h'],points);A=J[:,:,[axis,2]];G=np.einsum('pia,pib->pab',A,A);L=np.linalg.cholesky(G);U=np.linalg.solve(L,A.transpose(0,2,1)).transpose(0,2,1)
    H=np.stack((np.einsum('pi,kij,pj->pk',U[:,:,0],MB,U[:,:,0]),np.einsum('pi,kij,pj->pk',U[:,:,1],MB,U[:,:,1]),np.sqrt(2.)*np.einsum('pi,kij,pj->pk',U[:,:,0],MB,U[:,:,1])),axis=1)
    values=np.einsum('pik,pkn->pin',H,E);weight=np.array([a*b for a,b in itertools.product(w,w)])*np.sqrt(np.linalg.det(G));s=np.sqrt(weight/weight.sum());D=(s[:,None,None]*H).reshape(75,6);Y=(s[:,None,None]*values).reshape(75,87)
    coefficients,_,rank,_=lstsq(D,Y,cond=1e-10,lapack_driver='gelsd');assert rank>=3
    return coefficients
def first(state,points):
    points=np.asarray(points);J,E=fields(state['q'],state['h'],points);cov=np.einsum('pik,pkn->pin',transform(J),E);S,_,_=independent.interpolation(points);replacement=np.einsum('pai,in->pan',S,state['samples'])+np.einsum('pak,kn->pan',defect(state,points),state['mean'])
    for axis,sign in itertools.product([0,1],[-1,1]):
        boundary=points.copy();boundary[:,1-axis]=sign;K=defect(state,boundary)[:,axis,:];change=np.einsum('pk,kn->pn',K,state['edges'][axis,sign]-state['mean']);replacement[:,axis]+=((1+sign*points[:,1-axis])/2)[:,None]*change
    cov[:,4]=replacement[:,0];cov[:,3]=replacement[:,1];out=np.einsum('pik,pkn->pin',transform(np.linalg.inv(J)),cov)
    return J,E,out
def prepare(q,h):
    points,w=quad(9);S,ties,components=independent.interpolation(points);Js,Es=fields(q,h,ties);cov=np.einsum('pik,pkn->pin',transform(Js),Es);state=dict(q=q,h=h,samples=np.array([cov[i,4 if c==0 else 3] for i,c in enumerate(components)]),edges={})
    J,E=fields(q,h,points);measure=w*np.linalg.det(J);V=measure.sum();state['mean']=np.einsum('p,pkn->kn',measure,E)/V
    for key in itertools.product([0,1],[-1,1]):state['edges'][key]=edge(state,*key)
    _,_,out=first(state,points);state['moment']=np.einsum('p,pkn->kn',measure,out-E)/V;state['bubble']=float(measure@((1-points[:,0]**2)*(1-points[:,1]**2))/V)
    return state
def apply(state,points):
    points=np.asarray(points);J,E,R=first(state,points);b=(1-points[:,0]**2)*(1-points[:,1]**2)/state['bubble'];return J,E,R-b[:,None,None]*state['moment'][None]
def main():
    output=F/'q9-edge-moment-feasibility-audit.json'
    if output.exists():raise FileExistsError(output)
    start=time.monotonic();parent=F/'q9-edge-moment-feasibility-result.json';raw=json.loads(parent.read_text());assert raw['all_criteria_pass'];assert raw['source_sha256']==digest(ROOT/'tools/design-experiments/check_q9_edge_moment_feasibility.py');assert raw['candidate_sha256']==digest(ROOT/'docs/q9-edge-moment-candidate.md');assert raw['protocol_sha256']==digest(F/'q9-edge-moment-feasibility-experiment.md')
    for n,h in raw['dependency_sha256'].items():assert digest(ROOT/'tools/design-experiments'/n)==h
    for c in raw['cases']:assert digest(ROOT/c['artifact'])==c['artifact_sha256']
    records=[];pairs=[]
    for case in raw['cases']:
        if case['rotated'] or case['reparameterized']:continue
        d=np.load(ROOT/case['artifact']);q=d['q'];h=case['h'];state=prepare(q,h)
        coefficients=np.array([state['edges'][key] for key in itertools.product([0,1],[-1,1])]);coefferr=float(np.linalg.norm(coefficients-d['edge_coefficients'])/max(1.,np.linalg.norm(d['edge_coefficients'])));meanerr=float(np.linalg.norm(state['mean']-d['mean'])/max(1.,np.linalg.norm(d['mean'])));momenterr=float(np.linalg.norm(state['moment']-d['moment'])/max(1.,np.linalg.norm(d['moment'])));assert max(coefferr,meanerr,momenterr)<=1e-10
        points,w=quad(9);J,E,R=apply(state,points);measure=w*np.linalg.det(J);K=np.einsum('p,pik,ij,pjl->kl',measure,R[:,:,:81],C,R[:,:,:81],optimize=True);s=1/np.sqrt(d['rawK'].diagonal());kerr=float(np.linalg.norm(s[:,None]*(K-d['K'])*s[None,:])/np.linalg.norm(s[:,None]*d['K']*s[None,:]));assert kerr<=1e-10
        points,w=quad(11);J,E,R=apply(state,points);measure=w*np.linalg.det(J);work=np.einsum('p,pik->ki',measure,R[:,:,:81]);original=np.einsum('p,pik->ki',measure,E[:,:,:81]);workerr=float(np.max(np.linalg.norm(work-original,axis=0)/np.linalg.norm(original,axis=0)));assert workerr<=1e-10
        boundary=independent.boundary(q,h)@np.diag([1,1,1,1/np.sqrt(2),1/np.sqrt(2),1/np.sqrt(2)]);boundaryerr=float(np.linalg.norm(original-boundary)/np.linalg.norm(boundary));assert boundaryerr<=1e-10
        value=float(np.max(np.linalg.norm(R[:,:,81:]-np.eye(6)[None],axis=1)));assert value<=1e-10
        records.append(dict(geometry=case['geometry'],h=h,coefficient_difference=coefferr,mean_difference=meanerr,moment_difference=momenterr,stiffness_difference=kerr,work11_difference=workerr,boundary_traction_difference=boundaryerr,value_difference=value))
        assert time.monotonic()-start<=300 and memory.memory_bytes()<=1024**3
    for name,h in itertools.product(['RECT','SKEW','TRAP','SADDLE','NORMAL-C'],[.1,.01]):
        coords,ids=pair_geometry.two_elements(name);left=np.array([[1,t,z] for t,z in itertools.product([-.8,-.1,.6],[-.6,0,.7])]);right=left.copy();right[:,0]=-1;ops=[]
        for side,p in enumerate([left,right]):
            state=prepare(coords[ids[side]].reshape(81),h);J,_,R=apply(state,p);S=np.einsum('pi,kij,pj->pk',J[:,:,1],MB,J[:,:,2]);value=np.einsum('pi,pij->pj',S,R[:,:,:81]);op=np.zeros((9,135));op[:,(9*ids[side][:,None]+np.arange(9)).reshape(-1)]=value;ops.append(op)
        error=float(np.linalg.norm(ops[0]-ops[1])/max(np.linalg.norm(ops[0]),np.linalg.norm(ops[1])));old=next(c for c in raw['shared_pairs'] if c['geometry']==name and c['h']==h);assert error<=1e-10 and abs(error-old['trace_error'])<=1e-10;pairs.append(dict(geometry=name,h=h,trace_error=error))
        assert time.monotonic()-start<=300 and memory.memory_bytes()<=1024**3
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [independent,pair_geometry,material,memory]},cases=records,shared_pairs=pairs,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=bool(len(records)==10 and len(pairs)==10),product_shell_adopted=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_criteria_pass=result['all_criteria_pass'],cases=len(records),pairs=len(pairs),maximum_work11_error=max(c['work11_difference'] for c in records),maximum_stiffness_difference=max(c['stiffness_difference'] for c in records),seconds=result['seconds'])))
if __name__=='__main__':main()
