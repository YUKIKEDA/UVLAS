"""Independent adjoint stress work and boundary traction audit of curved screen."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from numpy.polynomial.legendre import leggauss
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;OFF=[(-1,-1),(1,-1),(1,1),(-1,1),(0,-1),(1,0),(0,1),(-1,0),(0,0)];A=1/np.sqrt(3.)
E=[];ST=[]
for k in range(6):
    e=np.zeros((3,3));s=np.zeros((3,3))
    if k<3:e[k,k]=s[k,k]=1
    else:
        i,j=[(1,2),(0,2),(0,1)][k-3];e[i,j]=e[j,i]=.5;s[i,j]=s[j,i]=1
    E.append(e);ST.append(s)
E=np.array(E);ST=np.array(ST)
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def lag(t,i):return (.5*t*(t-1),t-.5) if i==-1 else ((1-t*t,-2*t) if i==0 else (.5*t*(t+1),t+.5))
def data(q,h,points):
    phis=[];gs=[]
    for s,t,v in points:
        ls=[lag(s,i) for i,_ in OFF];lt=[lag(t,j) for _,j in OFF];n=np.array([x[0]*y[0] for x,y in zip(ls,lt)]);ns=np.array([x[1]*y[0] for x,y in zip(ls,lt)]);nt=np.array([x[0]*y[1] for x,y in zip(ls,lt)]);z=h*v/2
        phis.append(np.array([[r,r*z,r*z*z/2] for r in n]).reshape(27));gs.append(np.stack((np.array([[r,r*z,r*z*z/2] for r in ns]).reshape(27),np.array([[r,r*z,r*z*z/2] for r in nt]).reshape(27),np.array([[0.,r,r*z] for r in n]).reshape(27)*h/2),axis=1))
    phis=np.array(phis);g=np.array(gs);J=np.einsum('ia,pib->pab',q.reshape(27,3),g);grad=np.einsum('pia,pab->pib',g,np.linalg.inv(J))
    return phis,J,grad
def L(x,nodes,i):
    result=1.
    for j,v in enumerate(nodes):
        if j!=i:result*= (x-v)/(nodes[i]-v)
    return result
def interpolation(points):
    S=np.zeros((len(points),2,36));ties=[];comps=[];column=0
    for c,xs,ys in [(0,[-A,A],[-1.,0.,1.]),(1,[-1.,0.,1.],[-A,A])]:
        for i,j,k in itertools.product(range(len(xs)),range(len(ys)),range(3)):
            for p,(x,y,z) in enumerate(points):S[p,c,column]=L(x,xs,i)*L(y,ys,j)*L(z,[-1.,0.,1.],k)
            ties.append([xs[i],ys[j],[-1.,0.,1.][k]]);comps.append(c);column+=1
    return S,np.array(ties),comps
def strain_matrix(A):
    cov=np.einsum('pia,kij,pjb->pkab',A,E,A)
    return np.stack((cov[:,:,0,0],cov[:,:,1,1],cov[:,:,2,2],2*cov[:,:,1,2],2*cov[:,:,0,2],2*cov[:,:,0,1]),axis=1)
def boundary(q,h):
    p,w=leggauss(9);points=[];weights=[];axes=[];signs=[]
    for axis,sign in itertools.product(range(3),[-1.,1.]):
        for i,j in itertools.product(range(9),repeat=2):
            point=np.zeros(3);other=[k for k in range(3) if k!=axis];point[axis]=sign;point[other]=[p[i],p[j]];points.append(point);weights.append(w[i]*w[j]);axes.append(axis);signs.append(sign)
    N,J,_=data(q,h,points);area=np.array([s*np.cross(j[:,(a+1)%3],j[:,(a+2)%3]) for j,a,s in zip(J,axes,signs)]);traction=np.einsum('kab,pb->pak',ST,area)
    return np.einsum('p,pi,pak->iak',weights,N,traction).reshape(81,6)
def main():
    target=F/'q9-curved-strain-screen-audit.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();parent=F/'q9-curved-strain-screen-result.json';raw=json.loads(parent.read_text());assert raw['diagnostic_checks_pass'];assert raw['source_sha256']==digest(ROOT/'rd/structure/main-shell/q9-curved/experiments/check_q9_curved_strain_screen.py');assert raw['protocol_sha256']==digest(F/'q9-curved-strain-screen-experiment.md')
    for n,h in raw['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/n)==h
    for c in raw['cases']:assert digest(ROOT/c['artifact'])==c['artifact_sha256']
    p,w=leggauss(9);points=np.array(list(itertools.product(p,p,p)));weights=np.array([a*b*c for a,b,c in itertools.product(w,w,w)]);S,ties,components=interpolation(points);records=[]
    for case in raw['cases']:
        if case['rotated'] or case['reparameterized']:continue
        saved=np.load(ROOT/case['artifact']);q=saved['q'];h=case['h'];_,J,grad=data(q,h,points);_,Js,gs=data(q,h,ties);_,Jc,_=data(q,h,[[0.,0.,0.]]);Ag=J if case['scheme']=='G-COV' else np.broadcast_to(Jc,J.shape);At=Js if case['scheme']=='G-COV' else np.broadcast_to(Jc,Js.shape);inv=np.linalg.inv(Ag);weight=weights*np.linalg.det(J)
        # Pull physical stress into the dual of the covariant strain, keeping only retained components.
        dual=np.einsum('pai,kij,pbj->pabk',inv,ST,inv,optimize=True);retained=dual.copy();retained[:,0,2]=retained[:,2,0]=0;retained[:,1,2]=retained[:,2,1]=0
        physical=np.einsum('pia,pabk,pjb->pijk',Ag,retained,Ag,optimize=True);force=np.einsum('p,pib,pabk->iak',weight,grad,physical,optimize=True)
        lambdas=np.array([np.einsum('p,pk,p->k',weight,2*dual[:,c,2],S[:,c,i]) for i,c in enumerate(components)])
        moments=np.array([(np.outer(a[:,c],a[:,2])+np.outer(a[:,2],a[:,c]))/2 for a,c in zip(At,components)])
        force+=np.einsum('pib,pab,pk->iak',gs,moments,lambdas,optimize=True);force=force.reshape(81,6)
        direct=np.einsum('p,pib,kab->iak',weight,grad,ST,optimize=True).reshape(81,6);edge=boundary(q,h)
        defect=float(np.linalg.norm(force-saved['work'])/np.linalg.norm(saved['work']));boundary_difference=float(np.linalg.norm(direct-edge)/np.linalg.norm(edge));raw_difference=float(np.linalg.norm(direct-saved['raw_work'])/np.linalg.norm(direct));assert max(defect,boundary_difference,raw_difference)<=1e-10
        # Independent 6x6 engineering transformation matrices, without reducing displacement B.
        V=strain_matrix(Ag);Vs=strain_matrix(At);beta=np.array([Vs[i,4 if c==0 else 3] for i,c in enumerate(components)]);replacement=np.einsum('pai,ik->pak',S,beta,optimize=True);V[:,4]=replacement[:,0];V[:,3]=replacement[:,1];O=np.einsum('pij,pjk->pik',strain_matrix(inv),V,optimize=True)
        value_difference=float(np.max(abs(O-saved['O'])));assert value_difference<=1e-10
        records.append(dict(geometry=case['geometry'],h=h,scheme=case['scheme'],adjoint_force_difference=defect,boundary_traction_difference=boundary_difference,raw_force_difference=raw_difference,value_map_difference=value_difference))
        assert time.monotonic()-start<=300 and memory.memory_bytes()<=1024**3
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),memory_helper_sha256=digest(memory.__file__),cases=records,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=len(records)==20,operators_accepted=False,product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_criteria_pass=result['all_criteria_pass'],cases=len(records),maximum_adjoint_difference=max(c['adjoint_force_difference'] for c in records),maximum_boundary_difference=max(c['boundary_traction_difference'] for c in records),seconds=result['seconds'])))
if __name__=='__main__':main()
