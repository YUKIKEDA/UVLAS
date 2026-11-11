"""Prospective Q3 representation, rank, affine work and shared trace screen."""
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
from scipy.linalg import svdvals
import q3_curved_element as element
import q3_curved_independent as independent
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def energy(v):return np.einsum('...i,ij,...j->...',v,element.C,v)
def shapes():
    p=np.array(list(itertools.product(np.linspace(-1,1,9),repeat=2)));N,g=element.surface(p);Ni,gi=independent.shape(p);nodes=np.array([[s,t] for t,s in itertools.product(element.NODES,element.NODES)]);at,_=element.surface(nodes);error=max(np.max(abs(N-Ni)),np.max(abs(g-gi)),np.max(abs(at-np.eye(16))),np.max(abs(N.sum(axis=1)-1)),np.max(abs(g.sum(axis=1))))
    for i,j in itertools.product(range(4),repeat=2):
        values=nodes[:,0]**i*nodes[:,1]**j;expected=p[:,0]**i*p[:,1]**j;ds=i*p[:,0]**(i-1)*p[:,1]**j if i else np.zeros(len(p));dt=j*p[:,0]**i*p[:,1]**(j-1) if j else np.zeros(len(p));error=max(error,np.max(abs(N@values-expected)),np.max(abs(g[:,:,0]@values-ds)),np.max(abs(g[:,:,1]@values-dt)))
    return float(error)
def expected_geometry(k,h,points,corners):
    xy,g=element.bilinear(points,corners);x,y=xy.T;z=h*points[:,2]/2;position=np.column_stack((x-k*x*z,y,k*x*x/2+z));J=np.zeros((len(points),3,3));J[:,0,0]=1-k*z;J[:,1,1]=1;J[:,2,0]=k*x;J[:,0,2]=-k*x;J[:,2,2]=1;chart=np.zeros_like(J);chart[:,:2,:2]=g;chart[:,2,2]=h/2
    return xy,position,J@chart
def local(q,u,h,k,corners):
    observations={}
    for order in [9,11]:
        points,w=element.quadrature(order);phi,J,B=element.evaluate(q,h,points);weight=w*np.linalg.det(J);xy,pos,expectedJ=expected_geometry(k,h,points,corners);x=xy[:,0];z=h*points[:,2]/2;star=element.target(k,x,z);error=B@u-star;ref=float(weight@energy(star));norm=float(np.sqrt(weight@energy(error)/ref));observations[order]=(phi,J,B,weight,ref,norm)
    phi,J,B,weight,ref,norm=observations[11];points,_=element.quadrature(11);xy,pos,expectedJ=expected_geometry(k,h,points,corners);x=xy[:,0];z=h*points[:,2]/2;position=float(np.max(abs(phi@q.reshape(48,3)-pos)));jacobian=float(np.max(np.linalg.norm(J-expectedJ,axis=1)/np.maximum(np.linalg.norm(expectedJ,axis=1),1e-300)));value=np.column_stack((-k*x**3/3-x*z,np.zeros(len(x)),x*x/2-k*x*x*z));displacement=float(np.max(abs(phi@u.reshape(48,3)-value)))
    affine=np.einsum('kij,nj->nik',element.MB,q.reshape(48,3)).reshape(144,6);constant=float(np.max(abs(B@affine-np.eye(6))));work=np.einsum('p,pki->ik',weight,B);boundary=independent.boundary(q,h);work_difference=float(np.linalg.norm(work-boundary)/np.linalg.norm(boundary));_,Ji,grad=independent.data(q,h,points);Bi=np.einsum('kij,pnj->pkni',element.MB,grad).reshape(-1,6,144);independent_difference=float(np.linalg.norm(B-Bi)/np.linalg.norm(B))
    phi9,J9,B9,w9,ref9,norm9=observations[9];A=(np.einsum('ij,pjk->pik',element.CH,B9)*np.sqrt(w9)[:,None,None]).reshape(-1,144);scale=1/np.linalg.norm(A,axis=0);As=A*scale;s=svdvals(As);ranks=[int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]];rigids=[]
    for axis in np.eye(3):
        v=np.zeros((16,3,3));v[:,0]=axis;rigids.extend([v.reshape(144),np.cross(np.broadcast_to(axis,(16,3,3)),q.reshape(16,3,3)).reshape(144)])
    rigid=max(float(np.linalg.norm(As@(v/scale))/(s[0]*np.linalg.norm(v/scale))) for v in rigids);quadrature=abs(ref9/ref-1.)
    checks=dict(geometry=max(position,jacobian)<=1e-11,displacement=displacement<=1e-11,constant=constant<=1e-9,work=work_difference<=1e-9,rank=ranks==[138]*3,rigid=rigid<=1e-10,representation=max(norm,norm9)<=1e-8,quadrature=quadrature<=1e-10,independent_gradient=independent_difference<=1e-10)
    return dict(position=position,jacobian=jacobian,displacement=displacement,constant=constant,work_difference=work_difference,ranks=ranks,smallest_retained_ratio=float(s[137]/s[0]),largest_rigid_ratio=float(s[138]/s[0]),rigid=rigid,energy_error=norm,energy_error9=norm9,reference_norm=ref,quadrature=quadrature,independent_difference=independent_difference,checks=checks),work
def traces(cells,coords,n):
    rng=np.random.default_rng(20261102);u=rng.standard_normal((len(coords),9));lookup={(i,j):ids for i,j,ids in cells};p=np.array(list(itertools.product(np.linspace(-1,1,21),[-1.,-.5,0.,.5,1.])));maximum=0.;count=0
    for axis in [0,1]:
        for j,i in itertools.product(range(n),range(n)):
            if (i if axis==0 else j)==n-1:continue
            neighbor=(i+1,j) if axis==0 else (i,j+1);ids0=lookup[i,j];ids1=lookup[neighbor];outputs=[]
            for sign,ids in [(1.,ids0),(-1.,ids1)]:
                points=np.zeros((len(p),3));points[:,axis]=sign;points[:,1-axis]=p[:,0];points[:,2]=p[:,1];phi,g=element.tables(points,1.);J=np.einsum('na,pnb->pab',coords[ids].reshape(48,3),g);dJ=np.einsum('na,pnb->pab',u[ids].reshape(48,3),g);cov=(J.transpose(0,2,1)@dJ+dJ.transpose(0,2,1)@J)/2;idx=[1-axis,2];outputs.append((phi@u[ids].reshape(48,3),dJ[:,:,idx],cov[:,idx,:][:,:,idx]))
            maximum=max(maximum,max(float(np.max(abs(a-b))) for a,b in zip(outputs[0],outputs[1])));count+=1
    return maximum,count
def main():
    dest=ROOT/'docs/verification/q3-curved-feasibility-result.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();out=ROOT/'build/analysis'/('q3-curved-feasibility-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();basis=shapes();records=[];stop=None
    try:
        if basis>1e-11:raise ValueError('basis_gate')
        for k,amplitude,n,h in itertools.product([0.,.4],[0.,.2],[1,2],[.01,.0001]):
            xy,cells,coords,u=element.mesh(k,n,amplitude);volume=np.zeros((len(xy),9,6));exterior=np.zeros_like(volume);logs=[]
            for ci,(i,j,ids) in enumerate(cells):
                q=coords[ids].reshape(144);log,work=local(q,u[ids].reshape(144),h,k,xy[ids][[0,3,15,12]]);logs.append(dict(cell=ci,**log));volume[ids]+=work.reshape(16,9,6);faces=[(2,-1),(2,1)]
                for axis,index in [(0,i),(1,j)]:
                    if index==0:faces.append((axis,-1))
                    if index==n-1:faces.append((axis,1))
                exterior[ids]+=independent.boundary(q,h,faces).reshape(16,9,6)
            work=float(np.linalg.norm(volume-exterior)/np.linalg.norm(exterior));trace,count=traces(cells,coords,n);support=float(np.max(abs(u[np.arange(len(u))%(3*n+1)==0])));artifact=out/f'k{k:g}-a{amplitude:g}-n{n}-h{h:g}.npz';np.savez_compressed(artifact,xy=xy,coords=coords,u=u,ids=np.array([c[2] for c in cells]),ij=np.array([c[:2] for c in cells]));checks=dict(local=all(all(c['checks'].values()) for c in logs),shared_work=work<=1e-9,shared_trace=trace<=1e-10,support=support<=1e-14,resources=time.monotonic()-start<=300 and memory.memory_bytes()<=2*1024**3)
            records.append(dict(kappa=k,amplitude=amplitude,n=n,h=h,local=logs,shared_work_difference=work,shared_trace_difference=trace,internal_edges=count,support_error=support,checks=checks,artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=digest(artifact)));print(k,amplitude,n,h,checks,'rank',logs[-1]['ranks'],'error',max(c['energy_error'] for c in logs),flush=True)
            if not all(checks.values()):stop='required_gate_failed';break
    except Exception as exc:stop=repr(exc)
    sources=[Path(__file__),Path(element.__file__),Path(independent.__file__),Path(element.material.__file__),Path(element.old_grid.__file__),ROOT/'docs/q3-curved-bending-comparison-contract.md',ROOT/'docs/verification/q3-curved-feasibility-experiment.md']
    result=dict(candidate='compatible-Q3-P2-144',physical_adoption=False,basis_error=basis,records=records,all_checks_pass=len(records)==16 and stop is None,stop_reason=stop,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes());dest.write_text(json.dumps(result,indent=2),encoding='utf-8')
if __name__=='__main__':main()
