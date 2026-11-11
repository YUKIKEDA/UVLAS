"""Separate 2D derivatives and exact thickness mass for trace candidate gates."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from numpy.polynomial.legendre import leggauss
import check_q9_curved_load as original
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification';OFF=[(-1,-1),(1,-1),(1,1),(-1,1),(0,-1),(1,0),(0,1),(-1,0),(0,0)];A=1/np.sqrt(3.)
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def lag(t,i):return (.5*t*(t-1),t-.5) if i==-1 else ((1-t*t,-2*t) if i==0 else (.5*t*(t+1),t+.5))
def shapes(s,t):
    a=[lag(s,i) for i,_ in OFF];b=[lag(t,j) for _,j in OFF]
    return np.array([x[0]*y[0] for x,y in zip(a,b)]),np.array([x[1]*y[0] for x,y in zip(a,b)]),np.array([x[0]*y[1] for x,y in zip(a,b)])
def samples(xy,h):
    rows=[]
    for comp,ss,tt in [(0,[-A,A],[-1.,0.,1.]),(1,[-1.,0.,1.],[-A,A])]:
        for s in ss:
            for t in tt:
                N,ds,dt=shapes(s,t);der=ds if comp==0 else dt;tangent=der@xy
                for zz in [-1.,0.,1.]:
                    z=h*zz/2;r=np.zeros((9,9));r[:,2]=der;r[:,3]=tangent[0]*N;r[:,4]=tangent[1]*N;r[:,5]=z*der;r[:,6]=z*tangent[0]*N;r[:,7]=z*tangent[1]*N;r[:,8]=z*z/2*der;rows.append(r.reshape(81))
    return np.array(rows)
def metric(xy):
    points,weights=leggauss(9);M=np.zeros((12,12))
    def L(t,nodes,k):
        v=1.
        for j,x in enumerate(nodes):
            if j!=k:v*=(t-x)/(nodes[k]-x)
        return v
    for s,ws in zip(points,weights):
        for t,wt in zip(points,weights):
            _,ds,dt=shapes(s,t);J=xy.T@np.column_stack((ds,dt));S=np.zeros((2,12));col=0
            for comp,ss,tt in [(0,[-A,A],[-1.,0.,1.]),(1,[-1.,0.,1.],[-A,A])]:
                for i in range(len(ss)):
                    for j in range(len(tt)):S[comp,col]=L(s,ss,i)*L(t,tt,j);col+=1
            V=np.linalg.solve(J.T,S);M+=V.T@V*ws*wt*np.linalg.det(J)/2.6
    return M
def main():
    output=F/'q9-trace-shear-gates-audit.json'
    if output.exists():raise FileExistsError(output)
    start=time.monotonic();localpath=F/'q9-trace-shear-local-result.json';sharedpath=F/'q9-trace-shear-shared-result.json';local=json.loads(localpath.read_text());shared=json.loads(sharedpath.read_text());records=[];local_records=[]
    for data,script in [(local,'check_q9_trace_shear_local.py'),(shared,'check_q9_trace_shear_shared.py')]:
        assert data['all_criteria_pass'] and data['source_sha256']==digest(ROOT/'tools/design-experiments'/script)
        assert data['candidate_sha256']==digest(ROOT/'docs/q9-trace-shear-candidate.md') and data['protocol_sha256']==digest(F/'q9-trace-shear-gate-experiment.md')
        for n,h in data['dependency_sha256'].items():assert digest(ROOT/'tools/design-experiments'/n)==h
    for c in local['cases']:
        assert digest(ROOT/c['artifact'])==c['artifact_sha256'];d=np.load(ROOT/c['artifact']);xy=d['q'].reshape(9,3,3)[:,0,:2];M0=metric(xy);T=samples(xy,.1);M=np.kron(M0,np.array([[4.,2.,-1.],[2.,16.,2.],[-1.,2.,4.]])*.1/30)
        derr=float(np.linalg.norm(M-d['M'])/np.linalg.norm(M));terr=float(np.linalg.norm(T-d['T'])/np.linalg.norm(T));assert max(derr,terr)<=1e-10
        full=original.evaluate(d['q'],.1,9);_,_,w=original.volume_table(9,.1);Bs=full['b'][:,[4,3],:];Ks=np.einsum('p,pai,paj->ij',w*np.linalg.det(full['jac'])/2.6,Bs,Bs);K=full['stiffness']-Ks+T.T@M@T
        scale=1/np.sqrt(d['Kfull'].diagonal());err=float(np.linalg.norm(scale[:,None]*(K-d['Ktie'])*scale[None,:])/np.linalg.norm(scale[:,None]*K*scale[None,:]));assert err<=1e-10
        local_records.append(dict(case=c['case'],metric_difference=derr,sampling_difference=terr,stiffness_difference=err))
    for c in shared['cases']:
        assert digest(ROOT/c['artifact'])==c['artifact_sha256'];d=np.load(ROOT/c['artifact']);nodes=d['nodes'];n=c['n'];side=2*n+1;wm={int(v):i for i,v in enumerate(d['w_nodes'])};dm={tuple(map(int,v)):i for i,v in enumerate(d['d_node_components'])};D=np.zeros_like(d['D']);W=np.zeros_like(d['W']);tie_error=0.
        for j in range(n):
            for i in range(n):
                index=j*n+i;ids=np.array([(2*j+y+1)*side+2*i+x+1 for x,y in OFF]);xy=nodes[ids];T0=samples(xy,.1)[1::3];M0=metric(xy);C=np.linalg.cholesky(M0).T@T0;rows=slice(12*index,12*(index+1));tie_error=max(tie_error,float(np.linalg.norm(T0-d['ties'][index])/max(1.,np.linalg.norm(T0))))
                for k,node in enumerate(ids):
                    node=int(node)
                    if node in wm:W[rows,wm[node]]+=C[:,9*k+2]/n
                    for comp in [3,4]:
                        if (node,comp) in dm:D[rows,dm[node,comp]]+=C[:,9*k+comp]
        ed=float(np.linalg.norm(D-d['D'])/np.linalg.norm(D));ew=float(np.linalg.norm(W-d['W'])/np.linalg.norm(W));assert max(ed,ew,tie_error)<=1e-10
        U,s,_=np.linalg.svd(D,full_matrices=True);combined=np.linalg.svd(np.column_stack((D,W)),compute_uv=False);norm=max(s[0],combined[0]);dims={}
        for t in [1e-9,1e-10,1e-11]:
            rank=int(np.count_nonzero(s>t*norm));z=np.linalg.svd(U[:,rank:].T@W,compute_uv=False);dim=W.shape[1]-int(np.count_nonzero(z>t*norm));dims[str(t)]=dim;assert dim==len(wm)==c['ranks'][str(t)]['compatible_w_dimension']
        records.append(dict(n=n,amplitude=c['amplitude'],D_difference=ed,W_difference=ew,tie_difference=tie_error,compatible_w_dimensions=dims))
        assert time.monotonic()-start<=120 and memory.memory_bytes()<=512*1024**2
    result=dict(source_sha256=digest(__file__),local_result_sha256=digest(localpath),shared_result_sha256=digest(sharedpath),dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [original,memory]},local_cases=local_records,shared_cases=records,
        seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=True,product_shell_adopted=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_criteria_pass=True,local_cases=local_records,maximum_shared_difference=max(max(c['D_difference'],c['W_difference'],c['tie_difference']) for c in records),seconds=result['seconds'])))
if __name__=='__main__':main()
