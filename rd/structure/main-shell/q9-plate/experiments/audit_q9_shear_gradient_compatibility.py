"""Independent 2D moment assembly, versus leading rows of 3D P-S36."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve_triangular
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
OFF=[(-1,-1),(1,-1),(1,1),(-1,1),(0,-1),(1,0),(0,1),(-1,0),(0,0)]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def lag(t,i):return (.5*t*(t-1),t-.5) if i==-1 else ((1-t*t,-2*t) if i==0 else (.5*t*(t+1),t+.5))
def main():
    start=time.monotonic();path=F/'q9-shear-gradient-compatibility-result.json';data=json.loads(path.read_text());records=[]
    assert data['all_criteria_pass'] and data['source_sha256']==digest(ROOT/'rd/structure/main-shell/q9-plate/experiments/check_q9_shear_gradient_compatibility.py')
    assert data['protocol_sha256']==digest(ROOT/'rd/structure/main-shell/q9-plate/q9-shear-gradient-compatibility-contract.md')
    assert data['candidate_sha256']==digest(ROOT/'rd/structure/main-shell/q9-plate/q9-shear-projection-candidate.md')
    for n,h in data['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/n)==h
    p,w=leggauss(9);N=[];dx=[];dy=[];S=[];weight=[]
    for i in range(9):
        for j in range(9):
            a=[lag(p[i],x) for x,_ in OFF];b=[lag(p[j],y) for _,y in OFF]
            N.append([x[0]*y[0] for x,y in zip(a,b)]);dx.append([x[1]*y[0] for x,y in zip(a,b)]);dy.append([x[0]*y[1] for x,y in zip(a,b)])
            si=np.zeros((2,12));col=0
            P=lambda t,k:1. if k==0 else (t if k==1 else (3*t*t-1)/2)
            for component,ni,nj in [(0,2,3),(1,3,2)]:
                for ii in range(ni):
                    for jj in range(nj):si[component,col]=P(p[i],ii)*P(p[j],jj);col+=1
            S.append(si);weight.append(w[i]*w[j])
    N=np.array(N);dx=np.array(dx);dy=np.array(dy);S=np.array(S);weight=np.array(weight);mu=1/2.6
    for c in data['cases']:
        assert digest(ROOT/c['artifact'])==c['artifact_sha256'];arrays=np.load(ROOT/c['artifact']);n=c['n'];side=2*n+1;nodes=arrays['nodes'];wmap={int(node):i for i,node in enumerate(arrays['w_nodes'])};dmap={tuple(map(int,key)):i for i,key in enumerate(arrays['d_node_components'])}
        D=np.zeros_like(arrays['D']);W=np.zeros_like(arrays['W'])
        for j in range(n):
            for i in range(n):
                ids=np.array([(2*j+y+1)*side+2*i+x+1 for x,y in OFF]);xy=nodes[ids];jac=np.stack((dx@xy,dy@xy),axis=2);measure=mu*weight*np.linalg.det(jac)
                V=np.linalg.solve(jac.transpose(0,2,1),S);grad=np.linalg.solve(jac.transpose(0,2,1),np.stack((dx,dy),axis=1)).transpose(0,2,1)
                Bs=np.zeros((81,2,27));Bs[:,:,0::3]=grad.transpose(0,2,1);Bs[:,0,1::3]=N;Bs[:,1,2::3]=N
                M=np.einsum('p,pai,paj->ij',measure,V,V);H=np.einsum('p,pai,paj->ij',measure,V,Bs);C=solve_triangular(np.linalg.cholesky(M),H,lower=True);rows=slice(12*(j*n+i),12*(j*n+i+1))
                for k,node in enumerate(ids):
                    node=int(node)
                    if node in wmap:W[rows,wmap[node]]+=C[:,3*k]/n
                    for component in [3,4]:
                        if (node,component) in dmap:D[rows,dmap[node,component]]+=C[:,3*k+component-2]
        diffD=float(np.linalg.norm(D-arrays['D'])/np.linalg.norm(D));diffW=float(np.linalg.norm(W-arrays['W'])/np.linalg.norm(W));assert max(diffD,diffW)<=1e-10
        U,s,_=np.linalg.svd(D,full_matrices=True);combined=np.column_stack((D,W));sc=np.linalg.svd(combined,compute_uv=False);scale=max(s[0],sc[0]);dims={}
        for t in [1e-9,1e-10,1e-11]:
            rank=int(np.count_nonzero(s>t*scale));Z=U[:,rank:].T@W;sz=np.linalg.svd(Z,compute_uv=False);dim=W.shape[1]-int(np.count_nonzero(sz>t*scale));dims[str(t)]=dim
            assert dim==c['ranks'][str(t)]['compatible_w_dimension']
        smallest=float(sz[-1]);records.append(dict(n=n,amplitude=c['amplitude'],relative_D_difference=diffD,relative_W_difference=diffW,compatible_w_dimensions=dims,
            smallest_defect_singular_value=smallest,operator_scale=float(scale),smallest_relative_to_operator=smallest/scale))
        assert time.monotonic()-start<=120 and memory.memory_bytes()<=512*1024**2
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(path),memory_helper_sha256=digest(memory.__file__),cases=records,seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=True,product_shell_adopted=False)
    (F/'q9-shear-gradient-compatibility-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_criteria_pass=True,cases=records,seconds=result['seconds'])))
if __name__=='__main__':main()
