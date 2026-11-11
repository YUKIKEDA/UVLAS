from pathlib import Path
import numpy as np,json,hashlib
from scipy.linalg import solve
from scipy.spatial import cKDTree
from uvlm_frequency_fast import bound
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/multiple-surface-2-3';CONFIG=[(4.,1.5,.6,35.,[0.,0.,0.]),(2.5,1.,.4,20.,[5.,0.,.6])]
def geometry(m,n,config):
    s,cr,ct,sweep,offset=config;t=(1-np.cos(np.linspace(0,np.pi,n//2+1)))/2;y=np.r_[s*(t-1),s*t[1:]];xi=(1-np.cos(np.linspace(0,np.pi,m+1)))/2;xv=np.r_[.75*xi[:-1]+.25*xi[1:],1.];xc=.25*xi[:-1]+.75*xi[1:]
    def p(x,y):
        x,y=np.broadcast_arrays(x,y);return np.stack([abs(y)*np.tan(np.deg2rad(sweep))+x*(cr+(ct-cr)*abs(y)/s),y,np.zeros_like(x)],axis=-1)+offset
    X=p(xv[:,None],y[None]);P=p(xc[:,None],((y[:-1]+y[1:])/2)[None]).reshape(-1,3);return dict(X=X,P=P,left=X[:-1,:-1].reshape(-1,3),right=X[:-1,1:].reshape(-1,3))
def wrench(horse,left,right):
    F=1.2*horse[:,None]*np.cross([20.,0,0],right-left);return np.array([F[:,2].sum(),np.cross((left+right)/2,F)[:,1].sum()])
def main():
    records=[];artifacts={}
    for m,n in [(2,12),(4,24),(8,48)]:
        gs=[geometry(m,n,c) for c in CONFIG];P=np.vstack([g['P'] for g in gs]);N=np.tile([0.,0,1.],(len(P),1));blocks=[]
        for g in gs:
            X=g['X'].copy();X[-1,:,0]+=1024.;blocks.append(bound(P,N,X,0.))
        A=np.hstack(blocks);gamma=solve(A,-20*N[:,2]);res=float(np.linalg.norm(A@gamma+20*N[:,2])/np.linalg.norm(20*N[:,2]));assert res<1e-11;own=[];single=[];refw=[];refs=[]
        raw=OUT/f'reference-coupled-m{m}-n{n}.npz';ref=np.load(raw);dist,index=cKDTree(ref['collocation']).query(P);assert max(dist)<1e-12
        for i,g in enumerate(gs):
            sl=slice(i*m*n,(i+1)*m*n);G=gamma[sl].reshape(m,n);h=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel();own.append(wrench(h,g['left'],g['right']));sr=solve(A[sl,sl],-20*N[sl,2]).reshape(m,n);hs=(sr-np.vstack([np.zeros((1,n)),sr[:-1]])).ravel();single.append(wrench(hs,g['left'],g['right']));refw.append(wrench(ref['horse'][index[sl]],g['left'],g['right']));label=['main','tail'][i];rr=OUT/f'reference-{label}-m{m}-n{n}.npz';r=np.load(rr);refs.append(wrench(r['horse'],r['left'],r['right']));artifacts[rr.relative_to(ROOT).as_posix()]=hashlib.sha256(rr.read_bytes()).hexdigest()
        own=np.array(own);refw=np.array(refw);single=np.array(single);refs=np.array(refs);error=abs((own-refw)/refw);interference=own-single;iref=refw-refs;ie=abs((interference-iref)/iref);defect=abs(interference/single);assert max(error.ravel())<.01 and max(ie.ravel())<.01 and max(defect.ravel())>.01;records.append(dict(m=m,n=n,original_residual=res,per_wing_wrenches=own.tolist(),independent_wrenches=refw.tolist(),relative_errors=error.tolist(),coupling_change=defect.tolist(),interference_relative_errors=ie.tolist()));artifacts[raw.relative_to(ROOT).as_posix()]=hashlib.sha256(raw.read_bytes()).hexdigest();print(records[-1],flush=True)
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'tools/design-experiments/run_multiple_surface_reference.py',ROOT/'docs/verification/uvlm-multiple-surface-reference-protocol.md']};result=dict(passed=True,scope='physical multiple-surface mutual operator and first-order wrenches, not continuum accuracy',records=records,sources=sources,artifacts=artifacts,uvlm_2_3_complete=False);dest=ROOT/'docs/verification/multiple-surface-independent-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8')
if __name__=='__main__':main()
