"""Separate near-force and exact regularized Trefftz-strip diagnostics."""
from pathlib import Path
import numpy as np,json,time
from uvlm_frequency_fast import bound
from uvlm_steady_reference import geometry
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/finite-wing-2-3'
def far_drag(te,mu,eps,rho=1.2):
    yz=te[:,1:];kappa=np.r_[0.,mu]-np.r_[mu,0.];total=0.
    for j in range(0,len(mu),128):
        end=min(j+128,len(mu));aa=np.sum((yz[j:end,None]-yz[None])**2,axis=-1)+eps*eps;bb=np.sum((yz[j+1:end+1,None]-yz[None])**2,axis=-1)+eps*eps;total+=mu[j:end]@(np.log(bb/aa)@kappa)
    return -rho*total/(8*np.pi)
def main():
    for d in [0,5]:
        for m,n,key in [(16,96,'m16-n96'),(16,192,'m16-n192'),(16,384,'span384')]:
            name=f'drag-d{d}-{key}';dest=OUT/(name+'.json')
            if dest.exists():continue
            start=time.monotonic();own=np.load(OUT/f'candidate-d{d}-{key}.npz');ref=np.load(OUT/f'reference-d{d}-{key}.npz');g=geometry(m,n,d);X=g['X'].copy();X[-1,:,0]+=1024.;p=(g['left']+g['right'])/2;dl=g['right']-g['left'];normal=np.column_stack([np.zeros(len(dl)),-dl[:,2],dl[:,1]]);A=bound(p,normal,X,1.114285714285714e-5);near=-1.2*own['horse']@(A@own['gamma']);angle=np.deg2rad(1.);U=20*np.array([np.cos(angle),0,np.sin(angle)]);known=(1.2*ref['horse'][:,None]*np.sin(angle)*np.cross(U,ref['right']-ref['left'])).sum(axis=0);r=json.loads((OUT/f'reference-d{d}-{key}.json').read_text());refnear=(np.array(r['raw_F'])-known)[0]/np.sin(angle)**2;far=far_drag(g['X'][-1],own['gamma'][-n:],1.114285714285714e-5)
            result=dict(dihedral=d,m=m,n=n,near_drag=float(near),independent_near_drag=float(refnear),near_difference=float(abs(near/refnear-1)),regularized_far_drag=float(far),far_near_difference=float(far/near-1),eps=1.114285714285714e-5,seconds=time.monotonic()-start);dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result),flush=True);del A
if __name__=='__main__':main()
