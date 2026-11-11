import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np,json,time
from numpy.polynomial.legendre import leggauss
from uvlm_frequency_fast import bound
from uvlm_steady_reference import geometry
from probe_finite_wing_drag import far_drag
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/finite-wing-2-3'
def main():
    for m,n in [(4,24),(16,96)]:
        own=np.load(OUT/f'candidate-d0-m{m}-n{n}.npz');g=geometry(m,n);X=g['X'].copy();X[-1,:,0]+=1024.;left,right=g['left'],g['right'];dl=right-left;normal=np.column_stack([np.zeros(len(dl)),-dl[:,2],dl[:,1]]);records=[]
        for order in [1,2,4,8,16]:
            start=time.monotonic();z,w=leggauss(order);total=0.
            for t,weight in zip((z+1)/2,w/2):
                p=(1-t)*left+t*right;A=bound(p,normal,X,1.114285714285714e-5);total-=weight*1.2*own['horse']@(A@own['gamma'])
            records.append(dict(order=order,near_drag=float(total),seconds=time.monotonic()-start));print(m,n,records[-1],flush=True)
        (OUT/f'drag-quadrature-m{m}-n{n}.json').write_text(json.dumps(dict(m=m,n=n,records=records,regularized_far_drag=far_drag(g['X'][-1],own['gamma'][-n:],1.114285714285714e-5)),indent=2),encoding='utf-8')
    records=[]
    for eps in [1e-2,5e-3,2.5e-3]:
        for n in [256,512,1024,2048,4096]:
            y=-np.cos(np.linspace(0,np.pi,n+1));mu=np.sqrt(1-((y[:-1]+y[1:])/2)**2);te=np.column_stack([np.zeros(n+1),y,np.zeros(n+1)]);D=far_drag(te,mu,eps);records.append(dict(eps=eps,n=n,drag=float(D),ideal_drag=float(np.pi*1.2/8),relative_difference=float(D/(np.pi*1.2/8)-1)))
    (OUT/'drag-elliptic-operator.json').write_text(json.dumps(dict(records=records),indent=2),encoding='utf-8');print(records[-1],flush=True)
if __name__=='__main__':main()
