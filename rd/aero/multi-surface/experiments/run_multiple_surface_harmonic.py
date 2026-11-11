import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import sys,time,json,numpy as np
from scipy.linalg import solve
from uvlm_frequency_delta import baseline,delta
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/multiple-surface-harmonic';side=sys.argv[1];OUT.mkdir(exist_ok=True);start=time.monotonic()
def transform(g,d,tail=False):
    g={k:(v.copy() if isinstance(v,np.ndarray) else v) for k,v in g.items()}
    for key in ['X','P','left','right']:
        p=g[key];y=p[...,1].copy()
        if tail:p[...,0]=(p[...,0]-abs(y)*np.tan(np.deg2rad(35)))*2/3+abs(.625*y)*np.tan(np.deg2rad(20))+5;p[...,1]*=.625;p[...,2]=.6
        else:p[...,2]+=abs(y)*np.tan(np.deg2rad(d))
    if tail:
        if 'area' in g:g['area']*=2/3*.625
        g['normal'][:]=[0,0,1]
    else:
        a=np.deg2rad(d);g['normal']=np.column_stack([np.zeros(len(g['P'])),-np.sign(g['P'][:,1])*np.sin(a),np.full(len(g['P']),np.cos(a))])
        if 'area' in g:g['area']/=np.cos(a)
    return g
for d in [0,5]:
    for m,n in [(2,8),(4,16),(8,32)]:
        if side=='candidate':
            from uvlm_steady_reference import geometry
        else:
            from ptera_frequency_reference import geometry
        wings=[transform(geometry(m,n),d),transform(geometry(m,n),0,True)];nb=m*n;P=np.concatenate([g['P'] for g in wings]);N=np.concatenate([g['normal'] for g in wings]);A=[];W=[];omega=8/1.114285714285714;rhs=np.zeros((2*nb,4),complex)
        for j,g in enumerate(wings):
            allg=dict(g,P=P,normal=N);A.append(baseline(allg,side,1e-5,128.));W.append(delta(allg,side,omega,1e-5,128.,16));idx=slice(j*nb,(j+1)*nb);axis=.375 if j==0 else 5.25;rhs[idx,2*j]=1j*omega*g['normal'][:,2];rhs[idx,2*j+1]=-(20+1j*omega*(g['P'][:,0]-axis))*g['normal'][:,2]
        A=np.concatenate(A,axis=1).astype(complex)
        for j in range(2):A[:,(j+1)*nb-n:(j+1)*nb]+=W[j]
        gamma=solve(A,rhs);res=np.linalg.norm(A@gamma-rhs,axis=0)/np.linalg.norm(rhs,axis=0);assert max(res)<1e-11;dual=[];native=[]
        for j,g in enumerate(wings):
            G=gamma[j*nb:(j+1)*nb];horse=G.reshape(m,n,4)-np.concatenate([np.zeros((1,n,4)),G.reshape(m,n,4)[:-1]],axis=0);fc=1.2*horse.reshape(nb,4)[:,None,:]*np.cross([20.,0,0],g['right']-g['left'])[:,:,None];F=fc.sum(axis=0);M=np.cross(((g['left']+g['right'])/2-[.375,0,0])[:,:,None],fc,axisa=1,axisb=1,axisc=1).sum(axis=0);D=np.r_[F,M];nativeD=D.copy();X=g['X'];quads=np.stack([X[:-1,:-1],X[:-1,1:],X[1:,1:],X[1:,:-1]],axis=2).reshape(-1,4,3)
            for p,q in enumerate(quads):
                for ids in [[0,2,1],[0,3,2]]:
                    tri=q[ids];area=np.cross(tri[1]-tri[0],tri[2]-tri[0])/2;f=1.2j*omega*area[:,None]*G[p];D[:3]+=f;D[3:]+=np.cross((tri.mean(axis=0)-[.375,0,0])[:,None],f,axisa=0,axisb=0,axisc=0)
            if side=='reference':
                for p in range(nb):
                    f=1.2j*omega*g['area'][p]*g['normal'][p,:,None]*G[p];nativeD[:3]+=f;nativeD[3:]+=np.cross((g['P'][p]-[.375,0,0])[:,None],f,axisa=0,axisb=0,axisc=0)
            else:nativeD=D.copy()
            dual.append(D);native.append(nativeD)
        dest=OUT/f'{side}-d{d}-m{m}-n{n}.npz';assert not dest.exists();np.savez_compressed(dest,gamma=gamma,dual=np.array(dual),native=np.array(native),observation_difference=np.array(dual)-np.array(native),A=A,rhs=rhs);print(json.dumps(dict(side=side,d=d,m=m,n=n,residuals=res.tolist(),seconds=time.monotonic()-start)),flush=True);assert time.monotonic()-start<180
