"""Same physical operator: impose declared supports before all-column QR.

The scalar Zs contraction commutes with the Cartesian strain expansion.
Also evaluate only the position coefficients needed by physical load work.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import svd
import pg_physical_pde as pde
f=pde.f;m=pde.m;paired=pde.paired
original_components=pde.components

def components(a,cell,q,p):
    fs=paired.base.surface(a['n'],cell,p,a['geo']);R,D,C=[v for v,_,_ in fs];v,_,_=paired.pool(q,p,a['n'],a['geo']['metric']);e=a['geo']['enrich'][0];nr,nd,nc=a['nr'],a['nd'],a['nc'];T=e['V']/e['s'];er=(v[0]-R@e['coeff'][:nr])@T;ed=(v[1]-D@e['coeff'][nr:])@T;out=np.zeros((3,len(p),a['ns']));out[0,:,:nr]=R;out[1,:,nr:nr+nd]=D;out[2,:,nr+nd:nr+nd+nc]=C;out[0,:,nr+nd+nc:]=er;out[1,:,nr+nd+nc:]=ed;return out

def prepare(k,n,amp,h,order=13):
    geo=f.c.prepare_geometry(k,n,amp,order);fs=paired.base.surface(n,0,np.zeros((1,2)),geo);nr,nd,nc=[v.shape[1] for v,_,_ in fs];ns=nr+nd+nc+geo['enrich'][0]['rank'];a=dict(k=k,n=n,amp=amp,h=h,geo=geo,nr=nr,nd=nd,nc=nc,ns=ns);constraints=[]
    for cell,q in enumerate(geo['refs']):
        if cell%n==0:
            edge=np.array([[-1.,t,z] for t in np.linspace(-1,1,33) for z in [-1.,0.,1.]]);constraints.append(f.shape_values(a,cell,q,edge))
    C=np.vstack(constraints);_,s,V=svd(C,full_matrices=True,check_finite=False);ranks=[int(sum(s>tol)) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('support_rank:'+str(ranks))
    Zs=V[ranks[0]:].T;Z=np.kron(Zs,np.eye(3));rows=[];commutation=0.;load_values=0.
    for cell,q in enumerate(geo['refs']):
        probe=np.array([[-.71,.33,.43],[.47,-.29,-.61],[.23,.17,.91]]);_,G,B,_,_=f.c.kinematics(n,cell,q,h,probe,geo);Gz=np.einsum('pni,nj->pji',G,Zs,optimize=True);Bz=np.einsum('kab,pnb->pkna',m.MB,Gz).reshape(len(probe),6,3*Zs.shape[1]);commutation=max(commutation,float(np.linalg.norm(Bz-B@Z)/max(1.,np.linalg.norm(B@Z))));aa=components(a,cell,q,probe[:,:2]);bb=original_components(a,cell,q,probe[:,:2]);load_values=max(load_values,float(np.linalg.norm(aa-bb)/max(1.,np.linalg.norm(bb))))
        p,w=f.c.volume_rule(q,n,order,7);R=None
        for pp,ww in f.chunks(p,w):
            _,G,_,det,_=f.c.kinematics(n,cell,q,h,pp,geo);Gz=np.einsum('pni,nj->pji',G,Zs,optimize=True);Bz=np.einsum('kab,pnb->pkna',m.MB,Gz).reshape(len(pp),6,3*Zs.shape[1]);A=(np.einsum('ij,pjn->pin',m.CH,Bz)*np.sqrt(ww*det)[:,None,None]/h**1.5).reshape(-1,Z.shape[1]);R=f.old.qrtools.rfactor(np.asfortranarray(A if R is None else np.vstack((R,A))))
        rows.append(R)
    if max(commutation,load_values)>1e-12:raise ValueError('operator_or_load_equivalence')
    A=np.vstack(rows);scale=1/np.linalg.norm(A,axis=0);R=f.old.qrtools.rfactor(np.asfortranarray(A*scale));sv=svd(R,compute_uv=False,check_finite=False);fr=[int(sum(sv>tol)) for tol in [1e-9,1e-10,1e-11]]
    if fr!=[Z.shape[1]]*3:raise ValueError('stiffness_rank:'+str(fr))
    a.update(A_free=A,Zs=Zs,Z=Z,scale=scale,R=R[:Z.shape[1]],rcond=f.old.qrtools.rcond(R[:Z.shape[1]]),support_ranks=ranks,commutation_error=commutation,load_value_error=load_values);return a

pde.prepare=prepare;pde.components=components;f.BATCH=2048
