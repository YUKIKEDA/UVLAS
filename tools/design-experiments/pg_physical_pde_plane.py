"""Evaluate the same surface functions once for repeated thickness points."""
from pathlib import Path
import hashlib,json
import numpy as np
import pg_physical_pde_cached as cached
pde=cached.pde;f=pde.f;m=pde.m;paired=pde.paired
_kinematics=f.c.kinematics;_prepare=pde.prepare

def kinematics(n,cell,q,h,p,geo):
    plane,back=np.unique(p[:,:2],axis=0,return_inverse=True);fs=paired.base.surface(n,cell,plane,geo);raw=paired.pool(q,plane,n,geo['metric']);e=geo['enrich'][0];nr,nd,nc=[v.shape[1] for v,_,_ in fs];ns=nr+nd+nc+len(e['s']);T=e['V']/e['s'];parts=[]
    for deriv,value in enumerate(raw):
        R,D,C=[v[deriv] for v in fs];er=(value[0]-R@e['coeff'][:nr])@T;ed=(value[1]-D@e['coeff'][nr:])@T;part=np.zeros((3,len(plane),ns));part[0,:,:nr]=R;part[1,:,nr:nr+nd]=D;part[2,:,nr+nd:nr+nd+nc]=C;part[0,:,nr+nd+nc:]=er;part[1,:,nr+nd+nc:]=ed;parts.append(part[:,back])
    z=h*p[:,2]/2;N=parts[0][0]+z[:,None]*parts[0][1]+z[:,None]**2*parts[0][2]/2;G=np.empty((len(p),ns,3))
    for axis in [0,1]:G[:,:,axis]=parts[axis+1][0]+z[:,None]*parts[axis+1][1]+z[:,None]**2*parts[axis+1][2]/2
    G[:,:,2]=h/2*(parts[0][1]+z[:,None]*parts[0][2]);X,J,det,inv,_=f.f.base.geometry(q,h,p);G=np.einsum('pni,pij->pnj',G,inv,optimize=True);B=np.einsum('kij,pnj->pkni',m.MB,G).reshape(len(p),6,3*ns);return N,G,B,det,X

def prepare(*args,**kwargs):
    a=_prepare(*args,**kwargs);maximum=0.
    for cell,q in enumerate(a['geo']['refs']):
        p=np.array([[x,y,z] for x,y in [(-.71,.33),(.47,-.29),(.23,.17)] for z in [-1.,-.31,.47,1.]])
        actual=kinematics(a['n'],cell,q,a['h'],p,a['geo']);expected=_kinematics(a['n'],cell,q,a['h'],p,a['geo'])
        for x,y in zip(actual[:3],expected[:3]):maximum=max(maximum,float(np.linalg.norm(x-y)/max(1.,np.linalg.norm(y))))
    assert maximum<=1e-12;a['plane_kernel_error']=maximum
    parent=cached.ROOT/a['preparation_manifest'];meta=parent.with_name(parent.stem+'-plane.json');data=dict(parent_preparation=a['preparation_manifest'],plane_kernel_error=maximum,sources={parent.relative_to(cached.ROOT).as_posix():hashlib.sha256(parent.read_bytes()).hexdigest(),Path(__file__).relative_to(cached.ROOT).as_posix():hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    if meta.exists():assert json.loads(meta.read_text())==data
    else:meta.write_text(json.dumps(data,indent=2),encoding='utf-8')
    a['preparation_manifest']=meta.relative_to(cached.ROOT).as_posix();return a

f.c.kinematics=kinematics;pde.prepare=prepare
