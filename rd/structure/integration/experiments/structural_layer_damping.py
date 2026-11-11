"""Layer/point-dependent reference-Mandel viscosity; base kernels stay frozen."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from structural_load_damping import shell_damping as constant_damping
from pg_ritz_mechanics import MB

def shell_damping(data,u,v,factor,tangent=True):
    factor=np.asarray(factor)
    if factor.ndim==2:return constant_damping(data,u,v,factor,tangent)
    G=data['G'];w=data['weights'];npnt=len(G)
    if factor.ndim!=3 or factor.shape[0]!=npnt or factor.shape[2]!=6 or not np.isfinite(factor).all():raise ValueError('invalid_pointwise_viscosity_factor')
    F=np.eye(3)+np.einsum('ia,pib->pab',u,G);V=np.einsum('ia,pib->pab',v,G);Edot=(F.transpose(0,2,1)@V+V.transpose(0,2,1)@F)/2;e=np.einsum('kab,pab->pk',MB,Edot);D=np.swapaxes(factor,1,2)@factor;s=np.einsum('pkl,pl->pk',D,e);S=np.einsum('pk,kab->pab',s,MB);P=F@S;force=np.einsum('p,pib,pab->ia',w,G,P,optimize=True);rate=np.einsum('prk,pk->pr',factor,e);dissipation=float(np.real(np.sum(w[:,None]*rate*rate)));Kq=Kv=None
    if tangent:
        ns=len(u);Kq=np.zeros((3*ns,3*ns),dtype=np.result_type(u,v,float));Kv=Kq.copy()
        for start in range(0,npnt,128):
            sl=slice(start,start+128);B=np.einsum('kij,pai,pnj->pkna',MB,F[sl],G[sl],optimize=True).reshape(len(G[sl]),6,3*ns);Bv=np.einsum('kij,pai,pnj->pkna',MB,V[sl],G[sl],optimize=True).reshape(len(G[sl]),6,3*ns);DB=D[sl]@B;DBv=D[sl]@Bv;Kv+=np.einsum('p,pki,pkj->ij',w[sl],B,DB,optimize=True);Kq+=np.einsum('p,pki,pkj->ij',w[sl],B,DBv,optimize=True)
        geo=np.einsum('p,pni,pij,pmj->nm',w,G,S,G,optimize=True);Kq+=np.kron(geo,np.eye(3))
    return force,Kq,Kv,dict(dissipation=dissipation,strain_rate=e,viscous_second_PK=S,F=F,V=V)
