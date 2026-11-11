"""Joint assumed bulk strain and work-compatible weak position trace."""
import itertools
import numpy as np
import mortar_q4 as mortar
import tying_mid_projection_q4_v2 as tying
from mortar_q4 import broken_mesh,constraints,null_rows,rigid_scan,boundary,ranks
m=mortar.m;MB,C,CH=m.MB,m.C,m.CH

def linearized(q,h,points):
    order=int(np.sqrt(len(points)/5));_,_,B=m.raw(q,points,h*points[:,2]/2)
    return tying.project_fields(q,h,B,order=order,order_z=5)[2]

def mid_project(q,p,B0,independent=False):
    state=tying.make_basis(q,int(np.sqrt(len(p))),independent)
    return tying.project_mid(state,0.,B0)[0]

def raw_nonlinear(q,current,h,v1,v2,order=9,order_z=5,independent=False):
    plane,wp=m.quadrature(order,2);z,wz=np.polynomial.legendre.leggauss(order_z)
    points=np.array([[s,t,c] for (s,t),c in itertools.product(plane,z)]);w=np.array([a*b for a,b in itertools.product(wp,wz)])
    _,g,_=m.tables(points,h*points[:,2]/2,independent);J=np.einsum('na,pnb->pab',q.reshape(75,3),g);G=g@np.linalg.inv(J)
    weight=w*np.linalg.det(J)*h/2;F=np.eye(3)+np.einsum('na,pnb->pab',(current-q).reshape(75,3),G)
    E=np.einsum('kij,pij->pk',m.MB,(F.transpose(0,2,1)@F-np.eye(3))/2);B=np.einsum('kij,pai,pnj->pkna',m.MB,F,G).reshape(-1,6,225)
    d1=np.einsum('na,pnb->pab',v1.reshape(75,3),G);d2=np.einsum('na,pnb->pab',v2.reshape(75,3),G)
    second=np.einsum('kij,pij->pk',m.MB,(d1.transpose(0,2,1)@d2+d2.transpose(0,2,1)@d1)/2)
    state=tying.make_basis(q,order,independent) if independent else tying.prepare(q,order)
    def project(v):return tying.project_fields(q,h,v,order=order,order_z=order_z,state=state)[2]
    Er=project(E[:,:,None])[:,:,0];Br=project(B);second_r=project(second[:,:,None])[:,:,0];stress=Er@m.C.T
    energy=np.einsum('p,pk,pk->',weight,Er,stress)/2;force=np.einsum('p,pki,pk->i',weight,Br,stress)
    tangent=np.einsum('p,pk,kl,pl->',weight,Br@v1,m.C,Br@v2)+np.einsum('p,pk,pk->',weight,stress,second_r)
    return energy,force,tangent

def patch(mesh,h,degree,F,independent=False):
    G,meta=constraints(mesh,h,degree,11 if independent else 9,7 if independent else 5,independent)
    E=np.einsum('kij,ij->k',m.MB,(F.T@F-np.eye(3))/2);S=np.einsum('kij,k->ij',m.MB,m.C@E);P=F@S
    internal=np.zeros_like(mesh['ref']);outer=np.zeros_like(internal)
    for c,ids in enumerate(mesh['cells']):
        q=mesh['ref'][ids].reshape(225);current=(mesh['ref'][ids]@F.T).reshape(225)
        _,fe,_=raw_nonlinear(q,current,h,np.zeros(225),np.zeros(225),11 if independent else 9,7 if independent else 5,independent);internal[ids]+=fe.reshape(75,3)
        j,i=divmod(c,mesh['n']);faces=[(2,-1),(2,1)]
        if i==0:faces.append((0,-1))
        if i==mesh['n']-1:faces.append((0,1))
        if j==0:faces.append((1,-1))
        if j==mesh['n']-1:faces.append((1,1))
        outer[ids]+=boundary(q,h,P,faces,11 if independent else 9,independent)
    lm=np.zeros((len(G),3))
    for i,entry in enumerate(meta):
        if entry['kind']=='piola':lm[i]=-h*P[:,entry['component']]
    return dict(residual=float(np.linalg.norm(internal+G.T@lm-outer)/max(1e-30,np.linalg.norm(outer))),
                gap=float(np.linalg.norm(G@(mesh['ref']@F.T))/max(1.,np.linalg.norm(mesh['ref']))),
                force_norm=float(np.linalg.norm(internal)),reaction_norm=float(np.linalg.norm(G.T@lm)))

def shared(mesh,k,degree,order=9,independent=False):
    a=tying.assemble(mesh,k,order,'conservative');G,_=constraints(mesh,0.,degree,order,7 if independent else 5,independent)
    Zrd,rd=null_rows(G[:,a['free_rd']],independent);Zc,c=null_rows(G[:,a['free_c']],independent)
    return a,np.kron(Zrd,np.eye(3)),np.kron(Zc,np.eye(3)),dict(rd=rd,c=c)
