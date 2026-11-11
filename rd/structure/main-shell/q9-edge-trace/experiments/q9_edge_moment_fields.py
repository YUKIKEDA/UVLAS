"""Geometry-only G-EM1 field/adjoint research kernel, not a production API."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import itertools
import numpy as np
from numpy.polynomial.legendre import leggauss
import check_q9_curved_strain_screen as raw
from q9_edge_moment import MB

EDGES=list(itertools.product([0,1],[-1,1]))
M2=np.array([[[1.,0.],[0.,0.]],[[0.,0.],[0.,1.]],[[0.,1/np.sqrt(2.)],[1/np.sqrt(2.),0.]]])

def quadrature(n):
    p,w=leggauss(n)
    return np.array(list(itertools.product(p,p,p))),np.prod(np.array(list(itertools.product(w,w,w))),axis=1)

def congruence(A):
    return np.einsum('aij,pki,bkl,plj->pab',MB,A,MB,A,optimize=True)

def freeze(x):
    if isinstance(x,np.ndarray):x.setflags(write=False)
    elif isinstance(x,dict):
        for v in x.values():freeze(v)
    elif isinstance(x,(list,tuple)):
        for v in x:freeze(v)
    return x

class ReferenceReduction:
    def __init__(self,q,width,offset=0.,order=5):
        self.q=np.array(q,copy=True);self.width=width;self.offset=offset
        self.points,w=quadrature(order)
        self.ties,components=raw.trace.sampling_points()
        self.edge_points=[];self.edge_Linv=[];self.edge_U=[];self.pinv=[]
        J=self.jac(self.points);self.measure=w*np.linalg.det(J)
        if np.any(self.measure<=0):raise ValueError('invalid_reference_orientation')
        self.av=self.measure/self.measure.sum()
        self.Jinv=np.linalg.inv(J);self.tie_Jinv=np.linalg.inv(self.jac(self.ties))
        H=congruence(self.jac(self.ties))
        self.tie_map=H[np.arange(36),[4 if c==0 else 3 for c in components]]/np.sqrt(2.)
        for axis,sign in EDGES:
            p,gw=leggauss(5)
            ep=np.array([[x,sign,z] if axis==0 else [sign,x,z] for x,z in itertools.product(p,p)])
            A=self.jac(ep)[:,:,[axis,2]];G=np.einsum('pia,pib->pab',A,A)
            Li=np.linalg.inv(np.linalg.cholesky(G));U=A@Li.transpose(0,2,1)
            D=np.einsum('aij,pki,bkl,plj->pab',M2,U,MB,U,optimize=True).reshape(75,6)
            ew=np.array([a*b for a,b in itertools.product(gw,gw)])*np.sqrt(np.linalg.det(G))
            sw=np.repeat(np.sqrt(ew/ew.sum()),3);u,s,v=np.linalg.svd(sw[:,None]*D,full_matrices=False)
            ranks=[int(np.count_nonzero(s>tol*s[0])) for tol in [1e-9,1e-10,1e-11]]
            if len(set(ranks))!=1 or ranks[1]<3 or s[0]/s[ranks[1]-1]>1e8:raise ValueError('invalid_edge_rank')
            r=ranks[1];self.pinv.append((v[:r].T/s[:r])@(u[:,:r].T*sw[None,:]))
            self.edge_points.append(ep);self.edge_Linv.append(Li);self.edge_U.append(U)
        self.A,self.T,self.L,self.G=self.maps(self.points)
        b=(1-self.points[:,0]**2)*(1-self.points[:,1]**2);self.bubble=b/(self.av@b)
        for value in vars(self).values():freeze(value)

    def jac(self,points,q=None):
        _,g=raw.tables(points,self.width,self.offset)
        return np.einsum('ni,pnj->pij',(self.q if q is None else q).reshape(27,3),g)

    def defect(self,points):
        H=congruence(self.jac(points))
        return H[:,[4,3],:]/np.sqrt(2.)-np.einsum('pai,ik->pak',raw.trace.interpolation(points),self.tie_map)

    def maps(self,points):
        H=congruence(self.jac(points));Hi=congruence(np.linalg.inv(self.jac(points)))
        kept=H.copy();kept[:,[4,3]]=0
        A=Hi@kept;lift=Hi[:,:,[4,3]]*np.sqrt(2.)
        T=np.einsum('pka,pai->pki',lift,raw.trace.interpolation(points))
        L=lift@self.defect(points);G=[]
        for axis,sign in EDGES:
            ep=np.array(points,copy=True);ep[:,1-axis]=sign
            ge=lift[:,:,axis,None]*self.defect(ep)[:,axis,None,:]*(1+sign*points[:,1-axis,None,None])/2
            G.append(ge);L=L-ge
        return A,T,L,np.array(G)

    def validate(self,fields):
        for value,shape in zip(fields,[(len(self.points),6),(36,6),(4,25,3)]):
            if value.shape!=shape or not np.isrealobj(value) or not np.all(np.isfinite(value)):raise ValueError('invalid_field_values')

    def apply(self,fields):
        self.validate(fields)
        return self.kernel(fields)

    def kernel(self,fields):
        v,t,e=fields;mean=self.av@v;beta=np.einsum('pi,pi->p',self.tie_map,t)
        ce=np.array([P@x.reshape(75) for P,x in zip(self.pinv,e)])
        first=np.einsum('pij,pj->pi',self.A,v)+np.einsum('pki,i->pk',self.T,beta)+self.L@mean+np.einsum('epij,ej->pi',self.G,ce)
        moment=self.av@(first-v)
        return first-self.bubble[:,None]*moment

    def adjoint(self,y):
        if y.shape!=(len(self.points),6) or not np.all(np.isfinite(y)):raise ValueError('invalid_cotangent')
        moment_bar=-self.bubble@y;first_bar=y+self.av[:,None]*moment_bar
        v=np.einsum('pij,pi->pj',self.A,first_bar)-self.av[:,None]*moment_bar
        mean_bar=np.einsum('pij,pi->j',self.L,first_bar);v+=self.av[:,None]*mean_bar
        beta_bar=np.einsum('pki,pk->i',self.T,first_bar);t=self.tie_map*beta_bar[:,None]
        ce_bar=np.einsum('epij,pi->ej',self.G,first_bar)
        e=np.array([(P.T@x).reshape(25,3) for P,x in zip(self.pinv,ce_bar)])
        return v,t,e

    def deformation_maps(self):
        """F[p,i,j]=q[n,i] H[p,n,j]; edge uses A L0^-T."""
        result=[]
        for points,Ji in [(self.points,self.Jinv),(self.ties,self.tie_Jinv)]:
            _,g=raw.tables(points,self.width,self.offset);result.append(g@Ji)
        edges=[]
        for (axis,_),ep,Li in zip(EDGES,self.edge_points,self.edge_Linv):
            _,g=raw.tables(ep,self.width,self.offset);edges.append(g[:,:,[axis,2]]@Li.transpose(0,2,1))
        return result+[edges]

def green_sources(ref,q):
    out=[]
    for H in ref.deformation_maps()[:2]:
        F=np.einsum('ni,pnj->pij',q.reshape(27,3),H)
        E=(np.einsum('pki,pkj->pij',F,F)-np.eye(3))/2
        out.append(np.einsum('aij,pij->pa',MB,E))
    e=[]
    for H in ref.deformation_maps()[2]:
        A=np.einsum('ni,pnj->pij',q.reshape(27,3),H)
        e.append(np.einsum('aij,pij->pa',M2,(np.einsum('pki,pkj->pij',A,A)-np.eye(2))/2))
    return tuple(out+[np.array(e)])

def variations(ref,q,d):
    first=[];second=[]
    maps=ref.deformation_maps();entries=[(maps[0],MB),(maps[1],MB)]+[(H,M2) for H in maps[2]]
    for H,M in entries:
        F=np.einsum('ni,pnj->pij',q.reshape(27,3),H);D=np.einsum('ni,pnj->pij',d.reshape(27,3),H)
        first.append(np.einsum('aij,pki,pkj->pa',M,F,D))
        second.append((H,M,D))
    return (first[0],first[1],np.array(first[2:])),second

def pullback(ref,q,cot,d=None):
    maps=ref.deformation_maps();entries=[(maps[0],MB,cot[0]),(maps[1],MB,cot[1])]+[(H,M2,c) for H,c in zip(maps[2],cot[2])]
    force=np.zeros((27,3),dtype=np.result_type(q,*cot))
    for H,M,c in entries:
        F=np.einsum('ni,pnj->pij',(q if d is None else d).reshape(27,3),H)
        S=np.einsum('pa,aij->pij',c,M)
        force+=np.einsum('pij,pjk,pnk->ni',F,S,H,optimize=True)
    return force.reshape(81)

def energy_force(ref,q,initial,C):
    x=tuple(a-b for a,b in zip(green_sources(ref,q),initial));strain=ref.kernel(x)
    stress=strain@C.T;y=ref.measure[:,None]*stress
    return np.sum(strain*y)/2,pullback(ref,q,ref.adjoint(y))

def tangent_action(ref,q,initial,C,d):
    x=tuple(a-b for a,b in zip(green_sources(ref,q),initial));strain=ref.kernel(x)
    ds,_=variations(ref,q,d);reduced=ref.kernel(ds)
    material=pullback(ref,q,ref.adjoint(ref.measure[:,None]*(reduced@C.T)))
    geometric=pullback(ref,q,ref.adjoint(ref.measure[:,None]*(strain@C.T)),d=d)
    return material+geometric,material
