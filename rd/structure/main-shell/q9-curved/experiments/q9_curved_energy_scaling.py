"""Constitutive-energy weighted shear for shared approximation; research only."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import itertools
import numpy as np
from scipy.linalg import cho_factor,cho_solve,solve_triangular,svd
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
import check_q9_curved_bounded_lift as core

from q9_curved_best_approximation import ReducedProblem
from audit_q9_curved_manufactured_lift import values as analytic_values

D6=np.diag([1.,1.,1.,np.sqrt(2.),np.sqrt(2.),np.sqrt(2.)])
CM=D6@core.shapes.elastic@D6;CI=np.linalg.inv(CM);MU=core.shapes.mu

def energy_rows(H):
    D=np.linalg.inv(H@CI@H.transpose(0,2,1));metric=D/MU;chol=np.linalg.cholesky(metric);lift=CI@H.transpose(0,2,1)@D;P=lift@H;I=np.eye(6)[None];cross=(I-P).transpose(0,2,1)@CM@P;reconstructed=H.transpose(0,2,1)@D@H+(I-P).transpose(0,2,1)@CM@(I-P);error=max(float(np.max(np.linalg.norm(cross,axis=(1,2))))/np.linalg.norm(CM),float(np.max(np.linalg.norm(reconstructed-CM,axis=(1,2))))/np.linalg.norm(CM));ev=np.linalg.eigvalsh(metric)
    return chol.transpose(0,2,1)@H,error,float(ev.min()),float(ev.max())

def assemble(name,n,h,observation=9,position_order=9):
    xy,cells,coords,_,_,_=core.geometry.make_mesh(name,n);_,reference,analytic=core.fields(name,xy);side=2*n+1;interior=np.array([j*side+i for j in range(1,side-1) for i in range(1,side-1)]);index={node:i for i,node in enumerate(interior)}
    points,w=core.quadrature(observation);nr=2*len(points);A=np.zeros((nr*n*n,6*len(xy)));W=np.zeros((len(A),3*len(interior)));M=np.zeros((len(xy),len(xy)));M0=np.zeros_like(M);load=np.zeros((3*len(xy),3));norm=np.zeros(3);volume=0.;area=0.;area0=0.
    hp,hw=leggauss(5);mp,mw=leggauss(position_order);energy_error=0.;metric_min=np.inf;metric_max=0.;edge_valid=True;edge_condition=0.
    for cell,(_,_,ids,_) in enumerate(cells):
        q=coords[ids].reshape(81);state=core.method.prepare(q,h,0.,9);edge_valid=edge_valid and all(len(set(e['ranks']))==1 and e['condition']<=1e8 for e in state['edges'].values());edge_condition=max(edge_condition,max(e['condition'] for e in state['edges'].values()));J,E,B=core.method.evaluate(state,points);a=J[:,:,:2];dz=J[:,:,2]*2/h;G=a.transpose(0,2,1)@a;row=2*np.einsum('pia,kij,pj->pak',a,core.method.MB,dz);white=np.linalg.solve(np.linalg.cholesky(G),row)/np.linalg.norm(dz,axis=1)[:,None,None];white,err,lo,hi=energy_rows(white);energy_error=max(energy_error,err);metric_min=min(metric_min,lo);metric_max=max(metric_max,hi);weight=w*np.linalg.det(J);volume+=weight.sum();C=(np.einsum('pak,pkn->pan',white,B[:,:,:81])*np.sqrt(weight)[:,None,None]).reshape(nr,81);sl=slice(nr*cell,nr*(cell+1))
        for local,node in enumerate(ids):
            A[sl,6*node:6*node+3]+=C[:,9*local+3:9*local+6];A[sl,6*node+3:6*node+6]+=C[:,9*local+6:9*local+9]/core.LENGTH
            if node in index:W[sl,3*index[node]:3*index[node]+3]+=C[:,9*local:9*local+3]
        H=np.zeros((9,9));Q=np.zeros((9,9));f=np.zeros((9,3,3))
        for i,j in itertools.product(range(5),repeat=2):
            N,dx,dy=core.shapes.shape(hp[i],hp[j]);g=np.column_stack((dx,dy));Js=coords[ids,0,:].T@g;metric=Js.T@Js;aw=hw[i]*hw[j]*np.sqrt(np.linalg.det(metric));area+=aw;H+=aw*(np.outer(N,N)+core.LENGTH**2*g@np.linalg.inv(metric)@g.T)
        for i,j in itertools.product(range(position_order),repeat=2):
            N,dx,dy=core.shapes.shape(mp[i],mp[j]);Js=coords[ids,0,:].T@np.column_stack((dx,dy));aw=mw[i]*mw[j]*np.sqrt(np.linalg.det(Js.T@Js));area0+=aw;Q+=aw*np.outer(N,N);x,y=N@xy[ids];normal=core.geometry.surface(name,x,y)[3]
            for k,(wa,wb) in enumerate(core.WAVES):
                value=np.sin(wa*np.pi*(x+.5))*np.sin(wb*np.pi*(y+.5))*normal;f[:,:,k]+=aw*N[:,None]*value[None,:];norm[k]+=aw*(value@value)
        M[np.ix_(ids,ids)]+=H;M0[np.ix_(ids,ids)]+=Q
        dofs=(3*ids[:,None]+np.arange(3)).reshape(-1);load[dofs]+=f.reshape(27,3)
    rids=(3*interior[:,None]+np.arange(3)).reshape(-1);A/=np.sqrt(volume);W/=np.sqrt(volume);M/=area;M0/=area0;load/=area0;norm/=area0;Gr=np.kron(M0[np.ix_(interior,interior)],np.eye(3))
    return dict(A=A,W=W,M=M,G=Gr,b=load[rids],norm=norm,reference=reference,interior=interior,xy=xy,analytic=analytic,volume=volume,energy_identity_error=energy_error,metric_min=metric_min,metric_max=metric_max,edge_valid=bool(edge_valid),edge_condition=edge_condition)


def analytic_point(name,x,y):
    v=analytic_values(name,1,1,x,y);dx=[a.imag/1e-25 for a in analytic_values(name,1,1,x+1e-25j,y)];dy=[a.imag/1e-25 for a in analytic_values(name,1,1,x,y+1e-25j)];r,d,c,u,dd,dc=v
    J0=np.column_stack((dx[0],dy[0],d));J1=np.column_stack((dx[1],dy[1],c));J2=np.column_stack((dx[2],dy[2],np.zeros(3)));Q0=np.column_stack((dx[3],dy[3],dd));Q1=np.column_stack((dx[4],dy[4],dc));Q2=np.column_stack((dx[5],dy[5],np.zeros(3)));Ji=np.linalg.inv(J0);F0=Q0@Ji;F1=Q1@Ji-F0@J1@Ji
    em=np.einsum('kij,ij->k',core.method.MB,(F0+F0.T)/2);bend=np.einsum('kij,ij->k',core.method.MB,(F1+F1.T)/2)
    z=1e-25j;Fz=(Q0+z*Q1+z*z*Q2/2)@np.linalg.inv(J0+z*J1+z*z*J2/2);other=np.einsum('kij,ij->k',core.method.MB,(Fz+Fz.T)/2).imag/1e-25;derivative=float(np.linalg.norm(other-bend)/max(np.linalg.norm(bend),1.))
    return em,bend,J0,J1,J2,derivative

def reference(name,order,thicknesses):
    p,w=leggauss(order);z,zw=leggauss(3);Im=0.;Ib=0.;volumes={str(h):0. for h in thicknesses};derivative=0.
    for i,j in itertools.product(range(order),repeat=2):
        em,k,J0,J1,J2,e=analytic_point(name,p[i]/2,p[j]/2);weight=w[i]*w[j]/4;measure=weight*np.linalg.det(J0);Im+=measure*(em@CM@em);Ib+=measure*(k@CM@k);derivative=max(derivative,e)
        for h in thicknesses:
            volumes[str(h)]+=sum(float(weight*ww*h/2*np.linalg.det(J0+(h*zz/2)*J1+(h*zz/2)**2*J2/2)) for zz,ww in zip(z,zw))
    return dict(membrane=float(Im),bending=float(Ib),analytic_volume=volumes,derivative_error=derivative)
