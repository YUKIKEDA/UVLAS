"""T-S36-v1 research kernel: true reference covariant sample interpolation."""
import itertools
import numpy as np
from numpy.polynomial.legendre import leggauss
import check_q9_distorted_plate as base
import check_q9_curved_load as geometry
import check_q9_independent_load as material
import check_q9_shear_projection_gate as baseline
A=1/np.sqrt(3.)
AXES=[([-A,A],[-1.,0.,1.]),([-1.,0.,1.],[-A,A])]
def cardinal(x,nodes,index):
    value=np.ones_like(np.asarray(x,dtype=float))
    for j,t in enumerate(nodes):
        if j!=index:value*= (x-t)/(nodes[index]-t)
    return value
def interpolation(points,thickness=True):
    points=np.asarray(points);count=36 if thickness else 12;S=np.zeros((len(points),2,count));column=0
    for component,(xs,ys) in enumerate(AXES):
        for i,j in itertools.product(range(len(xs)),range(len(ys))):
            if thickness:
                for k in range(3):
                    S[:,component,column]=cardinal(points[:,0],xs,i)*cardinal(points[:,1],ys,j)*cardinal(points[:,2],[-1.,0.,1.],k);column+=1
            else:S[:,component,column]=cardinal(points[:,0],xs,i)*cardinal(points[:,1],ys,j);column+=1
    return S
def sampling_points():
    points=[];components=[]
    for component,(xs,ys) in enumerate(AXES):
        for x,y,z in itertools.product(xs,ys,[-1.,0.,1.]):points.append([x,y,z]);components.append(component)
    return np.array(points),np.array(components)
def kinematics(q,h,points):
    scalar,g=geometry.shape_tables(points,h);jac=np.einsum('ia,pib->pab',q.reshape(27,3),g)
    if np.min(np.linalg.det(jac))<=0:raise ValueError('nonpositive_reference_jacobian')
    grad=np.linalg.solve(jac.transpose(0,2,1),g.transpose(0,2,1)).transpose(0,2,1)
    B=np.zeros((len(points),2,81));B[:,0,0::3]=grad[:,:,2];B[:,0,2::3]=grad[:,:,0];B[:,1,1::3]=grad[:,:,2];B[:,1,2::3]=grad[:,:,1]
    return jac,B
def sampled(q,h,points,components,chain=None):
    jac,B=kinematics(q,h,points);rows=np.array([jac[i,:2,c]@B[i] for i,c in enumerate(components)])
    if chain is not None:rows*=np.asarray(chain)[:,None]
    return rows
def matrices(q,h,order):
    original=baseline.matrices(q,h,order);p,w=leggauss(order);z,wz=leggauss(3);points=np.array(list(itertools.product(p,p,z)));weights=np.array([a*b*c for a,b,c in itertools.product(w,w,wz)])
    jac,_=kinematics(q,h,points);S=interpolation(points);V=np.linalg.solve(jac[:,:2,:2].transpose(0,2,1),S);measure=material.mu*weights*np.linalg.det(jac)
    M=np.einsum('p,pai,paj->ij',measure,V,V);factor=baseline.factor(M)
    ties,components=sampling_points();T=sampled(q,h,ties,components);K=original['Kn']+T.T@M@T
    planar=np.array([[x,y,0.] for x,y in itertools.product(p,p)]);pw=np.array([a*b for a,b in itertools.product(w,w)]);j0,_=kinematics(q,h,planar);S0=interpolation(planar,False);V0=np.linalg.solve(j0[:,:2,:2].transpose(0,2,1),S0)
    M0=np.einsum('p,pai,paj->ij',material.mu*pw*np.linalg.det(j0[:,:2,:2]),V0,V0);T0=T[1::3]
    E=np.repeat(np.eye(12),3,axis=0);metric_difference=float(np.linalg.norm(E.T@M@E/h-M0)/np.linalg.norm(M0))
    return dict(M=M,T=T,Kn=original['Kn'],Kfull=original['Kfull'],Ktie=K,factor=factor,M0=M0,T0=T0,C0=np.linalg.cholesky(M0).T@T0,leading_metric_difference=metric_difference)
