"""Two-node constant-strain SE(3) beam; local chart AD and consistent inertia."""
import math
import numpy as np
from scipy.linalg import expm
import beam_jet2 as ad

def skew(w):
    x,y,z=w;return np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]],dtype=object)
def vee(W):return np.array([W[2,1],W[0,2],W[1,0]],dtype=object)
def eye():return np.eye(3,dtype=object)
def trig(t):
    if abs(ad.value(t))<1e-4:
        A=ad.polynomial(t,[(-1)**k/math.factorial(2*k+1) for k in range(7)])
        B=ad.polynomial(t,[(-1)**k/math.factorial(2*k+2) for k in range(7)])
        C=ad.polynomial(t,[(-1)**k/math.factorial(2*k+3) for k in range(7)])
    else:
        a=ad.sqrt(t);A=ad.sin(a)/a;B=(1-ad.cos(a))/t;C=(1-A)/t
    return A,B,C
def exp_rotation(w):
    W=skew(w);A,B,C=trig(sum(a*a for a in w));return eye()+W*A+(W@W)*B
def left_jacobian(w):
    W=skew(w);A,B,C=trig(sum(a*a for a in w));return eye()+W*B+(W@W)*C
def log_rotation(R):
    c=(np.trace(R)-3)/2
    if abs(ad.value(c))<1e-4:
        coeff=[1.]
        for k in range(1,10):coeff.append(-k/(2*k+1)*coeff[-1])
        factor=ad.polynomial(c,coeff)
    else:
        cosine=1+c;angle=ad.acos(cosine);factor=angle/ad.sin(angle)
    return vee((R-R.T)/2)*factor
def inverse_left(w):
    t=sum(a*a for a in w);W=skew(w)
    if abs(ad.value(t))<1e-4:D=ad.polynomial(t,[1/12,1/720,1/30240,1/1209600,1/47900160])
    else:
        a=ad.sqrt(t);D=(1-(a/2)*(1+ad.cos(a))/ad.sin(a))/t
    return eye()-W*.5+(W@W)*D
def relative(r1,R1,r2,R2):
    w=log_rotation(R1.T@R2);v=inverse_left(w)@(R1.T@(r2-r1));return np.r_[v,w]
def float_array(a):return np.vectorize(ad.value)(np.asarray(a))
def pose_at(r1,R1,d,t):
    w=t*d[3:];v=t*d[:3];return r1+R1@left_jacobian(w)@v,R1@exp_rotation(w)
def update(r,R,delta):
    out=r+np.asarray(delta).reshape(2,6)[:,:3];rot=[]
    for i in range(2):rot.append(R[i]@expm(np.array(skew(delta[6*i+3:6*i+6]),float)))
    return out,np.array(rot)
def validate(r,R,L):
    if not np.isfinite(L) or L<=0:raise ValueError('nonpositive_reference_length')
    if np.shape(r)!=(2,3) or np.shape(R)!=(2,3,3) or not np.isfinite(r).all() or not np.isfinite(R).all():raise ValueError('invalid_pose_shape_or_finite')
    for Q in R:
        if np.linalg.norm(Q.T@Q-np.eye(3))>1e-10 or np.linalg.det(Q)<=0:raise ValueError('invalid_rotation')
    cosine=np.clip((np.trace(R[0].T@R[1])-1)/2,-1,1)
    if math.acos(cosine)>=math.pi-1e-4:raise ValueError('relative_rotation_log_branch')
def jets(r,R):
    rs=[];Rs=[]
    for i in range(2):
        rs.append(np.array([ad.Jet.variable(r[i,j],6*i+j) for j in range(3)],dtype=object))
        theta=np.array([ad.Jet.variable(0.,6*i+3+j) for j in range(3)],dtype=object);Rs.append(np.array(R[i],dtype=object)@exp_rotation(theta))
    return rs,Rs
def constitutive(K,L):
    if np.shape(K)!=(6,6) or not np.isfinite(K).all() or np.linalg.norm(K-K.T)>1e-12*np.linalg.norm(K):raise ValueError('invalid_section_stiffness')
    s=np.r_[np.ones(3),np.full(3,1/L)];C=K*s[:,None]*s[None,:]
    if np.any(np.diag(C)<=0):raise ValueError('nonpositive_section_stiffness')
    scale=1/np.sqrt(np.diag(C))
    if np.linalg.eigvalsh(C*scale[:,None]*scale[None,:]).min()<=1e-12:raise ValueError('nonpositive_section_stiffness')
def energy_from_poses(rs,Rs,L,K,reference,e0):
    d=relative(rs[0],Rs[0],rs[1],Rs[1]);e=d/L-reference-e0;U=L*(e@(K@e))/2
    return U,e,d
def evaluate(r,R,L,K,reference,e0=None):
    validate(r,R,L);constitutive(K,L);e0=np.zeros(6) if e0 is None else np.asarray(e0)
    if np.shape(reference)!=(6,) or np.shape(e0)!=(6,) or not np.isfinite(reference).all() or not np.isfinite(e0).all():raise ValueError('invalid_reference_or_initial_strain')
    rs,Rs=jets(r,R);U,e,d=energy_from_poses(rs,Rs,L,K,reference,e0)
    B=np.array([a.g for a in e]);return dict(U=U.v,f=U.g,K=U.H,strain=float_array(e),d=float_array(d),B=B)
def spatial_mass(m,c,Jc):
    c=np.asarray(c);Jc=np.asarray(Jc)
    if not np.isfinite(m) or m<=0 or c.shape!=(3,) or Jc.shape!=(3,3) or not np.isfinite(c).all() or not np.isfinite(Jc).all():raise ValueError('invalid_mass_data')
    vals=np.linalg.eigvalsh(Jc)
    if np.linalg.norm(Jc-Jc.T)>1e-12*np.linalg.norm(Jc) or vals[0]<=0 or vals[-1]>sum(vals[:-1])*(1+1e-12):raise ValueError('nonphysical_section_inertia')
    C=np.array(skew(c),float);Jo=Jc+m*((c@c)*np.eye(3)-np.outer(c,c));return np.block([[m*np.eye(3),-m*C],[m*C,Jo]])
def point_jet(r,R,t,eta=None):
    rs,Rs=jets(r,R);d=relative(rs[0],Rs[0],rs[1],Rs[1]);x,Q=pose_at(rs[0],Rs[0],d,t)
    if eta is not None:x=x+Q@np.asarray(eta,dtype=object)
    return x,Q
def inertia(r,R,L,mass,order=12):
    validate(r,R,L);nodes,weights=np.polynomial.legendre.leggauss(order);M=np.zeros((12,12));dM=np.zeros((12,12,12));rs,Rs=jets(r,R);d=relative(rs[0],Rs[0],rs[1],Rs[1])
    for t,w in zip((nodes+1)/2,weights*L/2):
        x,Q=pose_at(rs[0],Rs[0],d,t);vals=float_array(Q);grad=np.array([[a.g for a in row] for row in Q]);hess=np.array([[a.H for a in row] for row in Q])
        gx=np.array([a.g for a in x]);hx=np.array([a.H for a in x]);V=vals.T@gx;W=np.empty((3,12));dV=np.empty((3,12,12));dW=np.empty_like(dV)
        for j in range(12):
            W[:,j]=np.array(vee(vals.T@grad[:,:,j]),float)
            for k in range(12):
                dV[:,j,k]=grad[:,:,k].T@gx[:,j]+vals.T@hx[:,j,k]
                dW[:,j,k]=np.array(vee(grad[:,:,k].T@grad[:,:,j]+vals.T@hess[:,:,j,k]),float)
        B=np.vstack((V,W));DB=np.concatenate((dV,dW),axis=0);M+=w*(B.T@mass@B)
        for k in range(12):dM[:,:,k]+=w*(DB[:,:,k].T@mass@B+B.T@mass@DB[:,:,k])
    return M,dM
def gyro(dM,v):return np.einsum('ijk,j,k->i',dM,v,v)-.5*np.einsum('jki,j,k->i',dM,v,v)
