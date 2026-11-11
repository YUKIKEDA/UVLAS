"""T-E24-CG1: constrain bubble-gradient modes to zero mean and trace."""
import numpy as np
import q9_edge_trace_moment as base
from q9_edge_moment_fields import quadrature
MB=base.MB
def modes(state,points,E):
    points=np.asarray(points);bubble=(1-points[:,0]**2)*(1-points[:,1]**2)/state['bubble_mean']
    psi=E[:,:,72:75]
    return psi-bubble[:,None,None]*state['mode_mean'][None,:,:]
def prepare(q,width,offset,order):
    state=base.prepare(q,width,offset,order);points,w=quadrature(order);J,E,R=base.evaluate(state,points);measure=w*np.linalg.det(J);V=measure.sum()
    state['mode_mean']=np.einsum('p,pki->ki',measure,E[:,:,72:75])/V
    phi=modes(state,points,E);G=np.einsum('p,pki,pkj->ij',measure,phi,phi,optimize=True);scale=1/np.sqrt(G.diagonal());gs=scale[:,None]*G*scale[None,:];sv=np.linalg.svd(gs,compute_uv=False)
    ranks=[int(np.count_nonzero(sv>t*sv[0])) for t in [1e-9,1e-10,1e-11]]
    if ranks!=[3,3,3] or sv[0]/sv[-1]>1e8:raise ValueError('invalid_constrained_mode_gram')
    rhs=np.einsum('p,pki,pkn->in',measure,phi,E-R,optimize=True)
    state['mode_alpha']=np.linalg.solve(G,rhs);state['mode_gram']=G;state['mode_gram_ranks']=ranks;state['mode_gram_condition']=float(sv[0]/sv[-1])
    return state
def evaluate(state,points):
    J,E,R=base.evaluate(state,points)
    return J,E,R+np.einsum('pki,in->pkn',modes(state,points,E),state['mode_alpha'],optimize=True)
