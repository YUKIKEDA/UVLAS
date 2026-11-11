"""T-E24-G1: flat extrusion only, restore the exact central bubble strain."""
import numpy as np
import q9_edge_trace_moment as base
from q9_edge_moment_fields import quadrature
MB=base.MB
def prepare(q,width,offset,order):
    if offset!=0:raise ValueError('symmetric_thickness_required')
    coords=q.reshape(9,3,3);J,_=base.raw.kinematics(q,width,offset,[[0.,0.,0.]])
    normal=np.cross(J[0,:,0],J[0,:,1]);normal/=np.linalg.norm(normal)
    if np.max(abs((coords[:,0]-coords[0,0])@normal))>1e-12 or np.max(abs(coords[:,1]-coords[0,1]))>1e-12 or np.linalg.norm(np.cross(coords[0,1],normal))>1e-12 or np.max(abs(coords[:,2]))>1e-12:raise ValueError('flat_extrusion_required')
    state=base.prepare(q,width,offset,order);u=np.zeros(81);u[72:75]=normal;points,w=quadrature(order);J,E,R=base.evaluate(state,points);measure=w*np.linalg.det(J);psi=E[:,:,:81]@u;den=np.einsum('p,pi,pi->',measure,psi,psi)
    if den<=0:raise ValueError('nonpositive_bubble_strain_norm')
    state['gradient_mode']=u;state['gradient_alpha']=np.einsum('p,pi,pin->n',measure,psi,E-R)/den
    return state
def evaluate(state,points):
    J,E,R=base.evaluate(state,points);psi=E[:,:,:81]@state['gradient_mode']
    return J,E,R+psi[:,:,None]*state['gradient_alpha'][None,None,:]
