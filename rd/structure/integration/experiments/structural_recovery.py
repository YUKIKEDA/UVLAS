"""Typed physical stresses and cut wrenches, without equilibrium correction."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from pg_ritz_mechanics import MB
import beam_se3_v2 as beam

def stress(F,C,e0,viscous_second_PK=None):
    F=np.asarray(F);det=np.linalg.det(F)
    if not np.isfinite(F).all() or np.any(det<=0):raise ValueError('invalid_recovery_deformation')
    E=(np.swapaxes(F,-1,-2)@F-np.eye(3))/2;e=np.einsum('kab,...ab->...k',MB,E)-e0;s=np.einsum('...kl,...l->...k',C,e);S=np.einsum('...k,kab->...ab',s,MB)
    if viscous_second_PK is not None:S=S+viscous_second_PK
    P=F@S;sigma=P@np.swapaxes(F,-1,-2)/det[...,None,None]
    return dict(green=E,second_PK=S,first_PK=P,cauchy=sigma,volume_J=det)

def cut_wrench(X,F,P,reference_area,weights,origin):
    """Reference cut area is oriented. Current area uses Nanson's map.

    Returns the actual stress integral; NEVER force-fit to boundary reactions.
    """
    det=np.linalg.det(F)
    if np.any(det<=0):raise ValueError('inverted_cut')
    area=np.einsum('pab,pb->pa',det[:,None,None]*np.swapaxes(np.linalg.inv(F),1,2),reference_area);traction=np.einsum('pab,pb->pa',P,reference_area);force=np.sum(weights[:,None]*traction,axis=0);moment=np.sum(weights[:,None]*np.cross(X-origin,traction),axis=0)
    return dict(force=force,moment=moment,current_area=area,reference_traction_density=traction)

def beam_cut(r,R,L,K,reference,e0,t,origin):
    state=beam.evaluate(r,R,L,K,reference,e0);x,Q=beam.pose_at(r[0],R[0],state['d'],t);x=beam.float_array(x);Q=beam.float_array(Q);mu=K@state['strain'];force=Q@mu[:3];moment=Q@mu[3:]+np.cross(x-origin,force);return dict(point=x,force=force,moment=moment,section_conjugate=mu)
