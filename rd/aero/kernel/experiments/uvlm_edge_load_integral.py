"""Candidate higher-order line-force quadrature; not yet an adopted load policy."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from numpy.polynomial.legendre import leggauss
import uvlm_wake_assembly_v2 as aero
import uvlm_analytic_induction_v2 as induction
def field_at_fraction(model,g,wakes,gamma,p,edge,frac,dp,dgamma,dwakes,dmu):
    kappa=model.D2@gamma;dkappa=model.D2@dgamma;value=np.zeros(3);der=value.copy()
    for i,(a,b) in enumerate(model.edges):
        v,dv=aero.segment_basis(p,g['X'][a],g['X'][b],model.eps,i==edge,dp,g['dX'][a],g['dX'][b],frac if i==edge else None);value+=kappa[i]*v;der+=dkappa[i]*v+kappa[i]*dv
    offset=np.cumsum([0]+[w.mu.size for w in wakes])
    for iw,idx,ids,x,name,r,j,tri in model.wake_faces(wakes):
        v,jac,_=induction.triangle(p,x,model.eps,False);dx=np.array([dwakes[iw][r,j],dwakes[iw][r,j+1],dwakes[iw][r+1,j+1],dwakes[iw][r+1,j]])[tri];dv=jac@np.r_[dp,dx.ravel(),0.];mu=wakes[iw].mu.ravel()[idx];dmm=dmu[offset[iw]+idx];value+=mu@v;der+=dmm@v+mu@dv
    return value,der
def loads(model,g,wakes,gamma,gamma_dot,flow,time,order=32,dgamma=None,dgamma_dot=None,dwakes=None,dmu=None,dU=None):
    if order<1 or not isinstance(order,int):raise ValueError('invalid_line_quadrature_order')
    dgamma=np.zeros_like(gamma) if dgamma is None else dgamma;dgamma_dot=np.zeros_like(gamma_dot) if dgamma_dot is None else dgamma_dot;dwakes=[np.zeros_like(w.vertices) for w in wakes] if dwakes is None else dwakes;dmu=np.zeros(sum(w.mu.size for w in wakes)) if dmu is None else dmu
    out=aero.loads(model,g,wakes,np.zeros_like(gamma),gamma_dot,flow,time,dgamma=np.zeros_like(gamma),dgamma_dot=dgamma_dot,dwakes=dwakes,dmu=dmu,dU=dU);base=out['circulatory_count'];Fv=out['Fv'].copy();dFv=out['dFv'].copy();positions=[];forces=[];dforces=[];kappa=model.D2@gamma;dkappa=model.D2@dgamma;z,w=leggauss(order)
    for edge in model.force_edges:
        a,b=model.edges[edge];length=g['X'][b]-g['X'][a];dl=g['dX'][b]-g['dX'][a]
        for s,weight in zip((z+1)/2,w/2):
            p=(1-s)*g['X'][a]+s*g['X'][b];dp=(1-s)*g['dX'][a]+s*g['dX'][b];ind,dind=field_at_fraction(model,g,wakes,gamma,p,edge,s,dp,dgamma,dwakes,dmu);relative=flow.value(p,time)+ind-(1-s)*g['V'][a]-s*g['V'][b];drel=flow.jvp(p,time,dp,dU)+dind-(1-s)*g['dV'][a]-s*g['dV'][b];f=weight*model.rho*kappa[edge]*np.cross(relative,length);df=weight*model.rho*(dkappa[edge]*np.cross(relative,length)+kappa[edge]*(np.cross(drel,length)+np.cross(relative,dl)));Fv[a]+=(1-s)*f;Fv[b]+=s*f;dFv[a]+=(1-s)*df;dFv[b]+=s*df;positions.append(p);forces.append(f);dforces.append(df)
    count=len(forces);positions=np.vstack([positions,out['positions'][base:]]);forces=np.vstack([forces,out['forces'][base:]]);dforces=np.vstack([dforces,out['dforces'][base:]])
    out.update(Q=model.Nv.T@Fv,dQ=model.Nv.T@dFv,Fv=Fv,dFv=dFv,positions=positions,forces=forces,dforces=dforces,circulatory_count=count,F=forces.sum(axis=0),M=np.cross(positions,forces).sum(axis=0),line_quadrature_order=order)
    if any(not np.isfinite(x).all() for x in out.values() if isinstance(x,np.ndarray)):raise ValueError('nonfinite_integrated_load')
    return out
