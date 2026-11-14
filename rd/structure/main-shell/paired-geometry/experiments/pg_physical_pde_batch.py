"""Same six independent physical RHSs share one geometry observation pass."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
import pg_physical_pde_plane as route
pde=route.pde;f=pde.f;m=pde.m;paired=pde.paired;loads=pde.loads
_curvature=loads.curvature

def curvature(x,mode):
    if mode!='combined':return _curvature(x,mode)
    a,b=_curvature(x,'rational');c,d=_curvature(x,'standard');return np.column_stack((a,c)),np.column_stack((b,d))

loads.curvature=curvature
def external(a,mode,order=15,count=14,beta=None):
    """Integrate only body/traction moments against scalar position functions."""
    nt=6 if mode=='combined' else 3;F=np.zeros((a['ns'],3,nt));fixed=np.zeros_like(F);coarse=np.zeros((8,nt));k,h,n=a['k'],a['h'],a['n']
    def point_data(q,p):
        X,J,jac=xy(q,p);x=X[:,0];v,vp=loads.curvature(x,mode);K,Kx,U=loads.global_basis(x,k)
        if beta is not None:v=v-K@beta;vp=vp-Kx@beta
        W,Fx,Fy=loads.moments(x,k,h,v,vp,count);return J,jac,U,W,Fx,Fy
    for cell,q in enumerate(a['geo']['refs']):
        p,w=f.c.plane_rule(q,n,order);J,jac,U,W,_,_=point_data(q,p);S=components(a,cell,q,p);weight=w*jac
        F+=np.einsum('p,lpn,lpat->nat',weight,S,W,optimize=True);coarse+=np.einsum('p,lpan,lpat->nt',weight,U,W,optimize=True)
        row,col=divmod(cell,n)
        for axis,sign,boundary in [(0,-1,col==0),(0,1,col==n-1),(1,-1,row==0),(1,1,row==n-1)]:
            if not boundary:continue
            # Split physical curvature knots on boundary edges, using the established rule.
            p3,ww=validation.face_rule(q,n,axis,sign,order,1);p=p3[:,:2];J,jac,U,_,Fx,Fy=point_data(q,p);S=components(a,cell,q,p);traction=sign*(Fx if axis==0 else Fy);weight=ww/2*J[:,1-axis,1-axis]
            term=np.einsum('p,lpn,lpat->nat',weight,S,traction,optimize=True)
            if axis==0 and sign<0:fixed+=term
            else:F+=term;coarse+=np.einsum('p,lpan,lpat->nt',weight,U,traction,optimize=True)
    return F.reshape(3*a['ns'],nt)/h**3,fixed.reshape(3*a['ns'],nt)/h**3,coarse/h**3


pde.external=external
