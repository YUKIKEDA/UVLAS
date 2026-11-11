"""Independent Ptera mesh/kernel and continuum-time limit of its load definitions."""
import numpy as np
from scipy.linalg import solve
from numpy.polynomial.legendre import leggauss
import pterasoftware as ps
from pterasoftware._aerodynamics_functions import expanded_velocities_from_ring_vortices
def geometry(m,n):
    foil=ps.geometry.airfoil.Airfoil('NACA0012');sections=[ps.geometry.wing_cross_section.WingCrossSection(foil,n//2,chord=1.5,spanwise_spacing='cosine',control_surface_symmetry_type='symmetric'),ps.geometry.wing_cross_section.WingCrossSection(foil,None,chord=.6,Lp_Wcsp_Lpp=[4*np.tan(np.deg2rad(35)),4,0],control_surface_symmetry_type='symmetric')];wing=ps.geometry.wing.Wing(sections,num_chordwise_panels=m,symmetric=True,symmetryNormal_G=[0,1,0],symmetryPoint_G_Cg=[0,0,0]);plane=ps.geometry.airplane.Airplane([wing],s_ref=8.4,c_ref=1.114285714285714,b_ref=8.)
    wing=plane.wings[0];panels=wing.panels;flat=list(panels.ravel());P=np.array([p.Cpp_G_Cg for p in flat]);normal=np.array([p.unitNormal_G for p in flat]);area=np.array([p.area for p in flat]);left=np.array([p.Flbvp_G_Cg for p in flat]);right=np.array([p.Frbvp_G_Cg for p in flat]);X=np.empty((m+1,n+1,3));X[:-1,:-1]=left.reshape(m,n,3);X[:-1,-1]=right.reshape(m,n,3)[:,-1];X[-1,:-1]=np.array([p.Blpp_G_Cg for p in panels[-1]]);X[-1,-1]=panels[-1,-1].Brpp_G_Cg
    return dict(X=X,P=P,normal=normal,area=area,left=left,right=right,m=m,n=n)
def rings(p,fl,fr,br,bl,eps):
    n=len(fl);counts=np.zeros(4,np.int64)
    return expanded_velocities_from_ring_vortices(p,br,fr,fl,bl,-np.ones(n),np.full(n,eps),counts,ages=np.zeros(n),nu=0.)
def project(g,fl,fr,br,bl,eps):
    B=np.empty((len(g['P']),len(fl)))
    for lo in range(0,len(B),64):
        K=rings(g['P'][lo:lo+64],fl,fr,br,bl,eps);B[lo:lo+64]=np.einsum('pij,pj->pi',K,g['normal'][lo:lo+64])
    return B
def run(m,n,k=.2,eps=1e-5,length=256.,order=12,refine_near=False):
    g=geometry(m,n);X=g['X'];omega=2*k*20/1.114285714285714;A=project(g,X[:-1,:-1].reshape(-1,3),X[:-1,1:].reshape(-1,3),X[1:,1:].reshape(-1,3),X[1:,:-1].reshape(-1,3),eps).astype(complex);te=X[-1];tail=te+[length,0,0];W=project(g,te[:-1],te[1:],tail[1:],tail[:-1],eps).astype(complex)
    T=length/20;nodes=[0.];a=0.;step=.0005
    if refine_near:
        e=te[1:]-te[:-1];r=g['P'][:,None]-te[:-1][None];s=np.clip(np.einsum('pij,ij->pi',r,e)/np.sum(e*e,axis=1)[None],0,1);distance=np.sqrt(np.min(np.sum((r-s[:,:,None]*e[None])**2,axis=-1)));step=min(step,distance/40)
    while a<T:a=min(T,a+min(step,np.pi/max(omega,1e-10)));nodes.append(a);step*=1.4
    gx,gw=leggauss(order);lo=np.array(nodes[:-1]);hi=np.array(nodes[1:]);ages=((lo+hi)[:,None]/2+(hi-lo)[:,None]*gx/2).ravel();weights=((hi-lo)[:,None]*gw/2).ravel()*(-1j*omega)*np.exp(-1j*omega*ages)
    for i in range(0,len(ages),8):
        tau=ages[i:i+8];front=te[None]+np.column_stack([20*tau,np.zeros((len(tau),2))])[:,None,:];back=np.broadcast_to(tail,front.shape);B=project(g,front[:,:-1].reshape(-1,3),front[:,1:].reshape(-1,3),back[:,1:].reshape(-1,3),back[:,:-1].reshape(-1,3),eps).reshape(len(g['P']),len(tau),n);W+=np.einsum('q,pqn->pn',weights[i:i+8],B)
    A[:,-n:]+=W;rhs=1j*omega*g['normal'][:,2];gamma=solve(A,rhs);G=gamma.reshape(m,n);horse=(G-np.vstack([np.zeros((1,n)),G[:-1]])).ravel();edgeF=1.2*horse[:,None]*np.cross([20.,0,0],g['right']-g['left']);unsteady=1.2j*omega*gamma[:,None]*g['area'][:,None]*g['normal'];F=edgeF.sum(axis=0)+unsteady.sum(axis=0);M=np.cross((g['left']+g['right'])/2-[.375,0,0],edgeF).sum(axis=0)+np.cross(g['P']-[.375,0,0],unsteady).sum(axis=0)
    return dict(g=g,gamma=gamma,F=F,M=M,W=W,residual=float(np.linalg.norm(A@gamma-rhs)/np.linalg.norm(rhs)))
