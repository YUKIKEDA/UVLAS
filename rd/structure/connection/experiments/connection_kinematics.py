"""Nonlinear attachment frame, rigid-body mass, and consistent derivatives."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
import beam_jet2 as ad
import beam_se3_v2 as group
from scipy.linalg import qr,svd
def cross(a,b):return np.array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],dtype=object)
def normalize(a):return a/ad.sqrt(sum(x*x for x in a))
def frame(y):
    p=y[:3];a=y[3:6];b=y[6:9];e1=normalize(a);e3=normalize(cross(a,b));e2=cross(e3,e1);return p,np.stack((e1,e2,e3),axis=1)
def validate(y):
    a=y[3:6];b=y[6:9]
    if np.shape(y)!=(9,) or not np.isfinite(y).all() or min(np.linalg.norm(a),np.linalg.norm(b))<=0 or np.linalg.norm(np.cross(a,b))<=1e-12*np.linalg.norm(a)*np.linalg.norm(b):raise ValueError('degenerate_attachment_frame')
def prepared(y,c):
    validate(y);a=np.array([ad.Jet.variable(x,i) for i,x in enumerate(y)],object);p,R=frame(a);center=p+R@np.asarray(c);Q=group.float_array(R)
    J=np.array([x.g[:9] for x in center]);H=np.array([x.H[:9,:9] for x in center]);JR=np.array([[x.g[:9] for x in row] for row in R]);HR=np.array([[x.H[:9,:9] for x in row] for row in R]);W=np.empty((3,9));DW=np.empty((3,9,9))
    for j in range(9):
        W[:,j]=np.array(group.vee(Q.T@JR[:,:,j]),float)
        for k in range(9):DW[:,j,k]=np.array(group.vee(JR[:,:,k].T@JR[:,:,j]+Q.T@HR[:,:,j,k]),float)
    return dict(center=group.float_array(center),R=Q,J=J,H=H,W=W,DW=DW)
def rigid_mass(y,m,c,Jc):
    group.spatial_mass(m,c,Jc);v=prepared(y,c);J=v['J'];H=v['H'];W=v['W'];DW=v['DW'];M=m*J.T@J+W.T@Jc@W;DM=np.empty((9,9,9))
    for k in range(9):DM[:,:,k]=m*(H[:,:,k].T@J+J.T@H[:,:,k])+DW[:,:,k].T@Jc@W+W.T@Jc@DW[:,:,k]
    return M,DM,v
def shoe(q0,translation,p0,p,R):
    eta=q0-translation[:,None]*p0;ns=len(q0);G=np.zeros((3*ns,6));H=np.zeros((3*ns,6,6));value=eta@R.T+translation[:,None]*p
    for i in range(ns):
        G[3*i:3*i+3,:3]=translation[i]*np.eye(3);G[3*i:3*i+3,3:]=-R@np.array(group.skew(eta[i]),float)
        for a in range(3):
            A=np.array(group.skew(np.eye(3)[a]),float)
            for b in range(3):
                B=np.array(group.skew(np.eye(3)[b]),float);H[3*i:3*i+3,a+3,b+3]=R@((A@B+B@A)/2)@eta[i]
    return value,G,H
def rigid_weld_graph(nodes,edges,tolerance=1e-9):
    # The global redundancy proof is group holonomy, not rank at one state.
    adjacency={n:[] for n in nodes};kept=[];redundant=[]
    for ident,a,b,H in edges:
        if a not in nodes or b not in nodes or np.shape(H)!=(4,4) or not np.isfinite(H).all() or np.linalg.norm(H[3]-[0,0,0,1])>1e-12:raise ValueError('invalid_weld_edge')
        Q=H[:3,:3]
        if np.linalg.norm(Q.T@Q-np.eye(3))>1e-10 or np.linalg.det(Q)<=0:raise ValueError('invalid_weld_rotation')
        if np.linalg.norm(nodes[a]@H-nodes[b])>tolerance:raise ValueError('reference_weld_mismatch')
        stack=[(a,np.eye(4))];visited=set();found=None
        while stack:
            u,T=stack.pop()
            if u==b:found=T;break
            if u in visited:continue
            visited.add(u)
            stack.extend((v,T@S) for v,S in adjacency[u] if v not in visited)
        if found is not None:
            if np.linalg.norm(found-H)>tolerance:raise ValueError('inconsistent_weld_cycle')
            redundant.append(ident)
        else:kept.append(ident);adjacency[a].append((b,H));adjacency[b].append((a,np.linalg.inv(H)))
    return kept,redundant
class MassLedger:
    def __init__(self):self.ids=set();self.sources=set();self.records=[]
    def add(self,ident,source,mode='add',host_already_included=False,removed_region=None):
        if not ident or ident in self.ids or not source or source in self.sources:raise ValueError('duplicate_mass_or_source')
        if mode not in ['add','replace'] or mode=='add' and host_already_included or mode=='replace' and not removed_region:raise ValueError('mass_ownership_or_replacement')
        self.ids.add(ident);self.sources.add(source);self.records.append((ident,source,mode,removed_region))
def reduce_linear_constraints(A,rhs_coefficients):
    A=np.asarray(A,float);rhs=np.asarray(rhs_coefficients,float)
    if A.ndim!=2 or rhs.ndim!=2 or len(A)!=len(rhs) or not np.isfinite(A).all() or not np.isfinite(rhs).all():raise ValueError('constraint_shape_or_finite')
    norms=np.linalg.norm(A,axis=1)
    if np.any(norms==0):raise ValueError('zero_constraint_row')
    B=A/norms[:,None];s=svd(B,compute_uv=False);ranks=[int(sum(s>tol*s[0])) for tol in [1e-9,1e-10,1e-11]]
    if len(set(ranks))!=1:raise ValueError('unstable_constraint_rank')
    _,_,order=qr(B.T,pivoting=True,mode='economic');keep=order[:ranks[0]];T=np.linalg.lstsq(A[keep].T,A.T,rcond=None)[0].T
    if np.linalg.norm(T@A[keep]-A)>1e-10*np.linalg.norm(A) or np.linalg.norm(T@rhs[keep]-rhs)>1e-10*max(1.,np.linalg.norm(rhs)):raise ValueError('inconsistent_dependent_constraint')
    return keep,T,ranks
