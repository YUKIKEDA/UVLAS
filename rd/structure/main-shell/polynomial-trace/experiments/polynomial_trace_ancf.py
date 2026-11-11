"""ANCF hybrid reconstruction with polynomial-consistent trace stabilization."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import itertools
import numpy as np
from scipy.linalg import svd,qr,solve_triangular,solve
import hybrid_gradient_ancf as old
m=old.m;qrtools=old.qrtools
FACES=old.FACES;mesh=old.mesh;edge_nodes=old.edge_nodes;reference_state=old.reference_state;volume_points=old.volume_points;ranks=old.ranks
POT_INDICES=list(itertools.product(range(6),range(6),range(3)))
LOW=np.array([i for i,(a,b,c) in enumerate(POT_INDICES) if a<=4 and b<=4]);HIGH=np.array([i for i in range(108) if i not in LOW])

def potential_basis(q,h,p,independent=False):
    value=[];derivative=[]
    for axis,degree in [(0,5),(1,5),(2,2)]:
        value.append(np.polynomial.legendre.legvander(p[:,axis],degree))
        d=[]
        for j in range(degree+1):
            c=np.eye(degree+1)[j]
            if independent:
                pc=np.polynomial.legendre.leg2poly(c);d.append(np.polynomial.polynomial.polyval(p[:,axis],np.polynomial.polynomial.polyder(pc)))
            else:d.append(np.polynomial.legendre.legval(p[:,axis],np.polynomial.legendre.legder(c)))
        derivative.append(np.array(d).T)
    phi=np.einsum('pi,pj,pk->pijk',*value).reshape(len(p),108)
    grad=np.stack([np.einsum('pi,pj,pk->pijk',*(derivative[j] if j==axis else value[j] for j in range(3))).reshape(len(p),108) for axis in range(3)],axis=2)
    grad[:,:,2]*=2/h;_,J,_=m.raw(q,p,h*p[:,2]/2,independent)
    return phi,(grad@np.linalg.inv(J)).transpose(0,2,1)

def face_legendre(points,axis):
    return np.einsum('pi,pj->pij',np.polynomial.legendre.legvander(points[:,1-axis],4),np.polynomial.legendre.legvander(points[:,2],2)).reshape(len(points),15)

def prepare(q,h,kind,order=13,zorder=7,independent=False):
    if kind not in ['weak','potential']:raise ValueError('unknown_reconstruction_variant')
    state=old.prepare(q,h,'S4',order,zorder,independent);phi,Gp=potential_basis(q,h,state['p'],independent);Gp=Gp[:,:,1:]
    W=(Gp*np.sqrt(state['mu'])[:,None,None]).reshape(-1,107);norm=np.linalg.norm(W,axis=0);Wn=W/norm
    rhs=(state['S']*np.sqrt(state['mu'])[:,None,None]).reshape(-1,state['S'].shape[2])
    if independent:
        singular=svd(Wn,compute_uv=False,lapack_driver='gesvd',check_finite=False);Q,R,pivot=qr(Wn,mode='economic',pivoting=True,check_finite=False)
        H=np.zeros((107,rhs.shape[1]));H[pivot]=solve_triangular(R,Q.T@rhs,check_finite=False);H/=norm[:,None]
    else:
        U,singular,V=svd(Wn,full_matrices=False,check_finite=False);H=((V.T/singular)@(U.T@rhs))/norm[:,None]
    rr=ranks(singular)
    if rr!=[107]*3:raise ValueError('potential_gradient_rank_failed:'+str(rr))
    low=phi[:,LOW];high=phi[:,HIGH];M=low.T@(state['mu'][:,None]*low)
    projection=solve(M,low.T@(state['mu'][:,None]*high),assume_a='pos',check_finite=False)
    Hhi=H[HIGH-1];state['H']=H;state['Hhi']=Hhi;state['potential_map']=H@state['D'];state['low_projection']=projection
    state['Gw_weak']=state['Gw'];state['S_body']=state['S'] if kind=='weak' else Gp@H
    state['Gw']=state['S_body']@state['D'];state['kind']=kind;state['potential_ranks']=rr
    for face,(axis,sign) in zip(state['faces'],FACES):
        phif,_=potential_basis(q,h,face['points'],independent);highf=phif[:,HIGH]-phif[:,LOW]@projection
        L=face_legendre(face['points'],axis);mu=face['measure']/face['measure'].sum();Mf=L.T@(mu[:,None]*L)
        coeff=solve(Mf,L.T@(mu[:,None]*highf),assume_a='pos',check_finite=False);E=L@coeff
        face['jump_raw']=face['jump'];face['lift_coeff']=coeff;face['E']=E;face['jump']=face['jump_raw']+E@Hhi@state['D']
    Gw=state['Gw'];B=np.einsum('kij,pnj->pkni',m.MB,Gw.transpose(0,2,1)).reshape(len(state['p']),6,405)
    elastic=(np.einsum('ij,pjn->pin',m.CH,B)*np.sqrt(state['weight'])[:,None,None]).reshape(-1,405)
    stabilization=np.vstack([np.kron(face['jump']*np.sqrt(state['tau']*face['measure'])[:,None],np.eye(3)) for face in state['faces']])
    state['B']=B;state['root']=np.vstack((elastic,stabilization));state['elastic_rows']=len(elastic)
    F0=np.einsum('na,pjn->paj',state['X0'],Gw);state['reference_gradient_error']=float(np.linalg.norm(F0-np.eye(3))/np.sqrt(3*len(F0)))
    return state

def nonlinear(state,current,v1,v2):
    u=current.reshape(135,3)-state['X0'];Gw=state['Gw'];w=state['weight'];F=np.eye(3)+np.einsum('na,pjn->paj',u,Gw)
    E=np.einsum('kij,pij->pk',m.MB,(F.transpose(0,2,1)@F-np.eye(3))/2);stress=E@m.C.T;S=np.einsum('kij,pk->pij',m.MB,stress);P=F@S
    B=np.einsum('kij,pai,pjn->pkna',m.MB,F,Gw).reshape(-1,6,405)
    d1=np.einsum('na,pjn->paj',v1.reshape(135,3),Gw);d2=np.einsum('na,pjn->paj',v2.reshape(135,3),Gw)
    second=np.einsum('kij,pij->pk',m.MB,(d1.transpose(0,2,1)@d2+d2.transpose(0,2,1)@d1)/2)
    U=np.einsum('p,pk,pk->',w,E,stress)/2;force=np.einsum('p,pjn,paj->na',w,Gw,P)
    tangent=np.einsum('p,pk,kl,pl->',w,B@v1,m.C,B@v2)+np.einsum('p,pk,pk->',w,stress,second)
    pcoef=np.einsum('p,pjr,paj->ar',state['mu'],state['S_body'],P,optimize=True);eta=np.zeros((3,len(HIGH)),dtype=F.dtype);jumps=[]
    for face in state['faces']:
        jump=face['jump']@u;jumps.append(jump);weighted=state['tau']*face['measure'];U+=np.einsum('p,pa,pa->',weighted,jump,jump)/2
        force+=np.einsum('p,pn,pa->na',weighted,face['jump'],jump)
        tangent+=np.einsum('p,pa,pa->',weighted,face['jump']@v1.reshape(135,3),face['jump']@v2.reshape(135,3))
        eta+=np.einsum('p,pr,pa->ar',weighted,face['E'],jump)
    dual=pcoef+eta@state['Hhi']/state['volume'];fluxes=[]
    for face,jump in zip(state['faces'],jumps):
        Pface=np.einsum('pjr,ar->paj',face['S'],dual);fluxes.append(np.einsum('paj,pj->pa',Pface,face['normal'])-state['tau']*jump)
    return U,force.reshape(-1),tangent,dict(min_F=float(np.min(np.linalg.det(F.real))),fluxes=fluxes,stabilization_dual_norm=float(np.linalg.norm(eta@state['Hhi']/state['volume'])))

def boundary_map(state,order=None,zorder=None,independent=False):return old.boundary_map(state,order,zorder,independent)

def evaluate_gradient(state,points,independent=True):
    if state['kind']=='weak':A,_,_=old.generators(state['q'],state['h'],points,'S4',independent);S=A@state['coeff']
    else:_,Gp=potential_basis(state['q'],state['h'],points,independent);S=Gp[:,:,1:]@state['H']
    return S@state['D']

def nodal_volume(coeff,h):
    xy=np.array([[s,t] for t,s in itertools.product(m.LOBATTO,repeat=2)]);Ls=np.polynomial.legendre.legvander(xy[:,0],4);Lt=np.polynomial.legendre.legvander(xy[:,1],4)
    z=[np.array([np.polynomial.legendre.legval(0.,np.polynomial.legendre.legder(np.eye(3)[j],a)) for j in range(3)])*(2/h)**a for a in range(3)]
    C=coeff.reshape(5,5,3,-1)
    return np.vstack([np.einsum('pi,pj,k,ijkn->pn',Ls,Lt,a,C) for a in z])

def nodal_face(coeff,h):
    L=np.polynomial.legendre.legvander(m.LOBATTO,4);C=coeff.reshape(5,3,-1)
    z=[np.array([np.polynomial.legendre.legval(0.,np.polynomial.legendre.legder(np.eye(3)[j],a)) for j in range(3)])*(2/h)**a for a in range(3)]
    return np.vstack([np.einsum('pi,k,ikn->pn',L,a,C) for a in z])

def polynomial_patch(state):
    phi,Gp=potential_basis(state['q'],state['h'],state['p']);ids=[i for i,(a,b,c) in enumerate(POT_INDICES) if a+b+c<=5]
    low=phi[:,LOW];test=phi[:,ids];mu=state['mu'];coeff=solve(low.T@(mu[:,None]*low),low.T@(mu[:,None]*test),assume_a='pos')
    u=np.zeros((135,len(ids)));u[:75]=nodal_volume(coeff,state['h'])
    for face,(axis,sign) in zip(state['faces'],FACES):
        phif,_=potential_basis(state['q'],state['h'],face['points']);L=face_legendre(face['points'],axis);mf=face['measure']/face['measure'].sum()
        coeff=solve(L.T@(mf[:,None]*L),L.T@(mf[:,None]*phif[:,ids]),assume_a='pos');u[face['offset']:face['offset']+15]=nodal_face(coeff,state['h'])
    truth=Gp[:,:,ids];observed=state['Gw']@u
    gradient=max(float(np.linalg.norm(observed[:,:,j]-truth[:,:,j])/max(1.,np.linalg.norm(truth[:,:,j]))) for j in range(len(ids)))
    jump=max(float(np.linalg.norm(face['jump']@u[:,j])/max(1.,np.linalg.norm(u[:25,j]))) for face in state['faces'] for j in range(len(ids)))
    return dict(polynomials=len(ids),gradient_error=gradient,stabilization_residual=jump,passed=bool(max(gradient,jump)<=1e-9))
