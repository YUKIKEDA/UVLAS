"""QR-compressed full constitutive strain approximation on a developable surface."""
from functools import lru_cache
import itertools
import numpy as np
from scipy.linalg import qr,svd,solve_triangular
import q9_curved_energy_scaling as material
import check_q9_curved_best_approximation_bounds as rounding
core=material.core;CH=np.linalg.cholesky(material.CM).T

@lru_cache(None)
def table(order):
    points,w=core.quadrature(order);values=[core.shapes.shape(x,y) for x,y,_ in points]
    return points,w,np.array([v[0] for v in values]),np.array([np.column_stack(v[1:]) for v in values])
@lru_cache(None)
def surface_table(order):
    p,w=np.polynomial.legendre.leggauss(order);points=list(itertools.product(p,p));values=[core.shapes.shape(*point) for point in points]
    return np.array([v[0] for v in values]),np.array([np.column_stack(v[1:]) for v in values]),np.array([a*b for a,b in itertools.product(w,w)])
def mesh(k,n,amplitude):
    _,xy,cells=core.geometry.grid.mesh(n,amplitude);xy=xy.copy();xy[:,1]-=.5;coords=np.zeros((len(xy),3,3));x,y=xy.T;coords[:,0]=np.column_stack((x,y,k*x*x/2));coords[:,1]=np.column_stack((-k*x,np.zeros(len(x)),np.ones(len(x))));u=np.zeros((len(xy),9));u[:,:3]=np.column_stack((-k*x**3/3,np.zeros(len(x)),x*x/2));u[:,3:6]=np.column_stack((-x,np.zeros(len(x)),-k*x*x));side=2*n+1;free=np.array([i for i in range(len(xy)) if i%side!=0]);return xy,cells,coords,u,free
def target(k,x,z):
    a=np.column_stack((np.ones_like(x),np.zeros_like(x),k*x));j=1+k*k*x*x-k*z;tensor=-z[:,None,None]*a[:,:,None]*a[:,None,:]/j[:,None,None];return np.einsum('aij,pij->pa',core.method.MB,tensor)
def compress(A):
    Q,R=qr(A,mode='economic',check_finite=False);product=Q@R;residual=float(np.linalg.norm(A-product));del product
    delta=residual+rounding.gamma(A.shape[1]+3)*float(np.linalg.norm(abs(Q)@abs(R)))/(1-rounding.gamma(A.shape[1]+3));orth=float(np.linalg.norm(Q.T@Q-np.eye(Q.shape[1])));eta=orth+rounding.gamma(A.shape[0]+3)*float(np.linalg.norm(abs(Q).T@abs(Q)))/(1-rounding.gamma(A.shape[0]+3));return R,dict(residual=residual,error_upper=delta,orthogonality=orth,orthogonality_upper=eta)
def geometry_error(k,h,xy,coords,J,points,N,g):
    base=N@xy;x,y=base.T;z=points[:,2]*h/2;expected=np.column_stack((x-k*x*z,y,k*x*x/2+z));phi,_=core.geometry.raw.tables(points,h,0.);actual=phi@coords.reshape(27,3);baseJ=np.einsum('ia,pib->pab',xy,g);analytic=np.zeros((len(points),3,3));analytic[:,0,0]=1-k*z;analytic[:,1,1]=1;analytic[:,2,0]=k*x;analytic[:,0,2]=-k*x;analytic[:,2,2]=1;chart=np.zeros_like(analytic);chart[:,:2,:2]=baseJ;chart[:,2,2]=h/2;expectedJ=analytic@chart;error=float(np.max(np.linalg.norm(J-expectedJ,axis=1)/np.maximum(np.linalg.norm(expectedJ,axis=1),1e-300)));return float(np.max(abs(actual-expected))),error
def assemble(k,n,h,amplitude,states):
    xy,cells,coords,u,free=mesh(k,n,amplitude);lookup={int(v):i for i,v in enumerate(free)};nr=len(free)*3;np_=len(free)*6;rows=82*len(cells);routes={name:dict(W=np.zeros((rows,nr)),D=np.zeros((rows,np_)),e=np.zeros(rows)) for name in ['raw','cg']};points,w,N,g=table(9);NS,gs,ws=surface_table(5);M=np.zeros((len(xy),len(xy)));area=0.;reference_norm=0.;interpolation={name:0. for name in routes};change=0.;positions=0.;jacobian=0.;edge=True;edge_condition=0.;locals_=[]
    for cell,(_,_,ids,_) in enumerate(cells):
        q=coords[ids].reshape(81)
        if cell not in states:states[cell]=core.method.prepare(q,h,0.,9)
        state=states[cell];J,E,B=core.method.evaluate(state,points);weight=w*np.linalg.det(J)
        if np.min(weight)<=0:raise ValueError('nonpositive_reference_measure')
        pe,je=geometry_error(k,h,xy[ids],coords[ids],J,points,N,g);positions=max(positions,pe);jacobian=max(jacobian,je);edge=edge and all(len(set(e['ranks']))==1 and e['condition']<=1e8 for e in state['edges'].values());edge_condition=max(edge_condition,max(e['condition'] for e in state['edges'].values()));x=(N@xy[ids])[:,0];star=target(k,x,h*points[:,2]/2);weighted_star=(star@CH.T*np.sqrt(weight)[:,None]).reshape(-1);reference_norm+=weighted_star@weighted_star;sl=slice(82*cell,82*(cell+1));local_matrices={}
        for name,field in [('raw',E[:,:,:81]),('cg',B[:,:,:81])]:
            weighted=(np.einsum('ij,pjk->pik',CH,field)*np.sqrt(weight)[:,None,None]).reshape(-1,81);aug=np.column_stack((weighted,weighted_star));R,qrlog=compress(aug);local_matrices[name]=weighted;scale=1/np.linalg.norm(local_matrices['raw'],axis=0);s=svd(R[:,:81]*scale[None,:],compute_uv=False);ranks=[int(np.count_nonzero(s>t*s[0])) for t in [1e-9,1e-10,1e-11]];modes=[]
            for axis in np.eye(3):
                v=np.zeros((9,3,3));v[:,0,:]=axis;modes.extend([v.reshape(81),np.cross(np.broadcast_to(axis,coords[ids].shape),coords[ids]).reshape(81)])
            rigid=max(float(np.linalg.norm((R[:,:81]*scale[None,:])@(v/scale))/(s[0]*np.linalg.norm(v/scale))) for v in modes);locals_.append(dict(cell=cell,route=name,ranks=ranks,rigid=rigid,qr=qrlog));routes[name]['e'][sl]=R[:,81]
            for local,node in enumerate(ids):
                if int(node) not in lookup:continue
                j=lookup[int(node)];routes[name]['W'][sl,3*j:3*j+3]=R[:,9*local:9*local+3];routes[name]['D'][sl,6*j:6*j+6]=R[:,9*local+3:9*local+9]/np.repeat([1.,core.LENGTH],3)[None,:]
            residual=weighted@u[ids].reshape(81)-weighted_star;interpolation[name]+=residual@residual
        diff=(local_matrices['cg']-local_matrices['raw'])@u[ids].reshape(81);change+=diff@diff
        surfaceJ=np.einsum('ia,pib->pab',coords[ids,0,:],gs);metric=surfaceJ.transpose(0,2,1)@surfaceJ;aw=ws*np.sqrt(np.linalg.det(metric));area+=aw.sum();localH=np.einsum('p,pi,pj->ij',aw,NS,NS)+core.LENGTH**2*np.einsum('p,pia,pab,pjb->ij',aw,gs,np.linalg.inv(metric),gs);M[np.ix_(ids,ids)]+=localH
    M=M[np.ix_(free,free)]/area;pref=u[free,3:9].copy();pref[:,3:]*=core.LENGTH;return dict(routes=routes,M=M,xy=xy,coords=coords,free=free,cells=cells,reference=u,reference_p=pref,reference_norm=float(reference_norm),interpolation={name:float(np.sqrt(value/reference_norm)) for name,value in interpolation.items()},reduction_change=float(np.sqrt(change/reference_norm)),geometry_position=positions,geometry_jacobian=jacobian,edge_valid=bool(edge),edge_condition=edge_condition,locals=locals_)
def solve_route(a,name,factors):
    data=a['routes'][name];normalization=np.sqrt(a['reference_norm']);aug=np.column_stack((data['W'],data['D'],data['e']))/normalization;R,log=compress(aug);nr=data['W'].shape[1];W=R[:nr,:nr];D=R[nr:,nr:-1];e=R[nr:,-1];singular_r=svd(W,compute_uv=False)
    if singular_r[-1]<=1e-11*singular_r[0]:raise ValueError('position_elimination_rank')
    A,chol=core.whiten(D,a['M']);U,s,Vh=svd(A,full_matrices=False,check_finite=False);reference=core.pnorm(a['reference_p'],a['M']);solutions=[]
    for factor in factors:
        fit=core.bounded(U,s,Vh,-e,factor*reference,A,chol);p=fit.pop('p');r=solve_triangular(W,R[:nr,-1]-R[:nr,nr:-1]@p,lower=False,check_finite=False);residual=R[:,:-1]@np.r_[r,p]-R[:,-1];gradient_r=R[:,:nr].T@residual;Hp=(a['M']@p.reshape(len(a['M']),6)).reshape(-1);gradient_p=R[:,nr:-1].T@residual+fit['multiplier']*Hp;denr=np.linalg.norm(R[:,:nr])*np.linalg.norm(R[:,-1]);denp=np.linalg.norm(R[:,nr:-1])*np.linalg.norm(R[:,-1])+fit['multiplier']*np.linalg.norm(Hp);fit.update(factor=factor,reference_h1=reference,full_error=float(np.linalg.norm(residual)),kkt_r=float(np.linalg.norm(gradient_r)/max(denr,1e-300)),kkt_p=float(np.linalg.norm(gradient_p)/max(denp,1e-300)),r=r,p=p)
        solutions.append(fit)
    return R,log,solutions,dict(position_condition=float(singular_r[0]/singular_r[-1]),direction_singular_values=s)
def full_field(a,r,p):
    u=np.zeros((len(a['xy']),9));u[a['free'],:3]=r.reshape(-1,3);u[a['free'],3:6]=p.reshape(-1,6)[:,:3];u[a['free'],6:9]=p.reshape(-1,6)[:,3:]/core.LENGTH;return u
def observe(k,h,a,states,solutions,order=11):
    points,w,N,g=table(order);total=0.;out=[dict(error=0.,energy=0.,error_lower_square=0.,error_upper_square=0.) for _ in solutions]
    for cell,(_,_,ids,_) in enumerate(a['cells']):
        J,E,B=core.method.evaluate(states[cell],points);weight=w*np.linalg.det(J);x=(N@a['xy'][ids])[:,0];star=target(k,x,h*points[:,2]/2);weighted=star@CH.T;total+=float(np.einsum('p,pi,pi->',weight,weighted,weighted))
        matrices={name:(np.einsum('ij,pjk->pik',CH,field[:,:,:81])*np.sqrt(weight)[:,None,None]).reshape(-1,81) for name,field in [('raw',E),('cg',B)]};star_weighted=(weighted*np.sqrt(weight)[:,None]).reshape(-1)
        for dest,(name,u) in zip(out,solutions):
            value,ve=rounding.dot_rows(matrices[name],u[ids].reshape(81));error=value-star_weighted;ee=ve+rounding.gamma(3)*(abs(value)+abs(star_weighted));low,high=rounding.norms(error,ee);dest['error_lower_square']+=low*low;dest['error_upper_square']+=high*high;dest['error']+=float(error@error);dest['energy']+=float(value@value)
    margin=rounding.gamma(len(points)*6*len(a['cells'])+10)
    for dest in out:
        dest['error']=float(np.sqrt(dest['error']/total));dest['energy_ratio']=dest['energy']/total;dest['error_lower']=float(np.sqrt(max(0.,dest['error_lower_square']*(1-margin)/(total*(1+margin)))));dest['error_upper']=float(np.sqrt(dest['error_upper_square']*(1+margin)/(total*(1-margin))))
    return total,out
def position_error(k,a,u,order):
    N,g,w=surface_table(order);error=0.;total=0.
    for _,_,ids,_ in a['cells']:
        xy=N@a['xy'][ids];x=xy[:,0];star=np.column_stack((-k*x**3/3,np.zeros(len(x)),x*x/2));Js=np.einsum('ia,pib->pab',a['coords'][ids,0,:],g);weight=w*np.sqrt(np.linalg.det(Js.transpose(0,2,1)@Js));value=N@u[ids,:3];error+=float(np.einsum('p,pi,pi->',weight,value-star,value-star));total+=float(np.einsum('p,pi,pi->',weight,star,star))
    return float(np.sqrt(max(0.,error/total)))
