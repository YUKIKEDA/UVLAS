"""Same exact PG space, with the global smooth curvature subspace protected.

This is a coordinate change, not an additional field or a dropped mode.
The coarse r image is formed directly from exact global polynomial fields;
its director graph is obtained before mixing it into sensitive fine modes.
"""
import numpy as np
from scipy.linalg import svd
import exact_paired_ancf as parent
import paired_coarse_polynomials as coarse_polynomials
polynomial=parent.polynomial
geometry=parent.geometry
lift=parent.lift
canonical=parent.canonical
rank_check=parent.rank_check
evaluate_legendre=parent.evaluate_legendre
compatible_particular=parent.compatible_particular
fields=parent.fields
kinematics=parent.kinematics

def prepare(n,parameters,order=15):
    C,record=canonical(n,parameters);nr=record['exact_r_rank'];nd=record['exact_d_kernel']
    a,b=np.meshgrid(np.arange(14),np.arange(14),indexing='ij')
    weight=np.broadcast_to(1/(2*n)/np.sqrt((2*a+1)*(2*b+1)),(2*n,2*n,14,14)).reshape(-1,1)
    R=C[0,:,:,:,:,:nr];DI=C[1,:,:,:,:,:nr];DK=C[1,:,:,:,:,nr:]
    RA=R.reshape(-1,nr)*weight;DA=DK.reshape(-1,nd)*weight
    U,s,V=svd(RA,full_matrices=False,check_finite=False);rr=rank_check(RA,'scaled_r_image_rank')
    if rr!=[nr]*3:raise ValueError('r_image_dimension')
    Tr=V.T/s;R=R@Tr;DI=DI@Tr
    U,s,V=svd(DA,full_matrices=False,check_finite=False);dr=rank_check(DA,'scaled_d_kernel_rank')
    if dr!=[nd]*3:raise ValueError('d_kernel_dimension')
    Td=V.T/s;DK=DK@Td
    DKvalues=DK.reshape(-1,nd)*weight;DIvalues=DI.reshape(-1,nr)*weight
    cross=DKvalues.T@DIvalues;DI-=DK@cross
    raw_coarse=coarse_polynomials.build(n,parameters)
    RC=raw_coarse[0];DC=raw_coarse[1];RCvalues=RC.reshape(-1,24)*weight
    U,sc,Vc=svd(RCvalues,full_matrices=False,check_finite=False)
    cranks=rank_check(RCvalues,'protected_coarse_r_rank');nc=cranks[0]
    Tc=Vc[:nc].T/sc[:nc];coarse_R=RC@Tc;coarse_D=DC@Tc
    coarse_D-=DK@(DKvalues.T@(coarse_D.reshape(-1,nc)*weight))
    overlap=(coarse_R.reshape(-1,nc)*weight).T@(R.reshape(-1,nr)*weight)
    RF=R-coarse_R@overlap;DF=DI-coarse_D@overlap
    RFvalues=RF.reshape(-1,nr)*weight
    U,sf,Vf=svd(RFvalues,full_matrices=False,check_finite=False)
    franks=rank_check(RFvalues,'protected_fine_r_rank');nf=franks[0]
    if nf!=nr-nc:raise ValueError('protected_r_decomposition')
    Tf=Vf[:nf].T/sf[:nf]
    R=np.concatenate((coarse_R,RF@Tf),axis=4);DI=np.concatenate((coarse_D,DF@Tf),axis=4)
    raw_coordinates=np.zeros((nr+nd,24));raw_coordinates[:nc]=sc[:nc,None]*Vc[:nc]
    difference=DC-coarse_D@raw_coordinates[:nc]
    raw_coordinates[nr:]=DKvalues.T@(difference.reshape(-1,24)*weight)
    # Convert the completed polynomial basis to shared C1/C0 controls for evaluation.
    hx,ht=polynomial.hermite_nodes(n);X,Y=np.meshgrid(hx,hx,indexing='ij');dx,dy=np.meshgrid(ht,ht,indexing='ij')
    hp=np.column_stack((X.ravel(),Y.ravel()));jets=evaluate_legendre(hp,n,R,True)
    index=np.where(dx.ravel(),np.where(dy.ravel(),4,1),np.where(dy.ravel(),2,0));controls=np.zeros_like(jets[0])
    for j in [0,1,2,4]:controls[index==j]=jets[j][index==j]
    controls*=((1/(4*n))**(dx+dy)).ravel()[:,None]
    Rc=controls.reshape(len(hx),len(hx),nr);Rc[0]=0
    gx=polynomial.nodes(n,13);X,Y=np.meshgrid(gx,gx,indexing='ij');p=np.column_stack((X.ravel(),Y.ravel()))
    Dc=evaluate_legendre(p,n,np.concatenate((DI,DK),axis=4))[0].reshape(len(gx),len(gx),nr+nd);Dc[0]=0
    p,w=geometry.quadrature(n,order);nativeC=polynomial.native(p,n)[2][0];indices=np.arange(nativeC.shape[1])[n+1:]
    Tc,cr,cc=polynomial.independent(nativeC[:,indices],np.sqrt(w))
    return dict(n=n,parameters=parameters,Rcoeff=Rc,Dcoeff=Dc,degree=13,Cindices=indices,Tc=Tc,
        nr=nr,nd=nd,nc=Tc.shape[1],ns=nr+nd+Tc.shape[1],lift=dict(n=n,parameters=parameters),
        ranks=dict(r=rr,d=dr,c=cr),canonical=record,nr_coarse=nc,coarse_coordinates=np.vstack((raw_coordinates,np.zeros((Tc.shape[1],24)))))

