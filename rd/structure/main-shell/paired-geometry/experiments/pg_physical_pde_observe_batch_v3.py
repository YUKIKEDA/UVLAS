"""Independent observation of actual physical-load solutions and preserved projection fields."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
import pg_physical_pde_batch_v3 as batch
pde=batch.pde
f=pde.f;m=pde.m
import curvature_work_profiles as curvature

def preserved_fields(q,p,h,k,n,beta):
    N,_=m.shape(p,0);x=(N@q.reshape(75,3)[:25])[:,0];z=h*p[:,2]/2;K,I=curvature.evaluate(x,n,0,[k,0.]);K=K@beta;I0=I[0]@beta;I1=I[1]@beta;I2=I[2]@beta;b=x[:,None]*I0-I1;C=(x[:,None]**2*I0-I2)/2;a=np.column_stack((np.ones(len(x)),np.zeros(len(x)),k*x));tensor=-z[:,None,None,None]*a[:,:,None,None]*a[:,None,:,None]*K[:,None,None,:]/(1+k*k*x*x-k*z)[:,None,None,None];E=np.einsum('kij,pijt->pkt',m.MB,tensor);U=np.zeros((len(p),3,beta.shape[1]));U[:,0]=-k*C-z[:,None]*I0;U[:,2]=b-k*x[:,None]*z[:,None]*I0;return E,U


def observe(a,s,mode,prior,order=15):
    reference=pde.validation;nt=6;error=np.zeros(nt);norm=np.zeros(nt);pe=np.zeros(nt);pn=np.zeros(nt);field_e=np.zeros(nt);field_u=np.zeros(nt);internal=np.zeros_like(s['F']);energy=np.zeros(nt);ns=a['ns']
    for cell,q in enumerate(a['geo']['refs']):
        p,w=f.c.volume_rule(q,a['n'],order,9)
        for pp,ww in f.chunks(p,w):
            N,G,B,det,_=f.c.kinematics(a['n'],cell,q,a['h'],pp,a['geo']);Eg,Ug=pde.loads.global_fields(q,pp,a['h'],a['k'],m);Er,Ur=reference.rational_reference(q,pp,a['h'],a['k']);Es,Us=reference.STANDARD_REFERENCE(q,pp,a['h'],a['k']);E=np.concatenate((Er,Es),axis=2);U=np.concatenate((Ur,Us),axis=2);El,Ul=preserved_fields(q,pp,a['h'],a['k'],a['n'],prior['beta']);strain=B@s['coef']+Eg@s['beta'];u=(N@s['coef'].reshape(ns,3*nt)).reshape(len(pp),3,nt)+Ug@s['beta'];previous_strain=B@prior['coef']+El;previous_u=(N@prior['coef'].reshape(ns,3*nt)).reshape(len(pp),3,nt)+Ul;delta=strain-E;stress=m.C@strain;weight=ww*det/a['h']**3
            error+=np.einsum('p,pkt,kl,plt->t',weight,delta,m.C,delta);norm+=np.einsum('p,pkt,kl,plt->t',weight,E,m.C,E);pe+=np.einsum('p,pat,pat->t',ww*det,u-U,u-U);pn+=np.einsum('p,pat,pat->t',ww*det,U,U);field_e+=np.einsum('p,pkt,kl,plt->t',weight,strain-previous_strain,m.C,strain-previous_strain);field_u+=np.einsum('p,pat,pat->t',ww*det,u-previous_u,u-previous_u);internal+=np.einsum('p,pkn,pkt->nt',weight,B,stress,optimize=True);energy+=np.einsum('p,pkt,pkt->t',weight,strain,stress)
    residual=internal-s['F'];g=a['scale'][:,None]*(a['Z'].T@residual);stationarity=np.linalg.norm(g,axis=0)/np.sqrt(a['Z'].shape[1]*norm);L=s['coarse_load'];work=np.sum(s['F']*s['coef'],axis=0)+np.sum(L*s['beta'],axis=0);work_defect=abs(work-energy)/np.maximum(abs(energy),1e-30)
    q0=reference.reference_coefficients(a);nb=(a['n']+2)**2
    def resultant(F):return np.vstack((F.reshape(ns,3,nt)[:nb].sum(axis=0),np.cross(q0[:,:,None],F.reshape(ns,3,nt),axisa=1,axisb=1,axisc=1).sum(axis=0)))
    balance=np.linalg.norm(resultant(internal),axis=0)/np.maximum(1.,np.linalg.norm(resultant(s['F']),axis=0));reaction_reference=np.linalg.norm(resultant(residual-s['fixed']),axis=0)/np.maximum(1.,np.linalg.norm(resultant(s['fixed']),axis=0))
    boundary=0.
    for cell,q in enumerate(a['geo']['refs']):
        if cell%a['n']:continue
        p=np.array([[-1.,t,z] for t in np.linspace(-1,1,57) for z in [-1.,-.31,0.,.47,1.]]);N=f.shape_values(a,cell,q,p);_,U=pde.loads.global_fields(q,p,a['h'],a['k'],m);value=(N@s['coef'].reshape(ns,3*nt)).reshape(len(p),3,nt)+U@s['beta'];boundary=max(boundary,float(np.max(abs(value))))
    return dict(error=np.sqrt(error/norm).tolist(),position=np.sqrt(pe/pn).tolist(),projection_field_energy_difference=np.sqrt(field_e/norm).tolist(),projection_field_position_difference=np.sqrt(field_u/pn).tolist(),stationarity=stationarity.tolist(),elastic_energy=energy.tolist(),external_work=work.tolist(),relative_work_defect=work_defect.tolist(),relative_force_moment_balance=balance.tolist(),relative_reaction_resultant_difference=reaction_reference.tolist(),fixed_boundary=boundary),residual
