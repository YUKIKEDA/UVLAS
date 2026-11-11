"""Independent observation of actual physical-load solutions and preserved projection fields."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
import pg_physical_pde as pde
f=pde.f;m=pde.m

def observe(a,s,mode,prior,order=15):
    reference=pde.validation;reference.select(mode);error=np.zeros(3);norm=np.zeros(3);pe=np.zeros(3);pn=np.zeros(3);field_e=np.zeros(3);field_u=np.zeros(3);internal=np.zeros_like(s['F']);energy=np.zeros(3);ns=a['ns']
    for cell,q in enumerate(a['geo']['refs']):
        p,w=f.c.volume_rule(q,a['n'],order,9)
        for pp,ww in f.chunks(p,w):
            N,G,B,det,_=f.c.kinematics(a['n'],cell,q,a['h'],pp,a['geo']);Eg,Ug=pde.loads.global_fields(q,pp,a['h'],a['k'],m);E,U=reference.f.old.reference(q,pp,a['h'],a['k']);El,Ul=f.bending(q,pp,a['h'],a['k'],a['n']);strain=B@s['coef']+Eg@s['beta'];u=(N@s['coef'].reshape(ns,9)).reshape(len(pp),3,3)+Ug@s['beta'];previous_strain=B@prior['coef']+El@prior['beta'];previous_u=(N@prior['coef'].reshape(ns,9)).reshape(len(pp),3,3)+Ul@prior['beta'];delta=strain-E;stress=m.C@strain;weight=ww*det/a['h']**3
            error+=np.einsum('p,pkt,kl,plt->t',weight,delta,m.C,delta);norm+=np.einsum('p,pkt,kl,plt->t',weight,E,m.C,E);pe+=np.einsum('p,pat,pat->t',ww*det,u-U,u-U);pn+=np.einsum('p,pat,pat->t',ww*det,U,U);field_e+=np.einsum('p,pkt,kl,plt->t',weight,strain-previous_strain,m.C,strain-previous_strain);field_u+=np.einsum('p,pat,pat->t',ww*det,u-previous_u,u-previous_u);internal+=np.einsum('p,pkn,pkt->nt',weight,B,stress,optimize=True);energy+=np.einsum('p,pkt,pkt->t',weight,strain,stress)
    residual=internal-s['F'];g=a['scale'][:,None]*(a['Z'].T@residual);stationarity=np.linalg.norm(g,axis=0)/np.sqrt(a['Z'].shape[1]*norm);L=s['coarse_load'];work=np.sum(s['F']*s['coef'],axis=0)+np.sum(L*s['beta'],axis=0);work_defect=abs(work-energy)/np.maximum(abs(energy),1e-30)
    q0=reference.reference_coefficients(a);nb=(a['n']+2)**2
    def resultant(F):return np.vstack((F.reshape(ns,3,3)[:nb].sum(axis=0),np.cross(q0[:,:,None],F.reshape(ns,3,3),axisa=1,axisb=1,axisc=1).sum(axis=0)))
    balance=np.linalg.norm(resultant(internal),axis=0)/np.maximum(1.,np.linalg.norm(resultant(s['F']),axis=0));reaction_reference=np.linalg.norm(resultant(residual-s['fixed']),axis=0)/np.maximum(1.,np.linalg.norm(resultant(s['fixed']),axis=0))
    boundary=0.
    for cell,q in enumerate(a['geo']['refs']):
        if cell%a['n']:continue
        p=np.array([[-1.,t,z] for t in np.linspace(-1,1,57) for z in [-1.,-.31,0.,.47,1.]]);N=f.shape_values(a,cell,q,p);_,U=pde.loads.global_fields(q,p,a['h'],a['k'],m);value=(N@s['coef'].reshape(ns,9)).reshape(len(p),3,3)+U@s['beta'];boundary=max(boundary,float(np.max(abs(value))))
    return dict(error=np.sqrt(error/norm).tolist(),position=np.sqrt(pe/pn).tolist(),projection_field_energy_difference=np.sqrt(field_e/norm).tolist(),projection_field_position_difference=np.sqrt(field_u/pn).tolist(),stationarity=stationarity.tolist(),elastic_energy=energy.tolist(),external_work=work.tolist(),relative_work_defect=work_defect.tolist(),relative_force_moment_balance=balance.tolist(),relative_reaction_resultant_difference=reaction_reference.tolist(),fixed_boundary=boundary),residual
