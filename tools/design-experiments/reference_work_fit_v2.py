"""Same RW-3 QR solve; evaluate its objective on original physical strains."""
import numpy as np
import reference_work_fit as prior
from reference_work_fit import f,m,KIND_TARGETS,target,assemble,observe,register,edge_keys

def solve(a):
    sol=prior.solve(a);compressed=sol['fit']['error'];elastic=np.zeros(2);jump=np.zeros(2);ref=np.zeros(2)
    trace_gradient=np.zeros_like(sol['trace']);trace_norm2=np.zeros(sol['trace'].shape[:2]);inside_g=[];target_norm2=0.;ncols=0
    for cell,values in zip(a['local'],sol['local']):
        s=f.prepare(cell['q'],a['h']);W=s['root']/a['h']**1.5;norm2=np.sum(W*W,axis=0);residuals=[]
        for j,kind in enumerate(KIND_TARGETS):
            star=target(cell['q'],s['p'],a['h'],a['grid']['kappa'],kind);E=s['B']@values[:,j]-star
            ep=(E@m.CH.T*np.sqrt(s['weight'])[:,None]).reshape(-1)/a['h']**1.5;elastic[j]+=ep@ep
            jr=[]
            for face in s['faces']:
                val=face['jump']@values[:,j].reshape(-1,3);jr.append((val*np.sqrt(s['tau']*face['measure'])[:,None]).reshape(-1)/a['h']**1.5)
            jp=np.concatenate(jr);jump[j]+=jp@jp;residuals.append(np.r_[ep,jp])
            r=(star@m.CH.T*np.sqrt(s['weight'])[:,None]).reshape(-1)/a['h']**1.5;ref[j]+=r@r
        residual=np.column_stack(residuals);g=W.T@residual;ni=cell['ni'];inside_g.append(g[:ni]/np.sqrt(norm2[:ni,None]));ncols+=ni
        trace_gradient[cell['trace_ids']]+=(cell['trace_transform'].T@g[ni:]).reshape(-1,3,2)
        # These benchmark faces use the canonical Q4 basis, so trace_transform=I.
        # Reject a different transform rather than silently using diagonal norms.
        T=cell['trace_transform']
        if np.linalg.norm(T-np.eye(len(T)))>1e-10:raise ValueError('direct_stationarity_requires_canonical_benchmark_faces')
        trace_norm2[cell['trace_ids']]+=norm2[ni:].reshape(-1,3)
    grad=np.vstack(inside_g+[trace_gradient[a['free']].reshape(-1,2)/np.sqrt(trace_norm2[a['free']].reshape(-1,1))]);ncols+=3*len(a['free'])
    direct_kkt=np.linalg.norm(grad,axis=0)/np.maximum(1e-30,np.sqrt(ncols)*np.sqrt(ref))
    sol['fit']['compressed_error']=compressed;sol['fit']['error']=np.sqrt((elastic+jump)/ref).tolist();sol['fit']['direct_kkt']=direct_kkt.tolist()
    sol['fit']['objective_evaluation']='primary_physical_strain_and_stabilization_same_functional';return sol
