"""Ball-support lower bound after QR position elimination, with QR drift allowance."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
import check_q9_curved_best_approximation_bounds as safe

def vector_residual(A,p,e):
    v,err=safe.dot_rows(A,p);out=v-e;err=err+safe.gamma(3)*(abs(v)+abs(e));return out,err
def dual_norm(g,error,M,hmin):
    z=np.linalg.solve(M,g.reshape(len(M),6));q=np.zeros_like(z);qe=np.zeros_like(z)
    for j in range(6):
        mz,me=safe.dot_rows(M,z[:,j]);q[:,j]=g.reshape(z.shape)[:,j]-mz;qe[:,j]=error.reshape(z.shape)[:,j]+me+safe.gamma(3)*(abs(g.reshape(z.shape)[:,j])+abs(mz))
    v,e=safe.metric_square(M,z.reshape(-1));root=np.nextafter(np.sqrt(hmin),0.);return float(np.nextafter((np.sqrt(max(0.,v+e))+safe.norms(q.reshape(-1),qe.reshape(-1))[1]/root)*(1+safe.gamma(12)),np.inf))
def certify(R,nr,M,r,p,B,global_qr,local_delta,local_eta,original_zero_upper):
    hmin=safe.positive_lower(M);D=R[nr:,nr:-1];e=R[nr:,-1];res,err=vector_residual(D,p,e);lo,hi=safe.norms(res,err);gp,gpe=safe.dot_rows(D.T,res);prop,pe=safe.dot_rows(abs(D.T),err);gpe=gpe+prop+pe;dn=dual_norm(gp,gpe,M,hmin);dot,de=safe.dot(gp,p);dotup=dot+de+safe.upper_sum(abs(p)*gpe);lb=lo*lo/2-dotup-B*dn;lb-=safe.gamma(32)*(lo*lo/2+abs(dotup)+B*dn);lb=float(np.nextafter(lb,-np.inf))
    rr=r.copy();pp=p.copy();scale=1.
    for _ in range(6):
        value,error=safe.metric_square(M,pp);pn=np.sqrt(max(0.,value+error))
        if pn<=B:break
        factor=(B/pn)*(1-1e-10);rr*=factor;pp*=factor;scale*=factor
    full,fe=vector_residual(R[:,:-1],np.r_[rr,pp],R[:,-1]);upper=safe.norms(full,fe)[1];target=safe.norms(R[:,-1]);target_lower=target[0];target_upper=target[1];qr_lower=np.sqrt(max(0.,2*lb));W=R[:nr,:nr];G=W.T@W;g=safe.gamma(nr+3);gram_error=g*(abs(W).T@abs(W))/(1-g);sigma_square=safe.positive_lower(G)-np.linalg.norm(gram_error)/(1-safe.gamma(nr*nr+3));sigma=float(np.sqrt(max(0.,sigma_square)))
    loc=max(0.,1-local_eta);glob=max(0.,1-global_qr['orthogonality_upper']);factor=np.sqrt(loc*glob);drift=np.sqrt(loc)*global_qr['error_upper']+local_delta;pnorm=B/np.sqrt(hmin);Dnorm=float(np.linalg.norm(R[:,nr:-1]))/(1-safe.gamma(R.shape[0]*(R.shape[1]-nr)+3));denominator=factor*sigma-drift
    if denominator<=0:raise ValueError('unbounded_qr_position_level_set')
    rnorm=(original_zero_upper+factor*(Dnorm*pnorm+target_upper)+drift*(pnorm+1))/denominator;tbound=np.sqrt(rnorm*rnorm+pnorm*pnorm+1);original_lower=max(0.,factor*qr_lower-drift*tbound);original_lower=max(0.,original_lower-safe.gamma(32)*(abs(factor*qr_lower)+abs(drift*tbound)))
    return dict(lower_objective=lb,qr_lower_relative=float(qr_lower/target_upper),qr_upper_relative=float(upper/target_lower),original9_lower_relative=float(original_lower/original_zero_upper),qr_error_allowance=float(drift*tbound),position_level_bound=float(rnorm),position_singular_lower=sigma,h1_gram_lower=hmin,scale=float(scale),feasible=bool(pn<=B),objective_order=bool(qr_lower<=upper),normalized_target_norm=target_upper,feasible_r=rr,feasible_p=pp)

def synthetic():
    R=np.array([[1.,0.,0.,0.,0.,0.,0.,0.],[0.,1.,0.,0.,0.,0.,0.,2.]])
    # One free position, six bounded direction coordinates; the reduced target is 2.
    record=certify(R,1,np.ones((1,1)),np.array([0.]),np.array([.5,0.,0.,0.,0.,0.]),.5,dict(error_upper=0.,orthogonality_upper=0.),0.,0.,2.)
    exact=1.5/2
    assert record['qr_lower_relative']<=exact<=record['qr_upper_relative'] and record['objective_order'] and record['feasible']
    return record
