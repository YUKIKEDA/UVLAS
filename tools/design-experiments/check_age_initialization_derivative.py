from pathlib import Path
import json,hashlib,numpy as np
from scipy.linalg import solve
from numpy.polynomial.legendre import leggauss
from check_coupled_bdf2_prescribed import operators,nodal_wake
from uvlm_age_scalar import age_integral
ROOT=Path(__file__).resolve().parents[2]
def main():
    g,A,C,H,Bq,Bv=operators();nb=len(A);n=g['n'];M=np.diag([40.,8.]);omega=10.;q=np.array([1.,.2],complex);v=1j*omega*q;acc=-omega**2*q;gx,gw=leggauss(12);errors=[];records=[]
    for R in [20,40,80,160]:
        h=1/R;W=nodal_wake(g['P'],g['normal'],g['X'][-1],1e-5,h,R,gx,gw);AC=A.astype(complex);AC[:,-n:]+=np.einsum('prn,r->pn',W,np.exp(-1j*omega*h*np.arange(R+1)));gamma=solve(AC,Bq@q+Bv@v);mu=np.exp(-1j*omega*h*np.arange(R+1))[:,None]*gamma[-n:];mdot=np.empty_like(mu);mdot[0]=0;mdot[1]=-(mu[2]-mu[0])/(2*h);mdot[2:]=-(3*mu[2:]-4*mu[1:-1]+mu[:-2])/(2*h);AE=A.copy();AE[:,-n:]+=W[:,0];block=np.block([[M,-H],[-Bv,AE]]);right=np.r_[M@acc-H@(1j*omega*gamma),Bq@v-np.einsum('prn,rn->p',W[:,1:],mdot[1:])];answer=solve(block,right);err=[np.linalg.norm(answer[:2]-acc)/np.linalg.norm(acc),np.linalg.norm(answer[2:]-1j*omega*gamma)/np.linalg.norm(omega*gamma)];errors.append(err);records.append(dict(intervals=R,relative_acceleration_and_GammaDot_errors=[float(x) for x in err]))
    orders=np.log2(np.array(errors[-2])/np.array(errors[-1]));print("AGE_ORDERS",orders.tolist(),"ERRORS",errors,flush=True);assert min(orders)>1.8 and max(orders)<2.2 and max(errors[-1])<.01
    steady=A.copy();steady[:,-n:]+=W.sum(axis=1);gam=solve(steady,-.1*g['normal'][:,2]);Ktail=age_integral(g['P'],g['normal'],g['X'][-1],1e-5,np.array([1.]),np.array([1.+0j]));naive=solve(block,np.r_[np.zeros(2),-Ktail@gam[-n:]]);assert np.linalg.norm(naive)>1e-8;paths=[Path(__file__),ROOT/'docs/j1-initialization-output-derivative-addendum.md',ROOT/'docs/verification/coupled-bdf2-prescribed-result.json'];result=dict(passed=True,scope='fixed-age initial derivative with actual aero operators and added-mass block; material-tail negative control',records=records,final_orders=orders.tolist(),naive_material_tail_spurious_initial_derivative=float(np.linalg.norm(naive)),sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'docs/verification/age-initialization-derivative-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
