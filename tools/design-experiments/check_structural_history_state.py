from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.linalg import solve
import json,hashlib,time,traceback,copy,sys
import structural_history_state as state
import structural_load_damping as law
import check_structural_combined as combined
import check_connection_beam_shell as coupled
ROOT=Path(__file__).resolve().parents[2]
def main():
    start=time.monotonic();rows=[];failures=[]
    try:
        knots=np.array([0.,.3,.8,1.4]);poly=lambda t:np.array([t**3-2*t*t+.3*t+.1,.4*t**3+.2*t]);p1=lambda t:np.array([3*t*t-4*t+.3,1.2*t*t+.2]);p2=lambda t:np.array([6*t-4,2.4*t]);h=state.History(knots,np.array([poly(t) for t in knots]),np.array([p1(t) for t in knots]),role='prescribed_motion');errors=[0.,0.,0.]
        for t in np.linspace(.001,1.399,31):
            for j,fn in enumerate([poly,p1,p2]):errors[j]=max(errors[j],float(np.max(abs(h.evaluate(t,j)-fn(t)))))
        assert max(errors)<1e-10;linear=state.History([0.,.3,1.],[0.,.6,.2]);assert abs(linear.evaluate(.15)-.3)<1e-15 and linear.events==(.3,);assert abs(linear.evaluate(.3,1,'left')-2)<1e-14 and abs(linear.evaluate(.3,1,'right')+4/7)<1e-14;rows.append(dict(case='explicit_linear_and_C1_Hermite_histories',polynomial_errors=errors,all_internal_knots_registered_as_events=True,no_implicit_extrapolation=True))
        rejected=0
        for callback in [lambda:state.History([0,0],[1,2]),lambda:state.History([0,1],[1,float('nan')]),lambda:state.History([0,1],[0,1],role='prescribed_motion'),lambda:linear.evaluate(-.1),lambda:linear.evaluate(.3,1),lambda:h.evaluate(.3,2)]:
            try:callback()
            except ValueError:rejected+=1
            else:raise AssertionError('invalid history accepted')
        assert rejected==6;rows.append(dict(case='history_negative_controls',rejected=rejected))
        # Independent steady harmonic oscillator and dissipated cycle power.
        mass=.3;k=55.;c=.2;omega=2.3;U=1/(k-mass*omega**2+1j*c*omega);K=np.diag([k,0,0,0,0,0]);D=np.diag([c,0,0,0,0,0]);J=np.zeros((6,1));J[0]=1.;H=np.zeros((6,1,1))
        def rhs(t,y):
            e=np.r_[y[0],np.zeros(5)];ed=np.r_[y[1],np.zeros(5)];f=law.strain_support(e,J,H,ed,K,D,np.zeros((6,1)))[0][0];return [y[1],(np.cos(omega*t)-f)/mass]
        tt=np.linspace(0,2*np.pi/omega,101);sol=solve_ivp(rhs,[tt[0],tt[-1]],[U.real,(1j*omega*U).real],t_eval=tt,method='DOP853',rtol=1e-11,atol=1e-13);exact=np.vstack(((U*np.exp(1j*omega*tt)).real,(1j*omega*U*np.exp(1j*omega*tt)).real));error=float(np.max(abs(sol.y-exact))/np.max(abs(exact)));inputpower=.5*np.real(np.conj(1j*omega*U));loss=.5*c*omega**2*abs(U)**2;power=abs(inputpower-loss)/loss;assert sol.success and max(error,power)<1e-8;rows.append(dict(case='finite_elastic_viscous_support_vs_harmonic_ODE',response_error=error,cycle_power_error=power,not_product_time_integrator_adoption=True))
        d=coupled.prepare();p=d['p0'].copy();R=np.eye(3);u=np.zeros_like(d['d']['q0']);lam=np.zeros(3*len(d['T']));ns=u.size;n=6+ns;spring=5.;probe=np.zeros(n);probe[6::3]=d['face'];error=1.
        for iteration in range(15):
            _,f,K0=coupled.internal(d,p,R,u);C,J,H=coupled.constraints(d,p,R,u);tip=float(probe[6:]@u.ravel());res=f-probe+spring*tip*probe+J.T@lam;KL=K0+spring*np.outer(probe,probe);KL[3:6,3:6]+=np.einsum('i,iab->ab',lam,H);A=np.block([[KL,J.T],[J,np.zeros((len(C),len(C)))]])
            error=float(np.linalg.norm(res)/np.linalg.norm(probe))
            if error<1e-9 and np.max(abs(C))<1e-10:break
            dx=solve(A,-np.r_[res,C],assume_a='sym');p,R=coupled.update(p,R,dx[:6]);u+=dx[6:n].reshape(-1,3);lam+=dx[n:]
        else:raise AssertionError('supported series equilibrium failed')
        def series(T):
            stretch=brentq(lambda x:50*x*(x*x-1)-T,1,1.1);return T/50+stretch-1
        Tref=brentq(lambda T:T+spring*series(T)-1,0,1);expected=series(Tref);tiperror=abs(tip/expected-1);assert tiperror<1e-8;rows.append(dict(case='actual_supported_beam_shell_vs_independent_series',tip_error=tiperror,original_residual=error,finite_spring=5.))
        # Save/restore the actual combined case, not a synthetic state vector.
        snap=ROOT/'build/analysis/structural-combined-snapshot.npz';a=np.load(snap);model=combined.prepare();initial={k:a[k].copy() for k in ['p','R','u','lam']};identity=hashlib.sha256((ROOT/'docs/verification/structural-combined-result.json').read_bytes()).hexdigest();store=state.Store(initial,identity)
        def evaluate(s):return combined.evaluate(model,s['p'],s['R'],s['u'],s['lam'])
        before=store.fingerprint();trial=store.begin();trial['state']['R'][0,0]=2.
        try:store.evaluate(trial,evaluate)
        except ValueError:pass
        else:raise AssertionError('bad SO3 trial accepted')
        assert store.fingerprint()==before
        try:store.accept(trial)
        except ValueError:pass
        else:raise AssertionError('unqualified state committed')
        good=store.begin();original=store.evaluate(good,evaluate);mutated=copy.deepcopy(good);mutated['state']['p'][0]+=.1
        try:store.accept(mutated)
        except ValueError:pass
        else:raise AssertionError('post-evaluation mutation accepted')
        assert store.fingerprint()==before;store.accept(good);payload=store.checkpoint();restored=state.Store.restore(payload,identity,evaluate);output=evaluate(restored.begin()['state']);difference=max(float(np.max(abs(original[i]-output[i]))) for i in [0,1]);assert difference==0 and store.fingerprint()==restored.fingerprint();before=store.fingerprint()
        try:state.Store.restore(payload,'wrong-model',evaluate)
        except ValueError:pass
        else:raise AssertionError('checkpoint generation mismatch accepted')
        try:store.accept(good)
        except ValueError:pass
        else:raise AssertionError('stale trial accepted')
        assert store.fingerprint()==before;rows.append(dict(case='actual_state_atomic_failure_and_restart',operator_difference=difference,bitwise_identity=True,negative_controls=5,scope='binary64-qualified structural physical state; coupled integrator/wake format is work 3/4'))
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    print(json.dumps(rows,indent=2),flush=True)
    if '--publish' in sys.argv:
        dest=ROOT/'docs/verification/structural-history-state-result.json';assert not dest.exists();paths=[Path(__file__),Path(state.__file__),Path(law.__file__),Path(combined.__file__),ROOT/'docs/verification/structural-combined-result.json',ROOT/'docs/verification/structural-final-integration-protocol.md'];dest.write_text(json.dumps(dict(passed=len(rows)==5 and not failures,records=rows,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print('STRUCTURAL_HISTORY_STATE_RESULT',len(rows)==5 and not failures,flush=True)
if __name__=='__main__':main()
