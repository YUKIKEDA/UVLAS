"""One actual adopted beam-shell system with all load paths and tangents."""
from pathlib import Path
import numpy as np
from scipy.linalg import solve,block_diag,null_space
import time,json,hashlib,traceback,sys
import check_connection_beam_shell as coupled
import structural_load_damping as law
import beam_se3_v2 as beam
import connection_shell as pg
from check_structural_load_damping import face_data,chart,fd,rel
ROOT=Path(__file__).resolve().parents[2]

def prepare():
    d=coupled.prepare();shell=d['d'];m=shell['model'];N,Da,Db,w,area=face_data(m,0,1.,21);npnt,_,_=pg.point(m,np.array([[1.,.8,1.]]));marker=np.kron(npnt,np.eye(3));ks=np.diag([.4,.3,.2]);cs=np.diag([.03,.02,.01]);g=np.array([0.,-.02,-.01]);shellgravity=(shell['M']@shell['translation'])[:,None]*g;nodes,weights=np.polynomial.legendre.leggauss(6);d.update(face=(N,Da,Db,w,area),marker=marker,Ks=ks,Cs=cs,g=g,shellgravity=shellgravity.ravel(),beam_quadrature=list(zip((nodes+1)/2,weights/2)),dead=(d['face'][:,None]*np.array([1.,.001,-.001])).ravel(),mass=.3,offset=np.array([0.,.03,.02]),Ic=np.diag([.0003,.0004,.0005]));return d

def evaluate(d,p,R,u,lam,load=1.):
    U,fi,Ki=coupled.internal(d,p,R,u);C,J,H=coupled.constraints(d,p,R,u);rs=np.array([[0.,0.,0.],p]);Rs=np.array([np.eye(3),R]);beamforce=np.zeros(12);beamK=np.zeros((12,12));wrenchF=np.zeros(3);wrenchM=np.zeros(3)
    for t,w in d['beam_quadrature']:
        f,k,info=law.point_wrench(rs,Rs,t,np.zeros(3),.1*w*d['g'],np.zeros(3));beamforce+=f;beamK+=k;wrenchF+=info['force_world'];wrenchM+=np.cross(info['point'],info['force_world'])
    f,k,info=law.point_wrench(rs,Rs,1.,d['offset'],d['mass']*d['g'],np.zeros(3));beamforce+=f;beamK+=k;wrenchF+=info['force_world'];wrenchM+=np.cross(info['point'],info['force_world']);N,Da,Db,w,area=d['face'];pressure,Kp,pressure_info=law.surface(N,Da,Db,d['d']['q0']+u,w,'pressure',.02,area);exshell=d['dead']+d['shellgravity']+pressure.ravel();ext=np.r_[beamforce[6:],exshell]*load;Kext=block_diag(beamK[6:,6:],Kp)*load;e=d['marker']@u.ravel();sf=d['marker'].T@d['Ks']@e;supportK=d['marker'].T@d['Ks']@d['marker'];fi[6:]+=sf;Ki[6:,6:]+=supportK;KL=Ki-Kext;KL[3:6,3:6]+=np.einsum('i,iab->ab',lam,H);force=fi-ext+J.T@lam;A=np.block([[KL,J.T],[J,np.zeros((len(C),len(C)))]])
    # Independent physical wrench ledger, including force at each actual point.
    qs=d['d']['q0']+u;fs=d['dead']+d['shellgravity'];wrenchF+=d['d']['translation']@fs.reshape(-1,3);wrenchM+=np.sum(np.cross(qs,fs.reshape(-1,3)),axis=0);wrenchF+=np.sum(w[:,None]*pressure_info['density'],axis=0);wrenchM+=np.sum(w[:,None]*np.cross(pressure_info['point'],pressure_info['density']),axis=0);supportforce=-d['Ks']@e;wrenchF=load*wrenchF+supportforce;wrenchM=load*wrenchM+np.cross((d['marker']@qs.ravel()),supportforce);kb,_=coupled.section();root=beam.evaluate(rs,Rs,1.,kb,np.array([1.,0,0,0,0,0]))['f'][:6]-load*beamforce[:6];balance=np.r_[wrenchF+root[:3],wrenchM+root[3:]]
    return np.r_[force,C],A,dict(KL=KL,J=J,external=ext,root_reaction=root,wrench_balance=balance,pressure_K=Kp,support_force=sf)

def main():
    start=time.monotonic();rows=[];failures=[];artifact=None
    try:
        d=prepare();p=d['p0'].copy();R=np.eye(3);u=np.zeros_like(d['d']['q0']);lam=np.zeros(3*len(d['T']));history=[]
        for load in [.25,.5,1.]:
            for iteration in range(20):
                rr,A,info=evaluate(d,p,R,u,lam,load);n=6+u.size;error=float(np.linalg.norm(rr[:n])/np.linalg.norm(info['external']));gap=float(np.max(abs(rr[n:])));history.append(dict(load=load,iteration=iteration,original_residual=error,constraint=gap));print('COMBINED',history[-1],flush=True)
                if error<1e-9 and gap<1e-10:break
                sq=1/np.sqrt(np.maximum(abs(np.diag(A)[:n]),1e-10));sl=1/np.maximum(np.linalg.norm(info['J']*sq[None,:],axis=1),1e-12);s=np.r_[sq,sl];dx=s*solve(A*s[:,None]*s[None,:],-s*rr,assume_a='gen');p,R=coupled.update(p,R,dx[:6]);u+=dx[6:n].reshape(-1,3);lam+=dx[n:]
            else:raise ValueError('combined_equilibrium_not_converged')
        balance=rel(info['wrench_balance'],info['external']);assert balance<1e-8;rows.append(dict(case='combined_actual_nonlinear_equilibrium',history=history,global_wrench_error=balance,elastic_support_retained=True,eccentric_mass_gravity=True,current_pressure=True))
        # Test all original residual blocks in ONE fixed exponential chart.
        rng=np.random.default_rng(1505);direction=rng.normal(size=len(rr))*.001;n=6+u.size;gradient=np.einsum('ia,pib->pab',direction[6:n].reshape(-1,3),d['d']['G']);direction[6:n]*=min(1.,.1/np.max(abs(gradient)));off=direction*.1;po,Ro=coupled.update(p,R,off[:6]);uo=u+off[6:n].reshape(-1,3);lo=lam+off[n:];base,Ak,_=evaluate(d,po,Ro,uo,lo)
        def residual(z):
            pp,RR=coupled.update(po,Ro,z[:6]);uu=uo+z[6:n].reshape(-1,3);ll=lo+z[n:];rr2=evaluate(d,pp,RR,uu,ll)[0];T=block_diag(np.eye(3),beam.float_array(beam.left_jacobian(-z[3:6])),np.eye(u.size));return np.r_[T.T@rr2[:n],rr2[n:]]
        numeric=fd(residual,direction,h=.002);tangent_error=rel(numeric-Ak@direction,numeric);assert tangent_error<1e-8;rows.append(dict(case='whole_original_residual_directional_tangent',relative_error=tangent_error,off_equilibrium_test=True,nonsymmetric_load_tangent_retained=True))
        Z=null_space(info['J']);probe=rng.normal(size=n);sol=solve(A,np.r_[probe,np.zeros(len(lam))]);independent=Z@solve(Z.T@info['KL']@Z,Z.T@probe);difference=rel(sol[:n]-independent,independent);raw=rel(A@sol-np.r_[probe,np.zeros(len(lam))],probe);assert max(difference,raw)<1e-8;rows.append(dict(case='coupled_linearized_KKT_vs_independent_elimination',state_difference=difference,original_residual=raw))
        factor=law.viscosity_factor(np.diag([.4,.7,.9,.2,.3,.5]));_,_,Cv,_=law.shell_damping(d['d'],u,np.zeros_like(u),factor);rs=np.array([[0.,0.,0.],p]);Rs=np.array([np.eye(3),R]);_,_,Cb,_=law.beam_damping(rs,Rs,1.,np.zeros(12),factor);kb,mb=coupled.section();Mb,_=beam.inertia(rs,Rs,1.,mb,8);x,Q=beam.point_jet(rs,Rs,1.,d['offset']);Q0=beam.float_array(Q);Jx=np.array([a.g for a in x]);JQ=np.array([[a.g for a in row] for row in Q]);W=np.column_stack([np.array(beam.vee(Q0.T@JQ[:,:,j]),float) for j in range(12)]);Mbody=d['mass']*Jx.T@Jx+W.T@d['Ic']@W;M=block_diag((Mb+Mbody)[6:,6:],np.kron(d['d']['M'],np.eye(3)));C=block_diag(Cb[6:,6:],Cv+d['marker'].T@d['Cs']@d['marker']);massmin=float(np.linalg.eigvalsh(Z.T@M@Z).min());dampmin=float(np.linalg.eigvalsh(Z.T@C@Z).min());assert massmin>0 and dampmin>-1e-10*np.linalg.norm(C);rows.append(dict(case='same_equilibrium_mass_stiffness_damping',reduced_mass_min=massmin,reduced_damping_min=dampmin,stiffness_asymmetry=rel(info['KL']-info['KL'].T,info['KL'])))
        if '--publish' in sys.argv:
            out=ROOT/'build/analysis/structural-combined-snapshot.npz';assert not out.exists();np.savez_compressed(out,p=p,R=R,u=u,lam=lam,M=M,K=info['KL'],C=C,J=info['J'],original_residual=rr,root_reaction=info['root_reaction']);artifact=dict(artifact=out.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    print(json.dumps(rows,indent=2),flush=True)
    if '--publish' in sys.argv:
        paths=[Path(__file__),Path(coupled.__file__),Path(law.__file__),ROOT/'tools/design-experiments/check_structural_load_damping.py',ROOT/'docs/verification/structural-load-damping-result.json',ROOT/'docs/verification/connection-adoption-audit.json',ROOT/'docs/verification/structural-final-integration-protocol.md'];dest=ROOT/'docs/verification/structural-combined-result.json';assert not dest.exists();dest.write_text(json.dumps(dict(passed=len(rows)==4 and not failures,records=rows,failures=failures,artifact=artifact,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print('STRUCTURAL_COMBINED_RESULT',len(rows)==4 and not failures,flush=True)
if __name__=='__main__':main()
