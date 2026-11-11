"""Independent work, pressure-volume, chart and dissipation verification."""
from pathlib import Path
import numpy as np
from scipy.linalg import block_diag
import json,hashlib,time,traceback,sys
import structural_load_damping as law
import connection_shell as pg
import connection_kinematics as kin
import beam_se3_v2 as beam
import beam_jet2 as ad
ROOT=Path(__file__).resolve().parents[2]
def rel(a,b):return float(np.linalg.norm(a)/max(np.linalg.norm(b),1e-30))
def fd(fun,v,h=1e-4):return (-fun(2*h*v)+8*fun(h*v)-8*fun(-h*v)+fun(-2*h*v))/(12*h)

def face_data(model,axis,side,order=29):
    p,w,area=pg.face(model,axis,side,order);N,G,_,_,X=pg.fields(model,p);J=pg.shell.geometry.reference(p[:,:2],p[:,2],model['h'],model['parameters'],False)[1];Da=(2*side-1)*np.einsum('pib,pb->pi',G,J[:,:,(axis+1)%3]);Db=np.einsum('pib,pb->pi',G,J[:,:,(axis+2)%3]);return N,Da,Db,w,area

def chart(delta):
    return block_diag(*[block_diag(np.eye(3),beam.float_array(beam.left_jacobian(-delta[6*i+3:6*i+6]))) for i in range(2)])

def support_map(y,t):
    yy=np.array([ad.Jet.variable(a,i) for i,a in enumerate(y)],object);tt=ad.Jet.variable(t,9);p,R=kin.frame(yy);pb=np.array([.02*ad.sin(tt),.03*ad.cos(tt),.01*tt],object);Rb=beam.exp_rotation(np.array([0.,.1*ad.sin(tt),.08*tt],object));e=np.r_[Rb.T@(p-pb),beam.log_rotation(Rb.T@R)];return np.array([a.v for a in e]),np.array([a.g for a in e]),np.array([a.H for a in e])

def main():
    start=time.monotonic();rows=[];failures=[];rng=np.random.default_rng(20261105)
    try:
        model=pg.shell.prepare(1,(.4,.2,.07),.1,False,'extended');model=dict(model,physical_scale=1.,physical_h=.1,span=1.,width=1.);d=pg.data(model,order=29);q0=d['q0'];t=d['translation'];Rot=beam.float_array(beam.exp_rotation([.7,-.3,.2]));A=Rot@np.array([[1.08,.03,0.],[0.,.94,.02],[.01,0.,1.03]]);q=q0@A.T+t[:,None]*np.array([.2,-.1,.3]);v=rng.normal(size=q.shape)*.001;faces=[face_data(model,a,s) for a in range(3) for s in [0.,1.]];N,Da,Db,w,area0=faces[-1]
        for kind,value in [('pressure',2.3),('traction_current',np.array([.7,-.4,.2])),('traction_reference',np.array([.7,-.4,.2]))]:
            f,K,info=law.surface(N,Da,Db,q,w,kind,value,area0);cs=law.surface(N,Da,Db,q+1e-25j*v,w,kind,value,area0,False)[0].imag/1e-25;derivative=rel(cs-(K@v.ravel()).reshape(q.shape),cs);external_power=float(np.sum(w[:,None]*info['density']*(N@v)));work=abs(float(np.sum(f*v))-external_power)/max(abs(external_power),1e-30);F=t@f;M=np.sum(np.cross(q,f),axis=0);Fref=np.sum(w[:,None]*info['density'],axis=0);Mref=np.sum(w[:,None]*np.cross(info['point'],info['density']),axis=0);wrench=max(rel(F-Fref,Fref),rel(M-Mref,Mref));assert max(derivative,work,wrench)<1e-8
            skew=rel(K-K.T,K);symmetrized_defect=rel((K-K.T)@v.ravel()/2,cs.ravel()) if kind=='pressure' else 0.;assert kind!='pressure' or skew>1e-3 and symmetrized_defect>1e-3
            rows.append(dict(case='curved_face_'+kind,tangent_error=derivative,work_error=work,wrench_error=wrench,tangent_asymmetry=skew,symmetrization_negative_control=symmetrized_defect))
        ftotal=np.zeros_like(q);Ktotal=np.zeros((q.size,q.size))
        def volume(qq):
            out=0.
            for N,Da,Db,w,area0 in faces:out+=np.sum(w*np.sum((N@qq)*np.cross(Da@qq,Db@qq),axis=1))/3
            return out
        for N,Da,Db,w,area0 in faces:
            f,K,_=law.surface(N,Da,Db,q,w,'pressure',2.3,area0);ftotal+=f;Ktotal+=K
        volume_work=float(-2.3*volume(q+1e-25j*v).imag/1e-25);work_error=abs(np.sum(ftotal*v)-volume_work)/max(abs(volume_work),1e-30);symmetry=rel(Ktotal-Ktotal.T,Ktotal);assert max(work_error,symmetry)<1e-8;rows.append(dict(case='closed_pressure_volume_potential',volume=float(volume(q)),volume_work_error=float(work_error),closed_tangent_asymmetry=symmetry))
        r=np.array([[.1,-.2,.05],[1.1,.08,-.07]]);R=np.array([beam.float_array(beam.exp_rotation([.2,-.1,.15])),beam.float_array(beam.exp_rotation([-.1,.25,.3]))]);delta=rng.normal(size=12)*.2
        for ff,mf in [('world','world'),('body','body')]:
            F=np.array([.8,-.3,.2]);M=np.array([.04,.08,-.03]);c=np.array([0.,.03,-.02]);f,K,info=law.point_wrench(r,R,.73,c,F,M,ff,mf)
            def load_change(z):
                rr,RR=beam.update(r,R,z);fn=law.point_wrench(rr,RR,.73,c,F,M,ff,mf)[0];return chart(z).T@fn
            err=rel(fd(load_change,delta)-K@delta,K@delta);physical=info['force_world']@(info['J']@delta)+info['moment_world']@(info['R']@(info['W']@delta));work=abs(f@delta-physical)/max(abs(physical),1e-30);assert max(err,work)<1e-8;rows.append(dict(case='eccentric_beam_wrench_'+ff,tangent_error=err,work_error=float(work),tangent_asymmetry=rel(K-K.T,K)))
        factor=law.viscosity_factor(np.diag([.4,.7,.9,.2,.3,.5]));u=q-q0;vel=rng.normal(size=q.shape)*.0003;fs,Kq,Kv,info=law.shell_damping(d,u,vel,factor);direction=rng.normal(size=q.shape)*.001
        du=law.shell_damping(d,u+1e-25j*direction,vel,factor,False)[0].imag/1e-25;dv=law.shell_damping(d,u,vel+1e-25j*direction,factor,False)[0].imag/1e-25;errors=[rel(du.ravel()-Kq@direction.ravel(),du),rel(dv.ravel()-Kv@direction.ravel(),dv)];power=abs(np.sum(fs*vel)-info['dissipation'])/max(abs(info['dissipation']),1e-30);Omega=np.array(beam.skew([.3,-.2,.4]),float);rigid_velocity=q@Omega.T+t[:,None]*np.array([.4,.1,-.2]);rigid=law.shell_damping(d,u,rigid_velocity,factor,False)[3];rigid_error=float(np.max(abs(rigid['strain_rate'])));eigen=float(np.linalg.eigvalsh(Kv).min());assert max(errors+[power,rigid_error])<1e-8 and eigen>-1e-9*np.linalg.norm(Kv);rows.append(dict(case='objective_shell_viscosity',state_tangent_error=errors[0],velocity_tangent_error=errors[1],power_error=float(power),rigid_motion_strain_rate=rigid_error,dissipation=info['dissipation'],velocity_tangent_min=eigen))
        velocity=rng.normal(size=12)*.02;fb,Kqb,Kvb,bi=law.beam_damping(r,R,1.,velocity,factor)
        def damp_change(z):
            rr,RR=beam.update(r,R,z);T=chart(z);fn=law.beam_damping(rr,RR,1.,T@velocity,factor)[0];return T.T@fn
        bq=rel(fd(damp_change,delta)-Kqb@delta,Kqb@delta);bv=rel(fd(lambda z:law.beam_damping(r,R,1.,velocity+z,factor)[0],delta)-Kvb@delta,Kvb@delta);bp=abs(fb@velocity-bi['dissipation'])/bi['dissipation'];rv=np.zeros(12)
        for i in range(2):rv[6*i:6*i+3]=np.cross([.3,-.2,.4],r[i])+[.4,.1,-.2];rv[6*i+3:6*i+6]=R[i].T@np.array([.3,-.2,.4])
        br=law.beam_damping(r,R,1.,rv,factor)[3];be=float(np.linalg.norm(br['strain_rate']));assert max(bq,bv,bp,be)<1e-8;rows.append(dict(case='objective_beam_viscosity',state_tangent_error=bq,velocity_tangent_error=bv,power_error=float(bp),rigid_motion_strain_rate=be))
        y=np.array([.3,-.2,.1,1.1,.2,-.1,-.1,1.,.15]);yv=rng.normal(size=9)*.1;ts=.4;Ks=np.diag([12.,9.,8.,3.,4.,5.]);Ds=np.diag([.4,.3,.2,.1,.2,.3])
        def support(yy,vv):
            e,J,H=support_map(yy,ts);edot=J[:,:9]@vv+J[:,9];eq=np.einsum('kij,j->ki',H[:,:9,:9],vv)+H[:,9,:9];out=law.strain_support(e,J[:,:9],H[:,:9,:9],edot,Ks,Ds,eq);return out,(e,J,H,edot)
        (sf,sq,sv,si),(e,J,H,edot)=support(y,yv);yd=rng.normal(size=9)*.1;sqe=rel(fd(lambda z:support(y+z,yv)[0][0],yd)-sq@yd,sq@yd);sve=rel(fd(lambda z:support(y,yv+z)[0][0],yd)-sv@yd,sv@yd);basepower=float((Ks@e+Ds@edot)@J[:,9]);identity=float(sf@yv-e@Ks@edot-si['dissipation']+basepower);assert max(sqe,sve,abs(identity))<1e-8;rows.append(dict(case='moving_elastic_viscous_frame_support',state_tangent_error=sqe,velocity_tangent_error=sve,power_balance_error=abs(identity),moving_base_power=basepower))
        # Superpose a time-dependent rigid observer rotation on a nonzero
        # deformation/rate. Internal viscous stress must transform objectively.
        observer=beam.float_array(beam.exp_rotation([-.3,.5,.1]));u2=q@observer.T-q0;vel2=vel@observer.T+q@observer.T@Omega.T;fo,_,_,oi=law.shell_damping(d,u2,vel2,factor,False);covariance=rel(fo-fs@observer.T,fs);assert covariance<1e-8;rows[-3]['superposed_observer_force_error']=covariance
        compressed=beam.evaluate(np.array([[0.,0,0],[.95,0,0]]),np.array([np.eye(3)]*2),1.,np.diag([50.,20.,20.,.2,.2,.2]),np.array([1.,0,0,0,0,0]));negative=float(np.linalg.eigvalsh(compressed['K'][6:,6:]).min());cf,_,cv,ci=law.beam_damping(np.array([[0.,0,0],[.95,0,0]]),np.array([np.eye(3)]*2),1.,velocity,factor);assert negative<0 and ci['dissipation']>=0;rows[-2]['compressed_elastic_minimum']=negative;rows[-2]['compressed_viscous_power']=ci['dissipation']
        rejected=0
        for invalid in [-np.eye(6),np.diag([1,1,1,1,1,-.01]),np.eye(6)+np.triu(np.ones((6,6)),1)]:
            try:law.viscosity_factor(invalid)
            except ValueError:rejected+=1
            else:raise AssertionError('invalid viscosity accepted')
        assert rejected==3;zero=law.viscosity_factor(np.diag([1.,0,0,0,0,0]));assert np.linalg.matrix_rank(zero)==1;rows.append(dict(case='viscosity_admissibility',negative_controls=rejected,semidefinite_rank=1,no_stiffness_based_artificial_viscosity=True))
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    print(json.dumps(rows,indent=2),flush=True)
    if '--publish' in sys.argv:
        dest=ROOT/'docs/verification/structural-load-damping-result.json';assert not dest.exists();paths=[Path(__file__),Path(law.__file__),Path(pg.__file__),Path(kin.__file__),Path(beam.__file__),Path(ad.__file__),ROOT/'tools/design-experiments/mass_paired_ancf.py',ROOT/'tools/design-experiments/pg_ritz_mechanics.py',ROOT/'docs/verification/structural-final-integration-protocol.md'];dest.write_text(json.dumps(dict(passed=len(rows)==10 and not failures,records=rows,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print('STRUCTURAL_LOAD_DAMPING_RESULT',len(rows)==10 and not failures,flush=True)
if __name__=='__main__':main()
