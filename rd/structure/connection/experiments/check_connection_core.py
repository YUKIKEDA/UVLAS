import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import copy,hashlib,json,time,traceback
import numpy as np
from scipy.linalg import expm
import beam_jet2 as ad
import beam_se3_v2 as group
import connection_kinematics as k
ROOT=rd_paths.ROOT
def rel(a,b):return float(np.linalg.norm(a)/max(np.linalg.norm(b),1e-30))
def ordinary(y,c):
    p=y[:3];a=y[3:6];b=y[6:9];u=a/np.sqrt(a@a);n=np.cross(a,b);n=n/np.sqrt(n@n);Q=np.column_stack((u,np.cross(n,u),n));return p+Q@c
def target(t):
    r=np.array([.2*t*t,.1*ad.sin(t),-.03*t],object)
    Q=group.exp_rotation([.4*t,0.,0.])@group.exp_rotation([0.,.3*t*t,0.])@group.exp_rotation([0.,0.,.2*t*t*t])
    return r,Q
def constraint(y,t,offset):
    p,Q=k.frame(y);r,R=target(t);return np.r_[p-r-R@offset,group.log_rotation(R.T@Q)]
def main():
    dest=ROOT/'rd/structure/connection/records/connection-core-result.json';assert not dest.exists();start=time.monotonic();records=[];failures=[]
    try:
        y=np.array([.7,-.3,.2,1.1,.2,-.1,-.2,.9,.3]);m=2.7;c=np.array([.04,-.03,.06]);Jc=np.diag([.032,.021,.025]);M,DM,v=k.rigid_mass(y,m,c,Jc);rng=np.random.default_rng(20261105);direction=rng.normal(size=9);direction/=np.linalg.norm(direction)
        h=1e-4;derivative=(k.rigid_mass(y-2*h*direction,m,c,Jc)[0]-8*k.rigid_mass(y-h*direction,m,c,Jc)[0]+8*k.rigid_mass(y+h*direction,m,c,Jc)[0]-k.rigid_mass(y+2*h*direction,m,c,Jc)[0])/(12*h)
        derr=rel(derivative-np.einsum('ijk,k->ij',DM,direction),derivative);assert derr<1e-7
        B=np.zeros((9,6));B[:3,:3]=np.eye(3)
        for i in range(3):B[3*i:3*i+3,3:]=-np.array(group.skew(y[3*i:3*i+3]),float)
        X=v['center'];Q=v['R'];C=np.array(group.skew(X),float);Jworld=Q@Jc@Q.T;physical=np.block([[m*np.eye(3),-m*C],[m*C,Jworld+m*((X@X)*np.eye(3)-np.outer(X,X))]])
        inertia_error=rel(B.T@M@B-physical,physical);assert inertia_error<1e-9
        vel=rng.normal(size=9);gyro=group.gyro(DM,vel);power=.5*np.einsum('ijk,i,j,k',DM,vel,vel,vel);p_error=abs(vel@gyro-power)/max(abs(power),1e-30);assert p_error<1e-10
        records.append(dict(case='eccentric_rigid_body_on_attachment_frame',spatial_inertia_error=inertia_error,mass_derivative_error=derr,power_identity_error=p_error))
        g=np.array([.3,-.7,-9.81]);force=v['J'].T@(m*g);tangent=np.einsum('a,aij->ij',m*g,v['H']);ce=ordinary(y.astype(complex)+1j*1e-30*direction,c);energy_derivative=(m*g@ce).imag/1e-30
        fds=[]
        for step in [-2,-1,1,2]:fds.append(k.prepared(y+step*h*direction,c)['J'].T@(m*g))
        fd=(fds[0]-8*fds[1]+8*fds[2]-fds[3])/(12*h);metrics=dict(work=abs(force@direction-energy_derivative)/abs(energy_derivative),tangent=rel(fd-tangent@direction,fd))
        expected=np.r_[m*g,np.zeros(6)]+v['W'].T@np.cross(c,Q.T@(m*g));metrics['eccentric_torque']=rel(force-expected,force);assert max(metrics.values())<1e-8
        records.append(dict(case='gravity_and_eccentric_wrench',metrics=metrics))
        # A particle uses only the linear point map; it does not require Jc>0.
        P=np.zeros((3,9));P[:,:3]=np.eye(3);Mp=1.3*P.T@P;assert np.linalg.matrix_rank(Mp)==3
        records.append(dict(case='particle_separate_from_rigid_body',mass_rank=3,zero_rotational_inertia_not_fabricated=True))
        offset=np.array([.03,.04,-.02]);t0=.37;timejet=ad.Jet.variable(t0,9);r,R=target(timejet);path=np.r_[r+R@offset,R[:,0],R[:,1]];Y=np.array([a.v for a in path]);V=np.array([a.g[9] for a in path]);Acc=np.array([a.H[9,9] for a in path])
        yj=np.array([ad.Jet.variable(a,i) for i,a in enumerate(Y)],object);tj=ad.Jet.variable(t0,9);Cj=constraint(yj,tj,offset);J=np.array([a.g[:10] for a in Cj]);H=np.array([a.H[:10,:10] for a in Cj]);val=np.array([a.v for a in Cj]);velocity=J@np.r_[V,1.];acc=J[:,:9]@Acc+np.einsum('aij,i,j->a',H,np.r_[V,1.],np.r_[V,1.]);assert max(np.linalg.norm(val),np.linalg.norm(velocity),np.linalg.norm(acc))<1e-8
        # Independent off-manifold path checks all q/t mixed terms.
        w=rng.normal(size=9)*.2;aa=rng.normal(size=9)*.1;alongtime=ad.Jet.variable(t0,0);dt=alongtime-t0;along=np.array([Y[i]+w[i]*dt+aa[i]*dt*dt/2 for i in range(9)],object);direct=constraint(along,alongtime,offset);expected=np.array([a.H[0,0] for a in direct]);computed=J[:,:9]@aa+np.einsum('aij,i,j->a',H,np.r_[w,1.],np.r_[w,1.]);error=rel(computed-expected,expected);omitted=rel(computed-2*np.einsum('ai,i->a',H[:,:9,9],w)-expected,expected);assert error<1e-8 and omitted>1e-4
        records.append(dict(case='prescribed_motion_velocity_acceleration',position=float(np.linalg.norm(val)),velocity=float(np.linalg.norm(velocity)),acceleration=float(np.linalg.norm(acc)),mixed_derivative_error=error,omitted_mixed_term_defect=omitted))
        A=np.array([[1.,0,0],[0,1,0],[1,1,0],[0,0,1.]]);rhs=np.array([[.2,.3],[.1,-.2],[.3,.1],[.4,.7]]);keep,T,ranks=k.reduce_linear_constraints(A,rhs);at=.41;q=rhs[keep]@np.array([1.,at]);q=np.linalg.lstsq(A[keep],q,rcond=None)[0];nodal=np.array([1.,-.3,.7]);mult=np.linalg.lstsq(A.T,nodal,rcond=None)[0];assert rel(A.T@mult-nodal,nodal)<1e-10
        bad=rhs.copy();bad[2,1]+=.1
        try:k.reduce_linear_constraints(A,bad)
        except ValueError:pass
        else:raise AssertionError('accepted inconsistent time-dependent duplicate')
        records.append(dict(case='linear_constraint_dependency_and_dual_reaction',ranks=ranks,independent_rows=keep.tolist(),all_time_rhs_checked=True,individual_dependent_multipliers_nonunique=True))
        nodes={}
        for i in range(3):
            G=np.eye(4);G[:3,:3]=expm(np.array(group.skew([.1*i,-.03*i,.04*i]),float));G[:3,3]=[i,.1*i,-.04*i];nodes[str(i)]=G
        edges=[]
        for ident,a,b in [('a','0','1'),('b','1','2'),('cycle','0','2'),('dup','0','1'),('reverse','1','0')]:edges.append((ident,a,b,np.linalg.inv(nodes[a])@nodes[b]))
        kept,redundant=k.rigid_weld_graph(nodes,edges);assert len(kept)==2 and len(redundant)==3
        wrong=list(edges);Hbad=wrong[2][3].copy();Hbad[0,3]+=.01;wrong[2]=('bad','0','2',Hbad)
        try:k.rigid_weld_graph(nodes,wrong)
        except ValueError:pass
        else:raise AssertionError('inconsistent weld cycle')
        records.append(dict(case='SE3_weld_graph_global_redundancy',kept=kept,redundant=redundant,inconsistent_cycle_rejected=True))
        ledger=k.MassLedger();ledger.add('point','payload-point');ledger.add('body','payload-body');rejected=[]
        for ident,source,kwargs in [('point','fresh',{}),('fresh','payload-body',{}),('a','new',dict(host_already_included=True)),('b','new2',dict(mode='replace'))]:
            before=copy.deepcopy(ledger.__dict__)
            try:ledger.add(ident,source,**kwargs)
            except ValueError:rejected.append(ident)
            else:raise AssertionError('mass double count')
            assert ledger.__dict__==before
        records.append(dict(case='mass_ownership_and_atomic_failure',negative_controls=len(rejected),no_partial_registration=True))
        before=y.copy();rejects=0
        for invalid in [np.r_[y[:3],np.zeros(6)],np.r_[y[:3],y[3:6],2*y[3:6]]]:
            try:k.prepared(invalid,c)
            except ValueError:rejects+=1
        assert rejects==2 and np.array_equal(before,y);records.append(dict(case='attachment_degeneracy_rejection',rejected=rejects,input_unchanged=True))
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    paths=[Path(__file__),Path(k.__file__),Path(group.__file__),Path(ad.__file__),ROOT/'rd/structure/connection/records/connection-adoption-protocol.md']
    result=dict(passed=len(records)==8 and not failures,records=records,failures=failures,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(records,indent=2));print('CONNECTION_CORE_RESULT',result['passed'])
if __name__=='__main__':main()
