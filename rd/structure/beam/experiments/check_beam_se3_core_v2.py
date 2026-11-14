"""SE(3) beam core, inertia, frame/work and independent Decimal checks."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,time,traceback
import numpy as np
from scipy.linalg import expm
import beam_se3_v2 as b
import beam_decimal_reference as independent
ROOT=rd_paths.ROOT
def rel(a,c):return float(np.linalg.norm(a)/max(np.linalg.norm(c),1e-300))
def make_pose(d):
    r,Q=b.pose_at(np.zeros(3),np.eye(3),np.array(d,dtype=object),1.)
    return np.array([np.zeros(3),b.float_array(r)]),np.array([np.eye(3),b.float_array(Q)])
def main():
    dest=ROOT/'rd/structure/beam/records/beam-se3-core-v2-result.json';assert not dest.exists();start=time.monotonic();records=[];failures=[];artifact=None
    try:
        L=1.2;ref=np.array([1.,.01,-.02,.12,-.08,.21]);e0=np.array([.001,-.002,.001,.003,-.001,.002]);e=np.array([.02,-.01,.015,.03,.07,-.04])
        rng=np.random.default_rng(20261105);A=rng.normal(size=(6,6));K=A.T@A+np.diag([20.,12.,10.,2.,3.,4.]);r,R=make_pose(L*(ref+e0+e));answer=b.evaluate(r,R,L,K,ref,e0)
        require=lambda condition,message:(_ for _ in ()).throw(AssertionError(message)) if not condition else None
        err=rel(answer['strain']-e,e);require(err<1e-10,'constant strain')
        r0,R0=make_pose(L*ref);zero=b.evaluate(r0,R0,L,K,ref,np.zeros(6));require(abs(zero['U'])<1e-24,'curved stress-free reference')
        records.append(dict(case='curved_twisted_reference_and_constant_strain',strain_error=err,reference_energy=zero['U']))
        d=rng.normal(size=12);d/=np.linalg.norm(d);fd,sd=independent.directional(r,R,L,K,ref,e0,d)
        metrics=dict(force=abs(fd-answer['f']@d)/abs(fd),tangent=abs(sd-d@answer['K']@d)/abs(sd),symmetry=rel(answer['K']-answer['K'].T,answer['K']))
        require(max(metrics.values())<1e-8,'independent Decimal80')
        records.append(dict(case='independent_Decimal80_energy_derivatives',metrics=metrics));print('DERIVATIVES',metrics,flush=True)
        Q=expm(np.array(b.skew([.4,-.7,.3]),float));translation=np.array([1.7,-.4,.8]);rr=r@Q.T+translation;RR=np.array([Q@a for a in R]);rot=b.evaluate(rr,RR,L,K,ref,e0)
        transform=np.eye(12)
        for i in [0,6]:transform[i:i+3,i:i+3]=Q
        metrics=dict(energy=abs(rot['U']-answer['U'])/abs(answer['U']),force=rel(transform.T@rot['f']-answer['f'],answer['f']),tangent=rel(transform.T@rot['K']@transform-answer['K'],answer['K']))
        require(max(metrics.values())<1e-10,'rigid objectivity');records.append(dict(case='finite_global_rigid_transform',metrics=metrics))
        # Fixed section-origin offset: conjugation, strain dual map and full AD.
        offset=np.array([.01,.04,-.03]);C=np.array(b.skew(offset),float);S=np.block([[np.eye(3),-C],[np.zeros((3,3)),np.eye(3)]]);Si=np.linalg.inv(S);Ks=Si.T@K@Si
        rs,Rs=b.jets(r,R);shifted=[rs[i]+Rs[i]@offset for i in range(2)];Us,es,ds=b.energy_from_poses(shifted,Rs,L,Ks,S@ref,S@e0)
        metrics=dict(energy=abs(Us.v-answer['U'])/abs(answer['U']),force=rel(Us.g-answer['f'],answer['f']),tangent=rel(Us.H-answer['K'],answer['K']))
        require(max(metrics.values())<1e-10,'section origin');records.append(dict(case='section_origin_dual_transform',metrics=metrics))
        # Unstressed free beam has only six rigid modes (strain square root).
        straight=np.array([[0.,0.,0.],[L,0.,0.]]);eye=np.array([np.eye(3)]*2);base=b.evaluate(straight,eye,L,K,np.array([1.,0,0,0,0,0]));root=np.linalg.cholesky(K).T@base['B']*np.sqrt(L)
        singular=np.linalg.svd(root,compute_uv=False);ranks=[int(sum(singular>tol*singular[0])) for tol in [1e-9,1e-10,1e-11]];require(ranks==[6]*3,'rank')
        rigid=np.zeros((12,6))
        for i in range(2):
            rigid[6*i:6*i+3,:3]=np.eye(3);rigid[6*i:6*i+3,3:]=-np.array(b.skew(straight[i]),float);rigid[6*i+3:6*i+6,3:]=np.eye(3)
        require(np.linalg.norm(root@rigid)<1e-10,'rigid six');records.append(dict(case='free_element_rigid_six',root_ranks=ranks,nullity=12-ranks[0]))
        # Linear density and physical centroid/inertia, not arbitrary SPD mass.
        mass=b.spatial_mass(2.3,np.array([0.,.037,-.024]),np.diag([.023,.009,.014]));M,dM=b.inertia(r,R,L,mass,12);M2,_=b.inertia(r,R,L,mass,18)
        require(np.linalg.eigvalsh(M).min()>0,'positive mass');quad=rel(M-M2,M2);require(quad<1e-10,'mass quadrature')
        rigid=np.zeros((12,6))
        for i in range(2):
            rigid[6*i:6*i+3,:3]=np.eye(3);rigid[6*i:6*i+3,3:]=-np.array(b.skew(r[i]),float);rigid[6*i+3:6*i+6,3:]=R[i].T
        # Integrate physical rigid-section velocities independently at every s.
        nodes,weights=np.polynomial.legendre.leggauss(20);direct=np.zeros((6,6));twist=b.relative(r[0],R[0],r[1],R[1])
        for t,w in zip((nodes+1)/2,weights*L/2):
            x,Qp=b.pose_at(r[0],R[0],twist,t);x=b.float_array(x);Qp=b.float_array(Qp);B=np.block([[Qp.T,-Qp.T@np.array(b.skew(x),float)],[np.zeros((3,3)),Qp.T]]);direct+=w*B.T@mass@B
        error=rel(rigid.T@M@rigid-direct,direct);require(error<1e-9,'physical inertia')
        records.append(dict(case='consistent_eccentric_mass_and_physical_inertia',quadrature_error=quad,spatial_inertia_error=error,mass_min=float(np.linalg.eigvalsh(M).min())))
        v=rng.normal(size=12)*.07;gyro=b.gyro(dM,v);rhs=.5*np.einsum('ijk,i,j,k',dM,v,v,v);identity=abs(v@gyro-rhs)/max(abs(rhs),1e-30)
        eps=2e-5;rp,Rp=b.update(r,R,eps*v);rm,Rm=b.update(r,R,-eps*v)
        # These mass matrices live in new body frames: pull them back to the
        # fixed anchor coordinate chart using the right SO(3) Jacobian.
        def pulled(rr,rot,delta):
            MM,_=b.inertia(rr,rot,L,mass,12);T=np.eye(12)
            for i in range(2):T[6*i+3:6*i+6,6*i+3:6*i+6]=b.float_array(b.left_jacobian(-delta[6*i+3:6*i+6]))
            return T.T@MM@T
        diff=(pulled(rp,Rp,eps*v)-pulled(rm,Rm,-eps*v))/(2*eps);analytic=np.einsum('ijk,k->ij',dM,v);err=rel(diff-analytic,analytic)
        require(identity<1e-10 and err<1e-7,'configuration mass and gyro');records.append(dict(case='configuration_mass_and_Christoffel',power_identity_error=identity,mass_derivative_error=err))
        eta=np.array([.012,-.04,.023]);x,Qp=b.point_jet(r,R,.37,eta);pointJ=np.array([q.g for q in x]);pointH=np.array([q.H for q in x]);force=np.array([.7,-.4,.3])
        nodal=pointJ.T@force;d=rng.normal(size=12);work=abs(nodal@d-force@(pointJ@d));require(work<1e-10,'transpose work')
        # A common rigid displacement yields a known wrench about the global origin.
        wrench=rigid.T@nodal;X=b.float_array(x);expected=np.r_[force,np.cross(X,force)];require(rel(wrench-expected,expected)<1e-10,'eccentric point wrench')
        records.append(dict(case='point_trace_J_H_and_transpose_work',work_error=work,wrench_error=rel(wrench-expected,expected)))
        invalid=[]
        def reject(name,fn):
            try:fn()
            except ValueError:invalid.append(name)
            else:raise AssertionError('accepted '+name)
        reject('zero_length',lambda:b.evaluate(r,R,0.,K,ref));reject('negative_length',lambda:b.evaluate(r,R,-1.,K,ref));bad=R.copy();bad[0,0,0]*=2;reject('rotation',lambda:b.evaluate(r,bad,L,K,ref))
        near=np.array([np.eye(3),expm(np.array(b.skew([np.pi-1e-5,0,0]),float))]);reject('relative_branch',lambda:b.evaluate(r,near,L,K,ref))
        bad=K.copy();bad[0,0]=-1;reject('nonSPD',lambda:b.evaluate(r,R,L,bad,ref));bad=K.copy();bad[0,1]+=.1;reject('asymmetric',lambda:b.evaluate(r,R,L,bad,ref))
        reject('mass',lambda:b.spatial_mass(0.,np.zeros(3),np.eye(3)));reject('inertia_triangle',lambda:b.spatial_mass(1.,np.zeros(3),np.diag([3.,1.,1.])));reject('initial_strain',lambda:b.evaluate(r,R,L,K,ref,[float('nan')]*6))
        records.append(dict(case='invalid_input_controls',rejected=invalid,count=len(invalid)))
        out=ROOT/'build/analysis'/('beam-se3-core-v2-'+time.strftime('%Y%m%d-%H%M%S')+'.npz');np.savez_compressed(out,r=r,R=R,L=L,K_section=K,reference=ref,e0=e0,force=answer['f'],tangent=answer['K'],mass=M,mass_derivative=dM,point_J=pointJ,point_H=pointH)
        artifact=dict(artifact=out.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    paths=[Path(__file__),ROOT/'rd/structure/beam/experiments/check_beam_se3_core.py',ROOT/'rd/structure/beam/records/beam-se3-core-result.json',Path(b.__file__),Path(b.ad.__file__),Path(independent.__file__),ROOT/'rd/structure/beam/records/beam-se3-adoption-protocol.md']
    result=dict(passed=len(records)==9 and not failures,records=records,failures=failures,artifact=artifact,seconds=time.monotonic()-start,
        sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('BEAM_CORE_RESULT',result['passed'],flush=True)
if __name__=='__main__':main()
