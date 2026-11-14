import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np
import json,hashlib,time,traceback,sys
import uvlm_analytic_induction as candidate
import uvlm_induction_reference as reference
import check_uvlm_wp1_three_dimensional as old
import check_uvlm_adaptive_quadrature as adaptive
import check_vortex_core_load as core
ROOT=rd_paths.ROOT
def err(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30))
def qualified(a,b):return np.linalg.norm(a-b)<=1e-10 or err(a,b)<=1e-8
def main():
    start=time.monotonic();rows=[];failures=[]
    try:
        tests=[]
        for shape,quad in old.SHAPES.items():
            for eps in [0.,.02,.08]:
                for point in [old.POINTS[0],old.POINTS[1]]:
                    v,J=candidate.patch(point,quad,eps);low=np.array([old.area_reference(quad,mu,point[None,:],eps,64)[0] for mu in [[1.,0.],[0.,1.]]]);high=np.array([old.area_reference(quad,mu,point[None,:],eps,96)[0] for mu in [[1.,0.],[0.,1.]]]);assert qualified(low,high) and qualified(v,high);tests.append(dict(shape=shape,eps=eps,relative_error=err(v,high),reference_change=err(low,high)))
        rows.append(dict(case='finite_3D_three_shapes_and_all_LR_columns',cases=tests));print('FAR_SHAPES_OK',len(tests),flush=True)
        # Actual 3D triangle, all THREE P1 columns, no constant-mu shortcut.
        x=np.array([[0.,0.,0.],[1.,.1,.07],[.2,1.,-.05]]);normal=np.cross(x[1]-x[0],x[2]-x[0]);normal/=np.linalg.norm(normal);near=[]
        for beta,label in [(np.array([.3,.4,.3]),'interior'),(np.array([.5,.5,0.]),'edge'),(np.array([1.,0,0.]),'vertex')]:
            for eps,d in [(0.,.01),(0.,.001),(0.,.00001),(.02,0.),(.08,0.)]:
                p=beta@x+d*normal;v,J,_=candidate.triangle(p,x,eps,incident=(d==0));r1=reference.integrate(p,x,eps,48,True);r2=reference.integrate(p,x,eps,64,True)
                if not qualified(r1,r2):print('REFERENCE_UNRESOLVED',label,eps,d,err(r1,r2),flush=True)
                assert qualified(r1,r2) and qualified(v,r2);near.append(dict(location=label,eps=eps,distance=d,relative_error=err(v,r2),absolute_error=float(np.linalg.norm(v-r2)),reference_change=err(r1,r2)))
        rows.append(dict(case='general_triangle_near_and_declared_self_each_P1_basis',cases=near));print('NEAR_OK',len(near),flush=True)
        jacrows=[]
        for beta,d in [(np.array([.3,.4,.3]),.04),(np.array([.3,.4,.3]),0.),(np.array([1.,0,0.]),0.)]:
            p=beta@x+d*normal;eps=.08;v,J,_=candidate.triangle(p,x,eps,incident=(d==0));lo=reference.jacobian(p,x,eps,128);hi=reference.jacobian(p,x,eps,192);change=err(lo,hi);difference=err(J,hi);assert qualified(lo,hi) and qualified(J,hi)
            # Derivative along an incident point p=beta*x must include BOTH
            # point and vertex contributions. Check using the separate area rule.
            rng=np.random.default_rng(2101);dx=rng.normal(size=(3,3))*.03;dp=beta@dx;dir=np.r_[dp,dx.ravel(),.01];actual=J@dir;z=np.r_[p,x.ravel(),eps].astype(complex)+1e-25j*dir;ind=reference.integrate(z[:3],z[3:12].reshape(3,3),z[-1],192).imag/1e-25;chain=err(actual,ind);assert qualified(actual,ind);jacrows.append(dict(distance=d,all_13_derivative_error=difference,reference_derivative_change=change,incident_chain_rule_error=chain))
        rows.append(dict(case='point_vertices_core_Jacobian',cases=jacrows));print('JACOBIAN_OK',flush=True)
        # Constant basis sum is the independently evaluated analytic ring.
        p=np.array([.3,.2,.1]);v,J,_=candidate.triangle(p,x,.02);r=sum(core.kernel(p,x[i],x[j],.02)[0,0] for i,j in [(0,2),(2,1),(1,0)]);constant=err(v.sum(axis=0),r);assert constant<1e-10
        incident=x[0];v0,_,_=candidate.triangle(incident,x,.02,True);self_value=np.linalg.norm(v0.sum(axis=0));assert self_value>1e-3;rows.append(dict(case='constant_ring_and_nonzero_whole_self',ring_error=constant,whole_ring_self_norm=float(self_value)))
        transform=[];R=np.array([[np.cos(.4),-np.sin(.4),0],[np.sin(.4),np.cos(.4),0],[0,0,1.]])@np.array([[np.cos(.3),0,np.sin(.3)],[0,1.,0],[-np.sin(.3),0,np.cos(.3)]]);shift=np.array([.3,-.2,.5]);eps=.02;v,J,_=candidate.triangle(p,x,eps)
        for scale in [1e-3,1.,1e3]:
            z,Z,_=candidate.triangle(scale*(R@p+shift),scale*(x@R.T+shift),scale*eps);cov=err(scale*z,v@R.T);T=np.zeros((13,13))
            for j in range(4):T[3*j:3*j+3,3*j:3*j+3]=R
            T[-1,-1]=1.;expected=np.einsum('ab,kbc,dc->kad',R,J,T);jc=err(scale*scale*Z,expected);assert max(cov,jc)<1e-8;transform.append(dict(scale=scale,velocity_error=cov,Jacobian_error=jc))
        rv,_,_=candidate.triangle(p,x[[0,2,1]],eps);reverse=err(rv[[0,2,1]],-v);assert reverse<1e-10;rows.append(dict(case='rigid_scale_orientation_covariance',cases=transform,orientation_error=reverse))
        # In-plane exterior limit and smooth far-field branch.
        exterior=x[0]+1.5*(x[1]-x[0])-.4*(x[2]-x[0]);v,J,diag=candidate.triangle(exterior,x,0.);r=reference.integrate(exterior,x,0.,96);jref=reference.jacobian(exterior,x,0.,128);assert qualified(v,r) and qualified(J,jref);far=[]
        for distance in [8.,1e3,1e6]:
            pp=np.array([distance,.3*distance,.2*distance]);v,J,diag=candidate.triangle(pp,x,.02);r=reference.integrate(pp,x,.02,48);assert qualified(v,r);far.append(dict(distance=distance,branch=diag['branch'],relative_error=err(v,r)))
        rows.append(dict(case='coplanar_exterior_and_far_branches',exterior_derivative_error=err(candidate.triangle(exterior,x,0.)[1],jref),cases=far))
        # Shared vertices, two spanwise traces and a shared diagonal. Compare
        # to separately aggregated edge integrals plus independent sheet rules.
        vertices=np.array([[1.+.4*r,-.6+.6*j,.03*r*j] for r in range(3) for j in range(3)]);coeff=np.array([[1.,.2],[.7,-.4],[.3,.5]]);cells=old.fixture_triangles(vertices,None,coeff);pp=np.array([.7,.1,.3]);eps=.08;total=np.zeros(3);DJ=np.zeros((3,3+vertices.size+1))
        for ids,mu in cells:
            v,j,_=candidate.triangle(pp,vertices[ids],eps);total+=mu@v;g=np.einsum('k,kaj->aj',mu,j);DJ[:,:3]+=g[:,:3];DJ[:,-1]+=g[:,-1]
            for local,globalid in enumerate(ids):DJ[:,3+3*globalid:6+3*globalid]+=g[:,3+3*local:6+3*local]
        edges=old.edge_assembly(cells,[eps]*len(cells));aggregated=np.zeros(3)
        for (i,j),edge in edges.items():aggregated+=old.line(pp[None,:],vertices[i],vertices[j],*edge['mu'],eps,96)[0]
        for ids,mu in cells:aggregated+=old.triangle_primary(pp[None,:],vertices[ids],mu,eps,surface_order=96,boundary=False)[0]
        shared=err(total,aggregated);assert shared<1e-8;rng=np.random.default_rng(2102);delta=rng.normal(size=DJ.shape[1])*.03;zz=np.r_[pp,vertices.ravel(),eps].astype(complex)+1e-25j*delta;ref=np.zeros(3)
        for ids,mu in cells:ref+=(mu@reference.integrate(zz[:3],zz[3:-1].reshape(-1,3)[ids],zz[-1],96)).imag/1e-25
        chain=err(DJ@delta,ref);assert chain<1e-8;averaged=np.repeat(coeff.mean(axis=1)[:,None],2,axis=1);badcells=old.fixture_triangles(vertices,None,averaged);bad=sum(mu@candidate.triangle(pp,vertices[ids],eps)[0] for ids,mu in badcells);averaging_defect=err(bad,total);assert averaging_defect>1e-3
        node_velocity=np.zeros(3);node_reference=np.zeros(3)
        for ids,mu in cells:
            v,_,_=candidate.triangle(vertices[4],vertices[ids],eps,incident=(4 in ids));node_velocity+=mu@v;node_reference+=mu@reference.integrate(vertices[4],vertices[ids],eps,64,True)
        node_error=err(node_velocity,node_reference);assert node_error<1e-8;rows.append(dict(case='shared_edge_spanwise_traces_and_incident_mesh_node',aggregated_velocity_error=shared,shared_vertex_chain_rule_error=chain,trace_averaging_negative_control=averaging_defect,incident_node_velocity_error=node_error))
        rejected=[];before=x.copy()
        for label,fun in [('ideal_self',lambda:candidate.triangle(x[0],x,0.,True)),('nonincident_contact',lambda:candidate.triangle(x[0],x,.02)),('false_incidence',lambda:candidate.triangle(p,x,.02,True)),('negative_core',lambda:candidate.triangle(p,x,-.02)),('degenerate',lambda:candidate.triangle(p,np.zeros((3,3)),.02)),('patch_false_incidence',lambda:candidate.patch(p+normal,old.SHAPES['flat'],.02,True)),('nonfinite_point',lambda:candidate.triangle(np.array([np.nan,0,0]),x,.02)),('core_underflow',lambda:candidate.triangle(x[0],x,1e-200,True)),('inconsistent_shared_core',lambda:old.edge_assembly(cells,[.08]+[.02]*(len(cells)-1)))]:
            try:fun()
            except ValueError:rejected.append(label)
            else:raise AssertionError('negative control accepted '+label)
        assert len(rejected)==9 and np.array_equal(x,before);rows.append(dict(case='classification_and_failure_guards',rejected=rejected,input_unchanged=True))
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    print(json.dumps(rows,indent=2),flush=True)
    if '--publish' in sys.argv:
        dest=ROOT/'rd/aero/kernel/records/uvlm-induction-adoption-result.json';assert not dest.exists();artifact=None
        if not failures:
            B,G,_=candidate.triangle(p,x,.02);out=ROOT/'build/analysis/uvlm-induction-snapshot.npz';assert not out.exists();np.savez_compressed(out,point=p,vertices=x,eps=.02,basis=B,Jacobian=G,mesh_vertices=vertices,mesh_coefficients=coeff,mesh_point=pp,mesh_velocity=total,mesh_Jacobian=DJ);artifact=dict(artifact=out.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
        paths=[Path(__file__),Path(candidate.__file__),Path(reference.__file__),ROOT/'rd/aero/kernel/experiments/uvlm_dual13.py',Path(old.__file__),Path(adaptive.__file__),Path(core.__file__),ROOT/'rd/aero/kernel/records/uvlm-induction-adoption-protocol.md'];dest.write_text(json.dumps(dict(passed=len(rows)==8 and not failures,records=rows,failures=failures,artifact=artifact,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print('UVLM_INDUCTION_ADOPTION_RESULT',len(rows)==8 and not failures,flush=True)
if __name__=='__main__':main()
