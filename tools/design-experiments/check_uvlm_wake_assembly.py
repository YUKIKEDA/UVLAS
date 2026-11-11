from pathlib import Path
from dataclasses import replace
import numpy as np
import json,hashlib,time,traceback,sys,copy
import uvlm_wake_assembly as aero
import uvlm_surface_bindings as binding
import uvlm_geometry_guard as guard
import uvlm_induction_reference as reference
import check_uvlm_wp1_three_dimensional as old
import connection_shell as pg
ROOT=Path(__file__).resolve().parents[2]
def rel(a,b):return float(np.linalg.norm(a)/max(np.linalg.norm(b),1e-30))
def fd(fun,h=.001):return (-fun(2*h)+8*fun(h)-8*fun(-h)+fun(-2*h))/(12*h)

def fixture():
    b1=binding.q1('main',2,2,[[0.,-.6,0.],[.08,.6,.03],[1.08,.6,.07],[1.,-.6,0.]])
    b2=binding.q1('tail',2,2,[[1.6,-.4,.25],[1.64,.4,.29],[2.34,.4,.27],[2.3,-.4,.25]])
    model=aero.Model([b1,b2],.025);q=np.vstack([b1['q0'],b2['q0']]);flow=aero.Flow([1.,0,.025],[0.,0,.004],[1.,0,0],1.);return model,q,flow

def independent_basis(model,g,wakes,p,order=96):
    B=np.zeros((3,model.nb));sizes=[w.mu.size for w in wakes];offset=np.cumsum([0]+sizes);W=np.zeros((3,sum(sizes)))
    for i,(a,b) in enumerate(model.edges):v=old.line(p[None,:],g['X'][a],g['X'][b],1.,1.,model.eps,order)[0];B+=v[:,None]*model.D2[i]
    for iw,idx,ids,x,name,r,j,tri in model.wake_faces(wakes):
        v=reference.integrate(p,x,model.eps,order,True)
        for k,col in enumerate(idx):W[:,offset[iw]+col]+=v[k]
    return B,W

def main():
    start=time.monotonic();rows=[];failures=[];artifact=None
    try:
        model,q,flow=fixture();parent=aero.cold_state(model,q);fingerprint=parent.fingerprint();workspace=aero.Workspace();trial=aero.begin_trial(model,parent,q,parent.v,.12,flow);aero.solve_trial(model,parent,trial,flow,workspace);g=model.geometry(q,parent.v);A,W=workspace.assemble(model,g,trial['wakes']);Bref=[];Wref=[]
        for p,n in zip(g['P'],g['n']):b,w=independent_basis(model,g,trial['wakes'],p,64);Bref.append(n@b);Wref.append(n@w)
        Bref=np.array(Bref);Wref=np.array(Wref);Aref=Bref.copy();offset=0;known=[]
        for w,info in zip(trial['wakes'],model.surface_info):Aref+=Wref[:,offset:offset+info['n']]@model.T[info['mu_slice']];mu=w.mu.copy();mu[0]=0;known.extend(mu.ravel());offset+=w.mu.size
        rhs=np.sum(g['n']*(g['CV']-flow.value(g['P'],trial['time'])),axis=1)-Wref@known;independent_gamma=np.linalg.solve(Aref,rhs);Ae=rel(A-Aref,Aref);We=rel(W-Wref,Wref);Ge=rel(trial['gamma']-independent_gamma,independent_gamma);assert max(Ae,We,Ge)<1e-8 and trial['residual']<1e-11;rows.append(dict(case='finite_two_surface_independent_assembly',matrix_error=Ae,wake_basis_error=We,circulation_error=Ge,original_residual=trial['residual']));print('ASSEMBLY_OK',rows[-1],flush=True)
        independent=np.zeros_like(A)
        for s in [slice(0,4),slice(4,8)]:independent[s,s]=A[s,s]
        wrong=np.linalg.solve(independent,np.sum(g['n']*(g['CV']-flow.value(g['P'],trial['time'])),axis=1));coupling=rel(wrong-trial['gamma'],trial['gamma']);wrongfront=np.linalg.solve(Bref,rhs);frontdefect=rel(wrongfront-trial['gamma'],trial['gamma']);assert min(coupling,frontdefect)>.001;assert not np.any(model.D1@model.D2);rows.append(dict(case='integer_topology_mutual_and_implicit_front',integer_boundary_zero=True,dropped_mutual_defect=coupling,old_front_defect=frontdefect))
        # Integrate the independently evaluated physical edge velocity and
        # triangular load density, without consuming the candidate load packet.
        Fref=np.zeros(3);Mref=np.zeros(3);Qref=np.zeros((model.ns,3));kappa=model.D2@trial['gamma'];mu=np.concatenate([w.mu.ravel() for w in trial['wakes']])
        for edge in model.force_edges:
            a,b=model.edges[edge];p=(g['X'][a]+g['X'][b])/2;B0,W0=independent_basis(model,g,trial['wakes'],p,96);velocity=flow.value(p,trial['time'])+B0@trial['gamma']+W0@mu-(g['V'][a]+g['V'][b])/2;f=model.rho*kappa[edge]*np.cross(velocity,g['X'][b]-g['X'][a]);Fref+=f;Mref+=np.cross(p,f);Qref+=((model.Nv[a]+model.Nv[b])/2)[:,None]*f
        for panel,face in enumerate(model.faces):
            for tri in [[0,2,1],[0,3,2]]:
                ids=face[tri];points,weights,bary,n=old.quadrature(g['X'][ids],5);density=model.rho*trial['gamma_dot'][panel]*n;Fref+=sum(weights)*density;Mref+=np.sum(weights[:,None]*np.cross(points,density),axis=0);Qref+=(bary@model.Nv[ids]).T@(weights[:,None]*density)
        packet=trial['loads'];errors=[rel(packet['F']-Fref,Fref),rel(packet['M']-Mref,Mref),rel(packet['Q']-Qref,Qref)];assert max(errors)<1e-8;rows.append(dict(case='independent_all_vector_load_and_material_integrals',force_error=errors[0],moment_error=errors[1],generalized_error=errors[2]));print('INDEPENDENT_LOAD_OK',errors,flush=True)
        repeated=aero.solve_trial(model,parent,trial,flow,workspace);again=aero.begin_trial(model,parent,q,parent.v,.12,flow);aero.solve_trial(model,parent,again,flow,workspace);assert repeated['state'].fingerprint()==again['state'].fingerprint() and parent.fingerprint()==fingerprint;rows.append(dict(case='repeat_from_same_parent',identical_state=True,parent_unchanged=True,cache_hits=workspace.hits))
        state=aero.accept(model,parent,trial);oldg=model.geometry(state.q,state.v);cv=aero.convection(model,state,flow);errors=[]
        for iw,w in enumerate(state.wakes):
            for r,j in [(0,1),(1,1)]:
                b,a=independent_basis(model,oldg,state.wakes,w.vertices[r,j],96);mu=np.concatenate([ww.mu.ravel() for ww in state.wakes]);expected=flow.value(w.vertices[r,j],state.time)+b@state.gamma+a@mu;errors.append(rel(cv[iw][r,j]-expected,expected))
        assert max(errors)<1e-8;prescribed=aero.convection(model,state,flow,'prescribed');difference=max(rel(a-b,a) for a,b in zip(cv,prescribed));assert difference>.001;next_trial=aero.begin_trial(model,state,q,state.v,.12,flow);assert all(np.array_equal(wnew.mu[1:],wold.mu) for wnew,wold in zip(next_trial['wakes'],state.wakes));assert all(wnew.vertex_ids[2:]==wold.vertex_ids[1:] and wnew.vertex_ids[1]!=wold.vertex_ids[0] for wnew,wold in zip(next_trial['wakes'],state.wakes));rows.append(dict(case='all_source_free_convection_and_material_transport',independent_node_velocity_errors=errors,omitted_induction_defect=difference,old_mu_exactly_transported=True,detached_ids_distinct=True));print('CONVECTION_OK',flush=True)
        aero.solve_trial(model,state,next_trial,flow,workspace);rng=np.random.default_rng(2201);dq=rng.normal(size=q.shape)*.01;dv=rng.normal(size=q.shape)*.002;dU=np.array([.003,-.001,.002]);jvp=aero.trial_jvp(model,state,next_trial,flow,dq,dv,dU,workspace)
        def response(s):
            ff=aero.Flow(flow.U+s*dU,flow.amplitude,flow.wave,flow.omega);tt=aero.begin_trial(model,state,q+s*dq,state.v+s*dv,.12,ff);aero.solve_trial(model,state,tt,ff);return np.r_[tt['gamma'],tt['loads']['Q'].ravel()]
        numeric=fd(response,.002);analytic=np.r_[jvp['dgamma'],jvp['dQ'].ravel()];derror=rel(numeric-analytic,numeric);assert derror<1e-8;rows.append(dict(case='full_implicit_circulation_and_load_directional_tangent',relative_error=derror,includes_normals_velocity_inflow_and_old_time_convection=True));print('TANGENT_OK',derror,flush=True)
        # All independent parent variables, with the TE incidence chains intact.
        dw=[rng.normal(size=w.vertices.shape)*.003 for w in state.wakes];dgamma=rng.normal(size=state.gamma.shape)*.003;dmus=[];old_direction=model.geometry(state.q,state.v,dq)
        for iw,(w,info) in enumerate(zip(state.wakes,model.surface_info)):dw[iw][0]=old_direction['dX'][info['te']];dmu=rng.normal(size=w.mu.shape)*.003;dmu[0]=dgamma[info['gamma_te']];dmus.append(dmu)
        dc=aero.convection_jvp(model,state,flow,dq,dgamma,dw,np.concatenate([a.ravel() for a in dmus]));
        def conv(s):
            wakes=tuple(aero.Wake(w.sid,w.vertices+s*d,w.mu+s*mu,w.vertex_ids,w.row_ids,w.birth) for w,d,mu in zip(state.wakes,dw,dmus));ss=replace(state,q=state.q+s*dq,gamma=state.gamma+s*dgamma,wakes=wakes);return np.concatenate([a.ravel() for a in aero.convection(model,ss,flow)])
        numeric=fd(conv,.001);analytic=np.concatenate([a.ravel() for a in dc]);derror=rel(numeric-analytic,numeric);assert derror<1e-8;rows.append(dict(case='free_convection_full_state_directional_tangent',relative_error=derror));print('PARENT_TANGENT_OK',derror,flush=True)
        # Same physical point maps supply work, force and moment, including
        # the already adopted PG-ANCF coordinate basis.
        pm=pg.flat_model(1,.08,1.,1.);data=pg.data(pm,order=13);pb=binding.pg_shell('pg',2,2,pm,data['q0']);pmodel=aero.Model([pb],.025);ps=aero.cold_state(pmodel,pb['q0']);pt=aero.begin_trial(pmodel,ps,ps.q,ps.v,.12,flow);aero.solve_trial(pmodel,ps,pt,flow);gg=pmodel.geometry(ps.q);force=pt['loads'];testv=data['q0']@np.array([[.001,.002,0],[0,-.001,.001],[.002,0,.001]]).T;physical=float(np.sum(force['Fv']*(pmodel.Nv@testv)));general=float(np.sum(force['Q']*testv));work=abs(physical-general)/max(abs(physical),1e-30);Fdiff=rel(data['translation']@force['Q']-force['F'],force['F']);Mdiff=rel(np.cross(ps.q,force['Q']).sum(axis=0)-force['M'],force['M']);assert max(work,Fdiff,Mdiff)<1e-10
        pdir=rng.normal(size=ps.q.shape);grad=np.einsum('ia,pib->pab',pdir,data['G']);pdir*=.01/np.max(abs(grad));pj=aero.trial_jvp(pmodel,ps,pt,flow,pdir)['dQ']
        def pgload(s):
            trialp=aero.begin_trial(pmodel,ps,ps.q+s*pdir,ps.v,.12,flow);aero.solve_trial(pmodel,ps,trialp,flow);return trialp['loads']['Q']
        pnum=fd(pgload,.002);pjerr=rel(pnum-pj,pnum);assert pjerr<1e-8;rows.append(dict(case='adopted_PG_material_work_force_and_moment',work_error=work,force_error=Fdiff,moment_error=Mdiff,full_PG_load_tangent_error=pjerr,all_load_points_derived_from_same_material_map=True));print('PG_WORK_OK',pjerr,flush=True)
        # Finite window records retired material identities rather than silently
        # treating the removed vortex field as zero-error.
        third=aero.accept(model,state,next_trial);trim=aero.begin_trial(model,third,q,third.v,.12,flow,max_intervals=1);assert len(trim['retired'])==4 and trim['tail']=='finite_closed_window' and all(len(w.vertices)==2 for w in trim['wakes']);aero.solve_trial(model,third,trim,flow);growing=aero.begin_trial(model,third,q,third.v,.12,flow);aero.solve_trial(model,third,growing,flow);tail_difference=rel(trim['gamma']-growing['gamma'],growing['gamma']);assert tail_difference>1e-5;trim_state=aero.accept(model,third,trim);trim_restart=aero.restore(model,aero.checkpoint(model,trim_state),flow);resumed=aero.begin_trial(model,trim_restart,q,trim_restart.v,.12,flow);assert resumed['max_intervals']==1 and resumed['mode']=='free' and all(len(w.vertices)==2 for w in resumed['wakes']);rows.append(dict(case='explicit_finite_tail_and_retirement',retired_rows=len(trim['retired']),circulation_difference_from_growing=tail_difference,tail_error_not_assumed_zero=True,closure_and_convection_policy_restored=True))
        # The actual two-surface state and its derivative history survive restart.
        encoded=aero.checkpoint(model,third);restored=aero.restore(model,encoded,flow);assert restored.fingerprint()==third.fingerprint();assert np.array_equal(restored.gamma_dot(),next_trial['gamma_dot']);rr=aero.begin_trial(model,restored,q,restored.v,.12,flow);aero.solve_trial(model,restored,rr,flow);assert rr['state'].fingerprint()==growing['state'].fingerprint();assert np.array_equal(rr['loads']['Q'],growing['loads']['Q']);rows.append(dict(case='restart_and_current_load_history',bitwise_state_and_next_load_match=True,both_bound_circulation_levels_retained=True))
        # Epoch owner makes acceptance once-only. Evaluate changed trial objects
        # and incompatible input/checkpoint requests without changing state.
        store=aero.Store(model,third);before=store.state.fingerprint();rejected=[]
        bad=copy.deepcopy(growing);bad['q'][0,0]+=.01
        for label,call in [('post_solve_change',lambda:store.accept(bad)),('wrong_input',lambda:aero.restore(model,encoded,aero.Flow([1,0,0]))),('missing_input',lambda:aero.restore(model,encoded)),('bad_dt',lambda:store.begin(q,third.v,-1,flow)),('fractional_tail',lambda:store.begin(q,third.v,.1,flow,max_intervals=1.2))]:
            try:call()
            except ValueError:rejected.append(label)
            else:raise AssertionError('accepted '+label)
            assert store.state.fingerprint()==before
        store.accept(growing);accepted=store.state.fingerprint()
        try:store.accept(growing)
        except ValueError:rejected.append('stale_acceptance')
        else:raise AssertionError('accepted stale trial')
        assert store.state.fingerprint()==accepted;rows.append(dict(case='trial_atomic_failure_and_epoch',negative_controls=rejected,state_unchanged_on_failure=True))
        # Exercise whole-facet intersections, including shared entities whose
        # contact extends beyond the permitted vertex/edge.
        a=np.array([[0.,0,0],[1.,0,0],[0,1.,0]]);base=('a',('A','B','C'),a);checks=[]
        allowed=[('shared_vertex',('A','D','E'),np.array([[0.,0,0],[-1.,0,0],[0,-1.,0]])),('shared_edge',('B','A','D'),np.array([[1.,0,0],[0,0,0],[0,-1.,0]]))]
        for name,ids,x in allowed:guard.validate_facets([base,(name,ids,x)]);checks.append(name)
        forbidden=[('coplanar_overlap',('D','E','F'),a*.4+np.array([.1,.1,0.])),('crossing',('D','E','F'),np.array([[.2,.2,-.3],[.2,.2,.3],[.8,.2,0.]])),('beyond_vertex',('A','D','E'),a*.4),('beyond_edge',('B','A','D'),np.array([[1.,0,0],[0,0,0],[.2,.2,0]]))]
        for name,ids,x in forbidden:
            try:guard.validate_facets([base,(name,ids,x)])
            except ValueError:checks.append(name)
            else:raise AssertionError('missed global geometry '+name)
        rows.append(dict(case='whole_facet_geometry_guards',cases=checks));print('ALL_GROUPS_OK',len(rows),flush=True)
        if '--publish' in sys.argv:
            dest=ROOT/'build/analysis/uvlm-wake-assembly-snapshot.npz';assert not dest.exists();np.savez_compressed(dest,A=A,W=W,A_independent=Aref,gamma=trial['gamma'],q=q,force=trial['loads']['F'],Q=trial['loads']['Q'],integer_D1=model.D1,integer_D2=model.D2,JVP_numeric=numeric,JVP_analytic=analytic);artifact=dict(artifact=dest.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    print(json.dumps(rows,indent=2),flush=True)
    if '--publish' in sys.argv:
        dest=ROOT/'docs/verification/uvlm-wake-assembly-result.json';assert not dest.exists();paths=[Path(__file__),Path(aero.__file__),Path(binding.__file__),Path(guard.__file__),Path(reference.__file__),Path(old.__file__),Path(pg.__file__),ROOT/'docs/verification/uvlm-assembly-adoption-protocol.md',ROOT/'docs/verification/uvlm-induction-adoption-audit.json'];dest.write_text(json.dumps(dict(passed=len(rows)==12 and not failures,records=rows,failures=failures,artifact=artifact,seconds=time.monotonic()-start,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print('UVLM_WAKE_ASSEMBLY_RESULT',len(rows)==12 and not failures,flush=True)
if __name__=='__main__':main()
