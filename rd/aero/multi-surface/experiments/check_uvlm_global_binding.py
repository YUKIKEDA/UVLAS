import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import copy,json,hashlib,numpy as np
from uvlm_global_binding_witness import compile_bindings
ROOT=rd_paths.ROOT
def main():
    left=dict(dofs=[0,1,2,3],vertices=[0,1,2,3],N=np.eye(4),faces=[(0,1,3,2)]);right=dict(dofs=[1,3,4,5],vertices=[1,3,4,5],N=np.eye(4),faces=[(0,2,3,1)]);patches=[left,right];model=compile_bindings(patches,[(2,3),(3,5)]);assert len(model['vertex_ids'])==6 and np.all(model['D1']@model['D2']==0);assert np.array_equal(model['T'],np.eye(2,dtype=int));interior=np.count_nonzero(model['D2'],axis=1)==2;assert np.all(model['D2'][interior]@np.ones(2)==0);records=[dict(case='global_shared_edge_and_oriented_TE')]
    qtable={0:[0,-1,0],1:[0,0,0],2:[1,-1,0],3:[1,0,0],4:[0,1,0],5:[1,1,0]};q=np.array([qtable[k] for k in model['keys']]);N=model['N'];X=N@q;rng=np.random.default_rng(511);F=rng.normal(size=X.shape);dq=rng.normal(size=q.shape);Q=N.T@F;work=abs(np.sum(dq*Q)-np.sum((N@dq)*F));assert work<1e-12;translation=rng.normal(size=3);omega=rng.normal(size=3);rigid=np.cross(omega,q)+translation;dual=float(np.sum(rigid*Q));wrench=float(translation@F.sum(axis=0)+omega@np.cross(X,F).sum(axis=0));assert abs(dual-wrench)<1e-12;records.append(dict(case='shared_position_field_work_and_rigid_wrench',work_difference=float(work),rigid_work_difference=abs(dual-wrench)))
    release={(key,7):j for j,key in enumerate(dict.fromkeys([2,3,3,5]))};assert len(release)==3 and list(release).count((3,7))==1;records.append(dict(case='shared_TE_vertex_single_release_owner'))
    bad=copy.deepcopy(patches);bad[1]['N'][0,0]=.9;bad[1]['N'][0,1]=.1
    for label,ps,te in [('different_position_field',bad,[(2,3),(3,5)]),('internal_TE',patches,[(1,3)]),('duplicate_face',patches+[left],[(2,3)])]:
        try:compile_bindings(ps,te)
        except ValueError:pass
        else:raise AssertionError(label)
    records.append(dict(case='alias_and_topology_negative_controls',rejected=3));paths=[Path(__file__),ROOT/'rd/aero/multi-surface/experiments/uvlm_global_binding_witness.py',ROOT/'docs/specs/aero/uvlm-integrated-specification.md'];result=dict(passed=True,groups=4,records=records,scope='global shared material binding/topology/dual map; full nonlinear fluid assembly remains the previously validated operator composition',sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'rd/aero/multi-surface/records/uvlm-global-binding-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
