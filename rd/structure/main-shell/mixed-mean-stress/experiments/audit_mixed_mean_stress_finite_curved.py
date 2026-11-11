"""Audit preserved curved experiments, numerical flags and independent fields."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import ast,hashlib,json,re
import numpy as np
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(name):return json.loads((rd_paths.RECORDS/name).read_text(encoding='utf-8'))
def main():
    dest=ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-curved-audit.json';assert not dest.exists()
    closure=ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-curved-disposition.json';data=json.loads(closure.read_text())
    checked={};artifacts={};todo=[closure]
    while todo:
        path=todo.pop()
        if path in checked:continue
        checked[path]=digest(path)
        if path.suffix!='.json':continue
        d=json.loads(path.read_text(encoding='utf-8'))
        for name,expected in d.get('sources',{}).items():
            p=ROOT/name;assert digest(p)==expected,('source_hash',name)
            if p.suffix=='.json':todo.append(p)
            else:checked[p]=expected
        for record in d.get('records',[]):
            if 'artifact' in record:
                p=ROOT/record['artifact'];assert digest(p)==record['artifact_sha256'],p;artifacts[p]=digest(p)
    primary=read('mixed-mean-stress-finite-curved-result.json');full=read('mixed-mean-stress-curved-full-result.json')
    assert not data['main_shell_adopted'] and not data['physical_PDE_acceptance']
    recomputed=[];lookup={(r['n'],r['amplitude'],r['h']):r for r in primary['records']}
    for record in primary['records']:
        with np.load(ROOT/record['artifact']) as saved:
            assert saved['local'].shape==(record['n']**2,558,3)
            assert np.isfinite(saved['local']).all() and np.isfinite(saved['trace']).all()
    for record in full['records']:
        parent=lookup[(2,record['amplitude'],record['h'])]
        with np.load(ROOT/parent['artifact']) as p,np.load(ROOT/record['artifact']) as q:
            difference=np.linalg.norm(q['local']-p['local'],axis=(0,1))/np.linalg.norm(p['local'],axis=(0,1))
            assert np.allclose(difference,record['field_difference'],rtol=1e-13,atol=1e-25)
            recomputed.append(difference.tolist())
    ratios=[]
    for amp in [0.,.2]:
        if (4,amp,.01) not in lookup or (4,amp,.0001) not in lookup:continue
        for field in ['error','recovered_position_error']:
            for target in [1,2]:
                values=[lookup[(4,amp,h)]['observed'][field][target] for h in [.01,.0001]]
                ratios.append(dict(amplitude=amp,field=field,target=primary['kinds'][target],thin_amplification=values[1]/values[0]))
    assert any(r['thin_amplification']>2 for r in ratios)
    scripts=['check_mixed_mean_stress_finite_curved','check_mixed_mean_stress_curved_full','audit_mixed_mean_stress_finite_curved']
    for name in scripts:ast.parse((rd_paths.EXPERIMENTS/(name+'.py')).read_text(encoding='utf-8'))
    links=0
    for name in ['README.md','rd/planning/preimplementation-plan.md','rd/planning/structural-model-adoption-review.md','rd/planning/requirements-evidence-map.md','docs/specs/design-decisions.md','rd/planning/preimplementation-history.md','rd/structure/main-shell/review/main-shell-completion-outlook-v10.md','rd/structure/main-shell/review/main-shell-adoption-exit-plan.md','rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-curved-results.md']:
        p=ROOT/name
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            assert (p.parent/target.split('#')[0]).exists(),(name,target);links+=1
    result=dict(audit_passed=True,solver_candidate_adopted=False,numerical_experiment_all_passed=primary['all_numerical_checks'],independent_n2_all_passed=full['all_checks_pass'],checked_source_files=len(checked),checked_artifacts=len(artifacts),checked_document_links=links,recomputed_field_differences=recomputed,recomputed_finest_thin_amplifications=ratios,sources={closure.relative_to(ROOT).as_posix():digest(closure),Path(__file__).relative_to(ROOT).as_posix():digest(Path(__file__))})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
