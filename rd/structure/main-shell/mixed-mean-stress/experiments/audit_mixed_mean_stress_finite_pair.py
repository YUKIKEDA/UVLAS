"""Audit frozen sources, saved fields, scope flags and current documentation links."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import ast,hashlib,json,re
import numpy as np
ROOT=rd_paths.ROOT
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    start=ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-pair-disposition.json'
    checked={};todo=[start];artifacts={}
    while todo:
        path=todo.pop()
        if path in checked:continue
        assert path.exists(),path
        checked[path]=digest(path)
        if path.suffix!='.json':continue
        data=json.loads(path.read_text(encoding='utf-8'))
        for name,expected in data.get('sources',{}).items():
            source=ROOT/name;assert digest(source)==expected,('source_hash',name)
            if source.suffix=='.json':todo.append(source)
            else:checked[source]=expected
        for record in data.get('records',[]):
            if 'artifact' in record:
                p=ROOT/record['artifact'];assert digest(p)==record['artifact_sha256'],p;artifacts[p]=digest(p)
    disposition=json.loads(start.read_text(encoding='utf-8'))
    assert disposition['initial_package_completed'] and not disposition['main_shell_adopted'] and not disposition['physical_PDE_acceptance']
    assert disposition['overall_completion_percentage'] is None
    read=lambda name:json.loads((rd_paths.RECORDS/name).read_text(encoding='utf-8'))
    primary=read('mixed-mean-stress-finite-pair-result.json');full=read('mixed-mean-stress-finite-full-v2-result.json');old=read('mixed-mean-stress-finite-full-result.json')
    assert primary['all_checks_pass'] and full['all_checks_pass'] and not old['all_checks_pass']
    differences=[]
    for a,b in zip(primary['records'],full['records']):
        assert a['h']==b['h']
        with np.load(ROOT/a['artifact']) as p,np.load(ROOT/b['artifact']) as q:
            assert p['local'].shape==q['local'].shape==(4,558,3)
            difference=np.linalg.norm(p['local']-q['local'],axis=(0,1))/np.linalg.norm(p['local'],axis=(0,1))
            assert np.allclose(difference,b['field_difference'],rtol=1e-13,atol=1e-25)
            assert max(difference)<=1e-8 and np.isfinite(q['trace']).all()
            differences.append(difference.tolist())
    parsed=[]
    for name in ['moment_position_recovery_action','check_moment_position_recovery_action','mixed_mean_stress_finite','check_mixed_mean_stress_finite_pair','check_mixed_mean_stress_finite_full','check_mixed_mean_stress_finite_full_v2','audit_mixed_mean_stress_finite_pair']:
        p=rd_paths.EXPERIMENTS/(name+'.py');ast.parse(p.read_text(encoding='utf-8'));parsed.append(p.relative_to(ROOT).as_posix())
    docs=['README.md','rd/planning/preimplementation-plan.md','rd/planning/structural-model-adoption-review.md','rd/planning/requirements-evidence-map.md','docs/specs/design-decisions.md','rd/planning/preimplementation-history.md','rd/structure/main-shell/review/main-shell-completion-outlook-v9.md','rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-pair-results.md','rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-full-v2-protocol.md','rd/structure/main-shell/mixed-mean-stress/mixed-mean-stress-finite-formulation.md']
    link_count=0
    for name in docs:
        p=ROOT/name;text=p.read_text(encoding='utf-8')
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if '://' in target or target.startswith('#'):continue
            target=target.split('#')[0];assert (p.parent/target).exists(),(name,target);link_count+=1
    dest=ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-pair-audit.json';assert not dest.exists()
    result=dict(all_checks_pass=True,checked_source_files=len(checked),checked_artifacts=len(artifacts),parsed_python=parsed,checked_document_links=link_count,recomputed_field_differences=differences,main_shell_adopted=False,physical_PDE_acceptance=False,sources={start.relative_to(ROOT).as_posix():digest(start),Path(__file__).relative_to(ROOT).as_posix():digest(Path(__file__))})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
