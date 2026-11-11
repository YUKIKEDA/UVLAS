import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,re
ROOT=rd_paths.ROOT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    sources={};artifacts={};seen=set()
    def add(name,h=None,kind='sources'):
        p=ROOT/name;actual=sha(p);assert h is None or h==actual,name;table=sources if kind=='sources' else artifacts
        if name in table:assert table[name]==actual,name
        table[name]=actual;return p
    roots=['numerical-foundation-adoption-closure','numerical-foundation-adoption-audit','practical-acceptance-case-catalog','s10-reference-metadata-result'];queue=[ROOT/rd_paths.record(f'{n}.json') for n in roots]
    while queue:
        p=queue.pop()
        if p in seen:continue
        seen.add(p);r=json.loads(p.read_text(encoding='utf-8'));add(p.relative_to(ROOT).as_posix())
        for key in ['sources','immutable','status_sync_source']:
            for name,h in r.get(key,{}).items():
                q=add(name,h)
                if q.suffix=='.json':queue.append(q)
        for name,h in r.get('artifacts',{}).items():add(name,h,'artifacts')
        if isinstance(r.get('artifact'),dict):a=r['artifact'];add(a['artifact'],a['artifact_sha256'],'artifacts')
    for name in ['docs/plans/implementation-start-review.md','docs/requirements/practical-acceptance-plan.md','rd/acceptance/plan/records/practical-acceptance-plan-completion.md',rd_paths.experiment(Path(__file__).name)]:add(name)
    for name in ['s10-public-reassessment.log','s10-case-input-reassessment.log','s10-reference-metadata.log']:add('build/analysis/'+name,kind='artifacts')
    catalog=json.loads((ROOT/'rd/acceptance/plan/records/practical-acceptance-case-catalog.json').read_text());assert catalog['passed'] and catalog['all_requirements_mapped'];assert catalog['case_count']==19 and catalog['requirement_count']==37;assert not catalog['all_product_acceptance_complete']
    assert all(not c['product_accepted'] for c in catalog['cases'])
    s10=json.loads((ROOT/'rd/acceptance/s10/records/s10-reference-metadata-result.json').read_text());assert s10['passed'] and not s10['experimental_reference_qualified']
    report=ROOT/'docs/plans/implementation-start-review.md';assert 'IM-00' in report.read_text(encoding='utf-8')
    for name in ['docs/plans/implementation-start-review.md','docs/requirements/practical-acceptance-plan.md','rd/acceptance/plan/records/practical-acceptance-plan-completion.md']:
        p=ROOT/name
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' not in target:assert (p.parent/target.split('#')[0]).resolve().exists(),(name,target)
    r=dict(passed=True,structural_review_complete=True,uvlm_review_complete=True,j1_review_complete=True,numerical_foundation_review_complete=True,practical_acceptance_plan_complete=True,initial_implementation_review_complete=True,implementation_start_allowed_packages=['IM-00'],product_implementation_started=False,all_product_acceptance_complete=False,experimental_s10_reference_qualified=False,product_scale_performance_qualified=False,meaning='readiness of the first dependency-closed product package; not completed solver, validation or release',case_count=19,requirement_count=37,risk_ids=[f'R-{i:02}' for i in range(1,8)],sources=sources,artifacts=artifacts)
    dest=ROOT/'rd/acceptance/plan/records/implementation-readiness-audit.json';assert not dest.exists();dest.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,sources=len(sources),artifacts=len(artifacts),start_allowed=r['implementation_start_allowed_packages'],product_implementation_started=False,product_accepted=False)))
if __name__=='__main__':main()
