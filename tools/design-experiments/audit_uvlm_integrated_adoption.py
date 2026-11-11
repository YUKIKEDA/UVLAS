"""Preserve parent evidence, failed outcomes, and all imported design sources."""
from pathlib import Path
import ast,hashlib,json
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    sources={};artifacts={};seen=set();statuses={}
    def add(name,h,kind='sources'):
        p=ROOT/name;assert sha(p)==h,name;group=sources if kind=='sources' else artifacts
        if name in group:assert group[name]==h,name
        group[name]=h;return p
    stems=['trefftz-midpoint-operator','drag-material-self','drag-common-limit','uvlm-error-certification','age-scalar','frequency-delta-candidate','frequency-delta-reference','uvlm-known-input','uvlm-global-binding','multiple-surface-harmonic','multiple-surface-harmonic-ideal','multiple-surface-harmonic-refined','ptera-scalar-reference','frequency-settings-and-pitch-time','ideal-pitch-frequency','high-frequency-final']
    queue=[ROOT/'docs/verification/uvlm-reference-campaign-audit.json',ROOT/'docs/verification/uvlm-reference-campaign-closure.json']+[ROOT/f'docs/verification/{s}-result.json' for s in stems];queue.extend([ROOT/'docs/verification/high-frequency-final-symmetry-proof.json',ROOT/'docs/verification/high-frequency-symmetry-residual-diagnostic.json'])
    while queue:
        p=queue.pop()
        if p in seen:continue
        seen.add(p);name=p.relative_to(ROOT).as_posix();add(name,sha(p));r=json.loads(p.read_text(encoding='utf-8'))
        if 'passed' in r:statuses[name]=r['passed']
        for group in ['sources','immutable','status_sync_source']:
            for name,h in r.get(group,{}).items():
                q=add(name,h)
                if q.suffix=='.json':queue.append(q)
        for name,h in r.get('artifacts',{}).items():add(name,h,'artifacts')
        if isinstance(r.get('artifact'),dict):a=r['artifact'];add(a['artifact'],a['artifact_sha256'],'artifacts')
    pending=[ROOT/n for n in sources if n.endswith('.py')];checked=set()
    while pending:
        p=pending.pop()
        if p in checked:continue
        checked.add(p);tree=ast.parse(p.read_text(encoding='utf-8-sig'));modules=[]
        for n in ast.walk(tree):
            if isinstance(n,ast.Import):modules.extend(a.name.split('.')[0] for a in n.names)
            elif isinstance(n,ast.ImportFrom) and n.module:modules.append(n.module.split('.')[0])
        for module in modules:
            q=p.parent/(module+'.py')
            if q.exists():add(q.relative_to(ROOT).as_posix(),sha(q));pending.append(q)
    required=['drag-common-limit','uvlm-error-certification','multiple-surface-harmonic-refined','ptera-scalar-reference','frequency-settings-and-pitch-time','high-frequency-final','uvlm-known-input','uvlm-global-binding']
    for s in required:assert statuses[f'docs/verification/{s}-result.json'] is True
    for s in ['multiple-surface-harmonic','multiple-surface-harmonic-ideal','ideal-pitch-frequency']:assert statuses[f'docs/verification/{s}-result.json'] is False
    # Include completed campaign records that are not individually referenced by a
    # publisher, without freezing still-running logs or the six live status pages.
    for folder in ['drag-common-limit-v2','drag-core-resolution','frequency-delta','frequency-ideal-reference','frequency-ideal-reference-q24','frequency-settings','pitch-frequency-2-3','multiple-surface-harmonic','multiple-surface-harmonic-ideal','multiple-surface-harmonic-ideal-q32','multiple-surface-harmonic-ideal-q48','high-frequency-final-qualified']:
        for p in (ROOT/'build/analysis'/folder).glob('*'):
            if p.suffix in ['.npz','.json']:add(p.relative_to(ROOT).as_posix(),sha(p),'artifacts')
    for name in ['uvlm-integrated-specification.md','verification/uvlm-adoption-completion.md']:
        p=ROOT/'docs'/name;add(p.relative_to(ROOT).as_posix(),sha(p))
    add(Path(__file__).relative_to(ROOT).as_posix(),sha(Path(__file__)))
    for group in [sources,artifacts]:
        for name,h in group.items():assert sha(ROOT/name)==h,name
    result=dict(passed=True,meaning='immutable source/artifact and scoped-outcome integrity; method adoption separate from product acceptance',structural_review_complete=True,uvlm_2_1_complete=True,uvlm_2_2_complete=True,uvlm_2_3_evidence_complete=True,uvlm_2_4_evidence_complete=True,product_implementation_start=False,recorded_result_statuses=statuses,sources=sources,artifacts=artifacts);dest=ROOT/'docs/verification/uvlm-integrated-adoption-audit.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,sources=len(sources),artifacts=len(artifacts),recorded_outcomes=len(statuses)),indent=2))
if __name__=='__main__':main()
