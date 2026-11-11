"""Integrity audit, explicitly distinct from all-physics acceptance."""
from pathlib import Path
import ast,hashlib,json,re
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    dest=ROOT/'docs/verification/uvlm-reference-campaign-audit.json';assert not dest.exists();sources={};artifacts={};seen=set();statuses={}
    def add(name,h,kind='sources'):
        p=ROOT/name;assert sha(p)==h,name;group=sources if kind=='sources' else artifacts
        if name in group:assert group[name]==h,name
        group[name]=h
        return p
    queue=[ROOT/'docs/verification/uvlm-wake-adoption-audit.json',ROOT/'docs/verification/uvlm-wake-adoption-closure.json']
    stems=['finite-wing-specialization','finite-wing-initial-campaign','finite-wing-extended-campaign','fixed-band-final-extension','fixed-band-core-plateau','harmonic-age-integral','frequency-acceleration','unsteady-initial-campaign','unsteady-continuum-estimate','symmetric-heave-final','multiple-surface-independent','finite-wing-drag-diagnostic','edge-load-integral']
    queue += [ROOT/'docs/verification'/f'{s}-result.json' for s in stems]
    while queue:
        p=queue.pop()
        if p in seen:continue
        seen.add(p);name=p.relative_to(ROOT).as_posix();add(name,sha(p));r=json.loads(p.read_text(encoding='utf-8'))
        if 'passed' in r:statuses[name]=r['passed']
        for key in ['sources','immutable','status_sync_source']:
            for name,h in r.get(key,{}).items():
                q=add(name,h)
                if q.suffix=='.json':queue.append(q)
        for name,h in r.get('artifacts',{}).items():add(name,h,'artifacts')
        if isinstance(r.get('artifact'),dict):a=r['artifact'];add(a['artifact'],a['artifact_sha256'],'artifacts')
    pending=[ROOT/n for n in sources if n.endswith('.py')];checked=set()
    while pending:
        p=pending.pop()
        if p in checked:continue
        checked.add(p)
        tree=ast.parse(p.read_text(encoding='utf-8-sig'))
        modules=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):modules.extend(a.name.split('.')[0] for a in node.names)
            elif isinstance(node,ast.ImportFrom) and node.module:modules.append(node.module.split('.')[0])
        for module in modules:
            q=p.parent/(module+'.py')
            if q.exists():add(q.relative_to(ROOT).as_posix(),sha(q));pending.append(q)
    required=['finite-wing-extended-campaign','fixed-band-core-plateau','unsteady-continuum-estimate','symmetric-heave-final','multiple-surface-independent']
    for s in required:assert statuses[f'docs/verification/{s}-result.json']
    for s in ['finite-wing-initial-campaign','fixed-band-final-extension','unsteady-initial-campaign','finite-wing-drag-diagnostic']:assert statuses[f'docs/verification/{s}-result.json'] is False
    report=ROOT/'docs/verification/uvlm-reference-campaign-summary-2026-11-05.md';sources[report.relative_to(ROOT).as_posix()]=sha(report);sources[Path(__file__).relative_to(ROOT).as_posix()]=sha(Path(__file__));links=0
    for target in re.findall(r'\]\(([^)]+)\)',report.read_text(encoding='utf-8')):
        if '://' in target:continue
        q=(report.parent/target.split('#')[0]).resolve();assert q.exists() or q==dest,target;links+=1
    for name,h in sources.items():assert sha(ROOT/name)==h,name
    for name,h in artifacts.items():assert sha(ROOT/name)==h,name
    result=dict(passed=True,meaning='integrity and selected scoped claims; not complete UVLM adoption',structural_review_complete=True,uvlm_2_1_complete=True,uvlm_2_2_complete=True,uvlm_2_3_complete=False,uvlm_review_complete=False,product_implementation_start=False,current_work='2.3 induced-drag common limit, unsteady pitch/frequency extension, general setting/error rules',selected_scoped_claims=required,recorded_result_statuses=statuses,sources=sources,artifacts=artifacts,checked_report_links=links);dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,sources=len(sources),artifacts=len(artifacts),scoped_claims=required,uvlm_2_3_complete=False),indent=2))
if __name__=='__main__':main()
