from pathlib import Path
import ast,json,hashlib
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    sources={};artifacts={};seen=set();statuses={}
    def add(name,h,kind='sources'):
        p=ROOT/name;assert sha(p)==h,name;d=sources if kind=='sources' else artifacts
        if name in d:assert d[name]==h,name
        d[name]=h;return p
    names=['coupled-added-mass','coupled-bdf2-prescribed','coupled-bdf2-resource','coupled-chart-birth','pg-uvlm-coupled-step','pg-uvlm-free-step','so3-time-history','material-birth-derivative','age-initialization-derivative-initial','age-initialization-derivative-refined','coupled-checkpoint','static-branch-continuation','flutter-search-control'];queue=[ROOT/'docs/verification/uvlm-integrated-adoption-audit.json',ROOT/'docs/verification/uvlm-integrated-adoption-closure.json']+[ROOT/f'docs/verification/{n}-result.json' for n in names]
    while queue:
        p=queue.pop()
        if p in seen:continue
        seen.add(p);add(p.relative_to(ROOT).as_posix(),sha(p));r=json.loads(p.read_text(encoding='utf-8'))
        if 'passed' in r:statuses[p.relative_to(ROOT).as_posix()]=r['passed']
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
        checked.add(p);modules=[]
        for n in ast.walk(ast.parse(p.read_text(encoding='utf-8-sig'))):
            if isinstance(n,ast.Import):modules.extend(a.name.split('.')[0] for a in n.names)
            elif isinstance(n,ast.ImportFrom) and n.module:modules.append(n.module.split('.')[0])
        for m in modules:
            q=p.parent/(m+'.py')
            if q.exists():add(q.relative_to(ROOT).as_posix(),sha(q));pending.append(q)
    for n in names:assert statuses[f'docs/verification/{n}-result.json'] is (False if n.endswith('-initial') else True),n
    docs=['coupled-adopted-specification.md','coupled-initialization-residual-contract.md','coupled-time-stability-specification.md','j1-discrete-equilibrium-and-chart-contract.md','j1-initialization-output-derivative-addendum.md','coupled-state-search-contract.md','j1-static-continuation-error-contract.md','verification/j1-adoption-completion.md']
    for n in docs:
        p=ROOT/'docs'/n;add(p.relative_to(ROOT).as_posix(),sha(p))
    logs=['coupled-added-mass','coupled-bdf2-prescribed','coupled-bdf2-resource','coupled-chart-birth','pg-uvlm-coupled-step','pg-uvlm-free-step','so3-time-history','material-birth-derivative','material-birth-derivative-shape-fixed','age-initialization-derivative','age-initialization-derivative-diagnostic','age-initialization-derivative-refined','coupled-checkpoint','static-branch-continuation','flutter-search-control']
    for n in logs:
        p=ROOT/f'build/analysis/{n}.log';add(p.relative_to(ROOT).as_posix(),sha(p),'artifacts')
    for n in ['check_material_birth_derivative-initial.py','check_age_initialization_derivative-initial.py']:
        p=ROOT/'build/analysis'/n;add(p.relative_to(ROOT).as_posix(),sha(p),'artifacts')
    add(Path(__file__).relative_to(ROOT).as_posix(),sha(Path(__file__)))
    result=dict(passed=True,meaning='J1 formulation and scoped verification adoption integrity, not product acceptance',j1_evidence_complete=True,product_implementation_start=False,selected_results=names,recorded_result_statuses=statuses,sources=sources,artifacts=artifacts);dest=ROOT/'docs/verification/j1-adoption-audit.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,sources=len(sources),artifacts=len(artifacts),selected_evidence=len(names))))
if __name__=='__main__':main()
