import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,re
ROOT=rd_paths.ROOT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    sources={};artifacts={};seen=set()
    def add(name,h=None,kind='sources'):
        p=ROOT/name;actual=sha(p);assert h is None or actual==h,name
        table=sources if kind=='sources' else artifacts
        if name in table:assert table[name]==actual,name
        table[name]=actual;return p
    queue=[ROOT/'rd/coupled/j1/records/j1-adoption-audit.json',ROOT/'rd/coupled/j1/records/j1-adoption-closure.json']
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
    for folder in ['cpp-numerics','cpp-io-binding']:
        for p in (rd_paths.EXPERIMENTS/folder).iterdir():
            if p.is_file():add(p.relative_to(ROOT).as_posix())
    for name in ['fetch_numerical_dependencies.py','fetch_io_dependencies.py','export_j1_numerical_oracles.py','run_numerical_capacity.py',Path(__file__).name]:add(rd_paths.experiment(name))
    for name in ['numerical-foundation-adoption-specification.md','io-execution-adoption-specification.md','verification/numerical-foundation-adoption-completion.md']:
        p=add(rd_paths.doc(name))
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' not in link:assert (p.parent/link.split('#')[0]).resolve().exists(),link
    tests={}
    for os in ['windows','linux']:
        for family,count in [('numeric-contracts',4),('io-binding',1)]:
            name=f'build/analysis/{family}-{os}-tests.log';p=add(name,kind='artifacts');s=p.read_text(encoding='utf-8');assert f'100% tests passed, 0 tests failed out of {count}' in s,name;tests[f'{family}-{os}']=count
        p=add(f'build/analysis/numerical-capacity-{os}.json',kind='artifacts');r=json.loads(p.read_text());assert r['passed'];assert len(r['cases'])==2
        exe=ROOT/('build/numeric-spike-windows/Release/capacity.exe' if os=='windows' else 'build/numeric-spike-linux/capacity');add(exe.relative_to(ROOT).as_posix(),r['executable_sha256'],'artifacts');add(f'build/analysis/numerical-capacity-{os}.log',kind='artifacts')
        for name in [f'numeric-contracts-{os}-build.log',f'io-binding-{os}-configure.log',f'io-binding-{os}-build.log',f'numeric-capacity-{os}-build.log']:add('build/analysis/'+name,kind='artifacts')
    add('build/analysis/numeric-spike-linux-build.log',kind='artifacts')
    for name in ['pg-free-step.txt','pg-free-oracle.npz','oracle.json']:add('build/analysis/numeric-oracles/'+name,kind='artifacts')
    oracle=json.loads((ROOT/'build/analysis/numeric-oracles/oracle.json').read_text());assert sha(ROOT/'build/analysis/numeric-oracles/pg-free-step.txt')==oracle['sha256'];assert oracle['dimension']==24
    dependencies=[]
    for name in ['header-sources.lock.json','io-sources.lock.json']:
        p=add('build/numerical-deps/'+name,kind='artifacts');lock=json.loads(p.read_text())
        for package in lock['packages']:add(package['archive'],package['sha256'],'artifacts');dependencies.append({k:package[k] for k in ['name','version','sha256','url']})
    result=dict(passed=True,numerical_foundation_review_complete=True,product_implementation_start=False,practical_acceptance_complete=False,current_work='5 practical cases, public references, acceptance/performance plan',scope='adopted algorithms, concrete representations and contracts; scoped C++23/CPython spikes, not implemented product',tests=tests,dependencies=dependencies,sources=sources,artifacts=artifacts)
    dest=ROOT/'rd/foundation/numerics/records/numerical-foundation-adoption-audit.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,sources=len(sources),artifacts=len(artifacts),test_groups=sum(tests.values()),current_work=result['current_work'])))
if __name__=='__main__':main()
