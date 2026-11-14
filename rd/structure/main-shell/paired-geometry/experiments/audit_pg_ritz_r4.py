"""Audit evidence integrity and rejection logic; audit success is not R4 success."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import ast,hashlib,json,re
import numpy as np
ROOT=rd_paths.ROOT


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-r4-audit.json';assert not dest.exists()
    disposition=ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-r4-disposition.json'
    decision=json.loads(disposition.read_text());errors=[]
    assert decision['accepted'] is False and decision['R4_complete'] is False
    follow=json.loads((ROOT/'rd/structure/main-shell/paired-geometry/records/pg-ritz-r4-followup-result.json').read_text())
    ratio=np.array(follow['records'][1]['position'])/follow['records'][0]['position']
    if not np.all(ratio>2.) or not np.allclose(ratio,follow['position_thin_ratio'],rtol=1e-12,atol=0):errors.append('thin_rejection_not_reproduced')
    thin=next(r for r in follow['records'] if r['case']=='independent_full_trace_lift_mass_root' and r['h']==.0001)
    if thin['rank']!=[153,159,159] or thin['trace']>1e-9:errors.append('mass_rejection_not_reproduced')
    for name in ['pg-ritz-mechanics-v2-result.json','pg-ritz-independent-state-result.json','pg-ritz-full-trace-work-result.json']:
        if not json.loads((rd_paths.RECORDS/name).read_text())['passed']:errors.append('bounded_result:'+name)
    queue=[disposition];seen=set();sources={};artifacts={}
    while queue:
        path=queue.pop()
        if path in seen:continue
        seen.add(path)
        if path.suffix!='.json':continue
        data=json.loads(path.read_text(encoding='utf-8-sig'))
        if isinstance(data,dict):
            for name,expected in data.get('sources',{}).items():
                p=ROOT/name
                if not p.exists() or sha(p)!=expected:errors.append('source_hash:'+name)
                if name in sources and sources[name]!=expected:errors.append('conflicting_source:'+name)
                sources[name]=expected;queue.append(p)
        def walk(value):
            if isinstance(value,dict):
                if isinstance(value.get('artifact'),str) and 'artifact_sha256' in value:
                    name=value['artifact'];p=ROOT/name;artifacts[name]=value['artifact_sha256']
                    if not p.exists() or sha(p)!=value['artifact_sha256']:errors.append('artifact_hash:'+name)
                for item in value.values():walk(item)
            elif isinstance(value,list):
                for item in value:walk(item)
        walk(data)
    for name in sources:
        if name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf-8-sig'),filename=name)
    documents=['README.md','rd/planning/preimplementation-history.md','rd/planning/preimplementation-plan.md','rd/planning/structural-model-adoption-review.md','rd/planning/requirements-evidence-map.md','docs/specs/design-decisions.md','rd/structure/main-shell/paired-geometry/pg-ritz-general-mechanics-contract.md','rd/structure/main-shell/paired-geometry/records/pg-ritz-r4-results.md']
    links=0
    for name in documents:
        p=ROOT/name
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',p.read_text(encoding='utf-8-sig')):
            target=target.strip('<>').split('#')[0]
            if not target or re.match('[a-zA-Z]+://',target):continue
            q=(p.parent/target).resolve();links+=1
            if q!=dest and not q.exists():errors.append('link:'+name+':'+target)
    result=dict(audit_passed=not errors,meaning='rejection and evidence integrity only',
        accepted=False,R4_complete=False,main_shell_adopted=False,checked_sources=len(sources),
        checked_artifacts=len(artifacts),checked_links=links,position_thin_ratio=ratio.tolist(),
        mass_root_ranks=thin['rank'],errors=errors,
        sources={disposition.relative_to(ROOT).as_posix():sha(disposition),Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print({k:v for k,v in result.items() if k!='sources'});assert not errors


if __name__=='__main__':main()
