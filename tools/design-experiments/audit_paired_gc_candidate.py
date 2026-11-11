"""Closure audit: candidate gates, preserved source ancestry, artifacts and entry links."""
from pathlib import Path
import ast,hashlib,json,os,platform,re,sys
from urllib.parse import unquote
import numpy,scipy
ROOT=Path(__file__).resolve().parents[2]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    dest=ROOT/'docs/verification/paired-geometry-ancf-candidate-audit.json';assert not dest.exists()
    disposition=ROOT/'docs/verification/paired-geometry-ancf-candidate-disposition.json'
    queue=[disposition];seen=set();sources={};artifacts={};errors=[]
    while queue:
        path=queue.pop()
        if path in seen:continue
        seen.add(path)
        if not path.exists():errors.append('missing:'+str(path));continue
        if path.suffix!='.json':continue
        data=json.loads(path.read_text(encoding='utf-8-sig'))
        for name,expected in (data.get('sources',{}) if isinstance(data,dict) else {}).items():
            if not isinstance(expected,str) or not re.fullmatch('[0-9a-f]{64}',expected):continue
            p=ROOT/name
            if not p.exists() or sha(p)!=expected:errors.append('source_hash:'+name)
            if name in sources and sources[name]!=expected:errors.append('conflicting_hash:'+name)
            sources[name]=expected;queue.append(p)
        def walk(value):
            if isinstance(value,dict):
                if 'artifact' in value and 'artifact_sha256' in value:
                    name=value['artifact'];expected=value['artifact_sha256'];p=ROOT/name;artifacts[name]=expected
                    if not p.exists() or sha(p)!=expected:errors.append('artifact_hash:'+name)
                for v in value.values():walk(v)
            elif isinstance(value,list):
                for v in value:walk(v)
        walk(data)
    parsed=0
    for name in sources:
        if name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf-8-sig'),filename=name);parsed+=1
    documents=['README.md','docs/roadmap.md','docs/preimplementation-plan.md','docs/structural-model-adoption-review.md','docs/requirements-evidence-map.md','docs/design-decisions.md','docs/paired-geometry-ancf-definition.md','docs/paired-geometry-ancf-invariants.md','docs/verification/paired-geometry-ancf-candidate-results.md']
    links=0
    for name in documents:
        p=ROOT/name
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',p.read_text(encoding='utf-8-sig')):
            target=target.strip('<>').split('#')[0]
            if not target or re.match(r'[a-zA-Z]+://',target):continue
            resolved=(p.parent/unquote(target)).resolve();links+=1
            if resolved!=dest and not resolved.exists():errors.append('link:'+name+':'+target)
    result=json.loads((ROOT/'docs/verification/paired-geometry-ancf-candidate-result.json').read_text())
    basis=json.loads((ROOT/'docs/verification/paired-geometry-ancf-basis-audit.json').read_text())
    spatial=json.loads((ROOT/'docs/verification/paired-geometry-ancf-spatial-result.json').read_text())
    gates=result['adoption_candidate_qualified'] and result['all_numerical_checks'] and len(result['records'])==12 and basis['passed'] and spatial['passed'] and len(spatial['records'])==4
    if not gates:errors.append('candidate_gates')
    expected={('standard',amp,n,h) for amp in [0.,.2] for n in [2,4] for h in [.01,.0001]}|{('rational',amp,4,h) for amp in [0.,.2] for h in [.01,.0001]}
    lookup={(a['mode'],a['amplitude'],a['n'],a['h']):a for a in result['records']}
    if set(lookup)!=expected:errors.append('primary_case_coverage')
    for a in result['records']:
        if not all(a['checks'].values()) or a['fit']['rcond']<1e-10 or max(a['fit']['linear'])>1e-8 or a['mass_scaled_min']<=1e-10 or max(a['primary']['original_kkt'])>1e-8 or max(a['traces'].values())>1e-9 or max(a['anchor_consistency'])>1e-8 or a['seconds']>300:errors.append('primary_scalar_gates')
        if a['mode']=='standard' and a['h']==.01 and a['fit']['free_nullities']!=[6]*3:errors.append('free_body_nullity')
        if a['mode']=='standard' and a['n']==2:
            if a['h']==.01 and (not a['mechanics'] or max(a['mechanics'].values())>1e-9):errors.append('mechanics_presence_or_tolerance')
            if a['h']==.0001 and (not a['svd'] or a['svd']['rank']!=a['svd']['columns']):errors.append('independent_svd_presence')
    if result['seconds']>1800 or result['peak_bytes']>2*1024**3 or spatial['seconds']>300 or spatial['peak_bytes']>2*1024**3:errors.append('series_resources')
    for a in spatial['records']:
        if not all(a['checks'].values()) or a['rcond']<1e-10 or a['seconds']>300:errors.append('supplement_scalar_gates')
        lookup['rational',a['amplitude'],2,a['h']]=a
    if len(lookup)!=16:errors.append('combined_case_coverage')
    for mode in ['standard','rational']:
        for amp in [0.,.2]:
            for field in ['error','position']:
                for target in range(3):
                    for n in [2,4]:
                        thick=lookup[mode,amp,n,.01]['independent'][field][target];thin=lookup[mode,amp,n,.0001]['independent'][field][target]
                        if thin>2*thick or (n==4 and max(thick,thin)>.01):errors.append('recomputed_accuracy_or_uniformity')
                    for h in [.01,.0001]:
                        coarse=lookup[mode,amp,2,h]['independent'][field][target];fine=lookup[mode,amp,4,h]['independent'][field][target]
                        if max(coarse,fine)>=1e-12 and fine>=coarse:errors.append('recomputed_spatial_decrease')
    for a in list(result['records'])+list(spatial['records']):
        with numpy.load(ROOT/a['artifact']) as data:
            if not all(numpy.isfinite(data[k]).all() for k in data.files):errors.append('nonfinite_saved_state')
            if data['coef'].shape!=(3*data['Zs'].shape[0],3) or data['beta'].shape!=(16*a['n'],3):errors.append('saved_state_dimensions')
    output=dict(audit_passed=not errors,adoption_candidate_qualified=bool(gates),main_shell_adopted=False,checked_source_files=len(sources),checked_artifacts=len(artifacts),parsed_python_files=parsed,checked_document_links=links,errors=errors,environment=dict(python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,platform=platform.platform(),openblas_threads=os.environ.get('OPENBLAS_NUM_THREADS')),sources={disposition.relative_to(ROOT).as_posix():sha(disposition),Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))})
    dest.write_text(json.dumps(output,indent=2),encoding='utf-8');print(json.dumps(output,indent=2));assert not errors

if __name__=='__main__':main()
