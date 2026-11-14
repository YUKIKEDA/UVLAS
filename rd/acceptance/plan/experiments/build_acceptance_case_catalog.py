import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,re
ROOT=rd_paths.ROOT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    dest=ROOT/'rd/acceptance/plan/records/practical-acceptance-case-catalog.json';assert not dest.exists()
    entries=[
      ('V-PG',['GEO-01','STR-01','MAT-01'],['paired-shell-r5-domain-v2']),
      ('V-MAT',['MAT-01','MAT-02','OUT-02'],['material-layer-kernel','material-shell-integration-v2','structural-layer-damping']),
      ('V-BEAM',['STR-03','STR-04','OUT-01'],['beam-se3-physical-v4','beam-section-interface']),
      ('V-JOINT',['STR-02','STR-03','STR-04','BC-01'],['connection-beam-shell','structural-combined']),
      ('V-LOAD',['BC-01','LOAD-01','IO-02'],['uvlm-known-input','structural-history-state']),
      ('V-RECOVER',['OUT-01','OUT-02'],['structural-recovery']),
      ('V-AERO',['PHY-01','GEO-01'],['high-frequency-final','drag-common-limit']),
      ('V-MULTI',['PHY-02','QUAL-05'],['multiple-surface-harmonic-refined','uvlm-global-binding']),
      ('V-J1',['SOL-01','SOL-02','SOL-03','QUAL-05','QUAL-06'],['pg-uvlm-free-step','coupled-bdf2-prescribed']),
      ('V-SEARCH',['SOL-03','RUN-02','QUAL-06'],['static-branch-continuation','flutter-search-control']),
      ('V-STATE',['RUN-01','RUN-02','QUAL-03','QUAL-04'],['coupled-checkpoint']),
      ('V-USE',['IO-01','IO-02','IO-03','API-01','API-02','ENV-01','ENV-02','DIST-01','DIST-02'],[]),
      ('P-CURVED-BOX',['GEO-01','STR-01','STR-02','STR-03','MAT-01','SOL-01','SOL-02','SOL-03','PHY-02','OUT-01','OUT-02'],[]),
      ('X-S10-GVT-LE',['STR-01','STR-03','STR-04'],[]),
      ('X-S10-GVT-TE',['STR-01','STR-03','STR-04'],[]),
      ('X-S10-STATIC',['SOL-01','GEO-01'],[]),
      ('X-S10-FLUTTER',['SOL-03','QUAL-06'],[]),
      ('X-S10-DECAY',['SOL-02','LOAD-01'],[]),
      ('PERF-FL-01',['PERF-01','PERF-02','PERF-03'],[])]
    all_ids=set(re.findall(r'^\| ([A-Z]+-\d+) \|', (ROOT/'docs/requirements/product-requirements.md').read_text(encoding='utf-8'),re.M));all_ids={k for k in all_ids if not k.startswith('START-')};cases=[];coverage={k:[] for k in all_ids};source_files={};artifacts={}
    for cid,req,roots in entries:
        req=req+['QUAL-01','QUAL-02','QUAL-03','QUAL-04','QUAL-07'];req=sorted(set(req));assert set(req)<=all_ids
        for k in req:coverage[k].append(cid)
        evidence=[]
        for root in roots:
            p=ROOT/rd_paths.record(f'{root}-result.json');d=json.loads(p.read_text());assert (d['all_checks_pass'] if root=='material-layer-kernel' else d['passed']) is True,root;name=p.relative_to(ROOT).as_posix();source_files[name]=sha(p)
            fixtures={k:h for k,h in d.get('sources',{}).items() if '/experiments/' in k and k.endswith('.py')};assert fixtures,root
            for k,h in fixtures.items():assert sha(ROOT/k)==h,k;source_files[k]=h
            evidence.append(dict(result=name,result_sha256=sha(p),fixture_sources=fixtures,meaning='frozen design/reference witness; product migration must execute the same inputs independently'))
        experimental=cid.startswith('X-');cases.append(dict(id=cid,requirements=req,specification='docs/requirements/practical-acceptance-plan.md',specified=True,reference_fixture_available=bool(evidence),experimental_reference_qualified=False if experimental else None,product_accepted=False,design_evidence=evidence))
    assert all(coverage.values()),[k for k,v in coverage.items() if not v]
    for name in ['docs/requirements/product-requirements.md','docs/requirements/acceptance-criteria.md','docs/requirements/practical-acceptance-plan.md','rd/acceptance/s10/records/s10-reference-metadata-result.json','rd/foundation/numerics/records/numerical-foundation-adoption-closure.json',Path(__file__).relative_to(ROOT).as_posix()]:source_files[name]=sha(ROOT/name)
    r=dict(passed=True,meaning='planning coverage and frozen input/reference availability, not product acceptance',case_count=len(cases),requirement_count=len(all_ids),all_requirements_mapped=True,all_product_acceptance_complete=False,cases=cases,requirement_coverage=coverage,sources=source_files,artifacts=artifacts)
    dest.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,cases=len(cases),requirements=len(all_ids),product_acceptance_complete=False)))
if __name__=='__main__':main()
