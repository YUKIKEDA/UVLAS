"""Verify prospective factor separation and preserve the non-diagnostic reference."""
from pathlib import Path
import hashlib,json,re
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    path=F/'avl-allocation-factorial-preparation.json';p=json.loads(path.read_text())
    assert p['source_sha256']==digest(ROOT/'tools/design-experiments/prepare_avl_allocation_factorial.py')
    assert p['protocol_sha256']==digest(ROOT/p['protocol']) and p['solve_audit_sha256']==digest(F/'avl-symmetric-solve-artifact-audit.json')
    assert p['parent_sha256']==digest(F/'avl-span1280-fine-ks-runs.json') and p['tip_parent_sha256']==digest(F/'avl-tip-redistribution-runs.json')
    base=ROOT/'build/reference-avl/source-symmetry16000-ks1000-spacing-pi64'
    assert p['executable_sha256']==digest(base/'bin/avl') and p['all_criteria_pass'] and p['new_runs']==0 and not p['reference_acceptance']
    dim=(base/'src/ADIMEN.INC').read_text();heap=(base/'src/avl_heap_inc.f90').read_text();make=(base/'src/amake.f').read_text()
    capacity=int(re.search(r'NVMAX=(\d+)',dim)[1]);stripcap=int(re.search(r'NSMAX=(\d+)',dim)[1]);localcap=int(re.search(r'KSMAX=(\d+)',make)[1]);assert int(re.search(r'NVX=(\d+)',heap)[1])==capacity
    assert len(p['decks'])==4 and {(d['case'],d['allocation']) for d in p['decks']}=={(c,a) for c in ['FW-V01','FW-V02'] for a in ['inner-only','tip-only']}
    for d in p['decks']:
        source=ROOT/d['donor'];target=ROOT/d['path'];assert digest(source)==d['donor_sha256'] and digest(target)==d['sha256']
        a=source.read_bytes().splitlines();b=target.read_bytes().splitlines();positions=[j+1 for j,line in enumerate(a) if line==b'SECTION'];assert len(positions)==11
        counts=[48]*8+[64]*2 if d['allocation']=='inner-only' else [64]*8+[128]*2
        assert counts==d['interval_counts'] and sum(counts)==d['native_strips'] and 2*sum(counts)==d['physical_n'] and 12*sum(counts)==d['native_unknowns']
        assert d['native_unknowns']<=capacity and d['native_strips']<=min(stripcap,localcap)
        expected=[]
        for j,count in zip(positions[:10],counts):
            assert a[j].split()[:5]==b[j].split()[:5] and b[j].split()[-2:]==[str(count).encode(),b'1.0']
            if count!=64:expected.append(j)
        assert [j for j,(x,y) in enumerate(zip(a,b)) if x!=y]==expected and a[positions[-1]]==b[positions[-1]]
    report=dict(source_sha256=digest(__file__),preparation_sha256=digest(path),cases=4,unknowns=[6144,9216],within_current_capacity=True,new_runs=0,all_criteria_pass=True)
    (F/'avl-allocation-factorial-preparation-audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
