"""Four independent-factor decks preserving all original physical section bytes."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=F/'avl-allocation-factorial-preparation.json'
    if output.exists():raise FileExistsError(output)
    audit=F/'avl-symmetric-solve-artifact-audit.json';assert json.loads(audit.read_text())['all_criteria_pass']
    parent=F/'avl-span1280-fine-ks-runs.json';data=json.loads(parent.read_text());assert data['completed']
    folder=F/'fixtures/avl-fw-allocation-factorial';folder.mkdir();decks=[]
    for case in ['FW-V01','FW-V02']:
        donor=next(r for r in data['runs'] if r['case']==case and r['m']==12);source=ROOT/donor['directory']/'wing.avl';lines=source.read_bytes().splitlines(keepends=True)
        positions=[j+1 for j,s in enumerate(lines) if s.strip()==b'SECTION'];assert len(positions)==11
        for label,counts in [('inner-only',[48]*8+[64]*2),('tip-only',[64]*8+[128]*2)]:
            changed=list(lines)
            for j,count in zip(positions[:10],counts):
                old=lines[j];assert old.endswith(b' 64 1.0\r\n');changed[j]=old[:-len(b' 64 1.0\r\n')]+f' {count} 1.0\r\n'.encode();assert changed[j].split()[:5]==old.split()[:5]
            n=2*sum(counts);p=folder/f'{case}-M12-N{n}-{label}.avl';p.write_bytes(b''.join(changed))
            assert changed[positions[-1]]==lines[positions[-1]] and sum(counts) in [512,768]
            decks.append(dict(case=case,allocation=label,m=12,physical_n=n,native_strips=sum(counts),native_unknowns=12*sum(counts),interval_counts=counts,
                path=p.relative_to(ROOT).as_posix(),sha256=digest(p),donor=source.relative_to(ROOT).as_posix(),donor_sha256=digest(source)))
    exe=ROOT/'build/reference-avl/source-symmetry16000-ks1000-spacing-pi64/bin/avl'
    result=dict(source_sha256=digest(__file__),protocol='rd/aero/avl-reference/records/avl-allocation-factorial-experiment.md',protocol_sha256=digest(F/'avl-allocation-factorial-experiment.md'),
        solve_audit_sha256=digest(audit),parent_sha256=digest(parent),tip_parent_sha256=digest(F/'avl-tip-redistribution-runs.json'),executable_sha256=digest(exe),
        decks=decks,new_runs=0,all_criteria_pass=True,reference_acceptance=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(decks=4,unknowns=[6144,9216],new_runs=0)))
if __name__=='__main__':main()
