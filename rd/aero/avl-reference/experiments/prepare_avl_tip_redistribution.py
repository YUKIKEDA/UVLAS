"""Four fixed-count allocation decks; no AVL execution or precision acceptance."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
ALLOCATIONS={'tip-rich':[48]*8+[128]*2,'tip-poor':[72]*8+[32]*2}
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=F/'avl-tip-redistribution-preparation.json'
    if output.exists():raise FileExistsError(output)
    parent=F/'avl-span1280-fine-ks-runs.json';runs=json.loads(parent.read_text());assert runs['completed']
    folder=F/'fixtures/avl-fw-tip-redistribution';folder.mkdir();decks=[]
    for case in ['FW-V01','FW-V02']:
        donor=next(r for r in runs['runs'] if r['case']==case and r['m']==12);source=ROOT/donor['directory']/'wing.avl';lines=source.read_bytes().splitlines(keepends=True)
        indices=[i+1 for i,s in enumerate(lines) if s.strip()==b'SECTION'];assert len(indices)==11
        for label,allocation in ALLOCATIONS.items():
            assert sum(allocation)==640 and len(allocation)==10
            modified=list(lines)
            for index,count in zip(indices[:10],allocation):
                old=lines[index];assert old.endswith(b' 64 1.0\r\n');modified[index]=old[:-len(b' 64 1.0\r\n')]+f' {count} 1.0\r\n'.encode()
                assert modified[index].split()[:5]==old.split()[:5]
            p=folder/f'{case}-M12-N1280-{label}.avl';p.write_bytes(b''.join(modified))
            assert modified[indices[-1]]==lines[indices[-1]]
            assert [i for i,(a,b) in enumerate(zip(lines,modified)) if a!=b]==indices[:10]
            decks.append(dict(case=case,m=12,physical_n=1280,native_strips=640,native_unknowns=7680,allocation=label,
                interval_counts=allocation,unused_final_count=64,path=p.relative_to(ROOT).as_posix(),sha256=digest(p),
                donor=source.relative_to(ROOT).as_posix(),donor_sha256=digest(source),baseline_directory=donor['directory']+'/out'))
    exe=ROOT/'build/reference-avl/source-symmetry16000-ks1000-spacing-pi64/bin/avl'
    result=dict(source_sha256=digest(__file__),protocol='rd/aero/avl-reference/records/avl-tip-redistribution-experiment.md',protocol_sha256=digest(F/'avl-tip-redistribution-experiment.md'),
        parent_sha256=digest(parent),executable_sha256=digest(exe),decks=decks,criteria_pass=True,new_avl_runs=0,reference_acceptance=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(decks=4,all_native_unknowns=7680,all_native_strips=640,no_runs=True)))
if __name__=='__main__':main()
