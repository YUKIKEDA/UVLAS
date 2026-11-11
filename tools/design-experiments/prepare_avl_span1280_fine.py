"""N1280 decks from the validated capacity16000 N640 inputs, preserving section bytes."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=FOLDER/'avl-span1280-fine-preparation.json'
    if output.exists():raise FileExistsError(output)
    audit=FOLDER/'avl-span1280-capacity-artifact-audit.json';assert json.loads(audit.read_text())['all_criteria_pass']
    parent=FOLDER/'avl-span1280-capacity-runs.json';data=json.loads(parent.read_text());assert data['completed']
    folder=FOLDER/'fixtures/avl-fw-span1280';folder.mkdir();decks=[]
    for case in ['FW-V01','FW-V02']:
        donor=next(r for r in data['runs'] if r['case']==case and r['physical_n']==640)
        source=ROOT/donor['directory']/'wing.avl';raw=source.read_bytes()
        assert raw.count(b'SECTION')==11 and raw.count(b' 32 1.0\r\n')==11 and raw.count(b'24 1.0\r\n')==1
        for m in [3,6,12,24]:
            content=raw.replace(b' 32 1.0\r\n',b' 64 1.0\r\n').replace(b'24 1.0\r\n',f'{m} 1.0\r\n'.encode())
            p=folder/f'{case}-M{m}-N1280-half.avl';p.write_bytes(content)
            decks.append(dict(case=case,m=m,physical_n=1280,native_strips=640,native_unknowns=m*640,
                path=p.relative_to(ROOT).as_posix(),sha256=digest(p),donor=source.relative_to(ROOT).as_posix(),donor_sha256=digest(source)))
    exe=ROOT/'build/reference-avl/source-symmetry16000-spacing-pi64/bin/avl'
    result=dict(protocol='docs/verification/avl-span1280-fine-experiment.md',protocol_sha256=digest(FOLDER/'avl-span1280-fine-experiment.md'),
        source_sha256=digest(__file__),capacity_audit_sha256=digest(audit),parent_runs_sha256=digest(parent),
        capacity_preparation_sha256=digest(FOLDER/'avl-span1280-capacity-preparation.json'),executable_sha256=digest(exe),decks=decks)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(decks=len(decks),unknowns=[1920,3840,7680,15360])))
if __name__=='__main__':main()
