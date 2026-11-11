"""Reproducible fine-stage decks, using byte-identical positive-side sections."""
from pathlib import Path
import hashlib,json

ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parent_path=FOLDER/'avl-symmetry-pi-preparation.json';parent=json.loads(parent_path.read_text());decks=list(parent['decks'])
    old=json.loads((FOLDER/'avl-band-aligned-grid-runs.json').read_text())['runs']
    for case in ['FW-V01','FW-V02']:
        for n in [160,640]:
            donor=next(d for d in parent['decks'] if d['case']==case and d['m']==24 and d['physical_n']==n)
            text=(ROOT/donor['path']).read_bytes();assert text.count(b'24 1.0\r\n')==1
            for m in [3,6,12]:
                p=FOLDER/'fixtures/avl-fw-symmetric-pi'/f'{case}-M{m}-N{n}-half.avl';content=text.replace(b'24 1.0\r\n',f'{m} 1.0\r\n'.encode())
                if p.exists():assert p.read_bytes()==content
                else:p.write_bytes(content)
                rec=dict(donor,m=m,path=p.relative_to(ROOT).as_posix(),sha256=digest(p),fine_role='chord_anchor')
                if n==160:
                    r=next(r for r in old if (r['case'],r['m'],r['n'],r['alpha'])==(case,m,n,1))
                    rec['comparison_directory']=r['directory'];rec['full_input']=(ROOT/r['directory']).parent.joinpath('wing.avl').relative_to(ROOT).as_posix()
                decks.append(rec)
    result=dict(parent,decks=decks,source_sha256=digest(Path(__file__)),parent_preparation_sha256=digest(parent_path),
        protocol='docs/verification/avl-symmetry-fine-observation-addendum.md')
    (FOLDER/'avl-symmetry-pi-fine-preparation.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
