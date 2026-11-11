"""Reproduce the separately frozen band-aligned fixtures; --check never writes them."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib
import json
import math
import sys

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    check='--check' in sys.argv
    old_path=FOLDER/'finite-wing-reference-inputs.json';old=json.loads(old_path.read_text())
    out=FOLDER/'fixtures/avl-fw-aligned';meshes=sorted(set([(m,160) for m in [3,6,12,24]]+[(24,n) for n in [20,40,80,160]]+list(zip([3,6,12,24],[20,40,80,160]))))
    bounds=[-4+.4*j for j in range(21)];cases=[];checked=0
    for original in old['cases']:
        case={k:v for k,v in original.items() if k not in ['decks','sections']};sections=[]
        for y in bounds:
            sections.append(dict(x_le=abs(y)*math.tan(math.radians(35)),y=y,z=abs(y)*math.tan(math.radians(original['dihedral_deg'])),chord=1.5-.9*abs(y)/4,incidence_deg=0.))
        area=sum((b['y']-a['y'])*(a['chord']+b['chord'])/2 for a,b in zip(sections[:-1],sections[1:]));assert abs(area-8.4)<1e-13
        decks=[]
        for m,n in meshes:
            text=f'{case["id"]} thin swept tapered wing\n0.0\n0 0 0.0\n8.4 1.1142857142857143 8.0\n0.375 0.0 0.0\n0.0\nSURFACE\nMainWing\n{m} 1.0\nCOMPONENT\n1\n'
            for s in sections:text+='SECTION\n'+' '.join(format(s[k],'.17g') for k in ['x_le','y','z','chord','incidence_deg'])+f' {n//20} 1.0\n'
            data=text.replace('\n','\r\n').encode('ascii');p=out/f'{case["id"]}-M{m}-N{n}.avl'
            if check:assert p.read_bytes()==data,p
            else:
                if p.exists():raise FileExistsError('preserve existing fixtures; use --check')
                p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
            decks.append(dict(m=m,n=n,path=p.relative_to(FOLDER).as_posix(),sha256=hashlib.sha256(data).hexdigest()));checked+=1
        case.update(sections=sections,decks=decks,computed_area_xy=area);cases.append(case)
    result=dict(kind='separate_band_aligned_design_fixture',protocol='avl-band-aligned-grid-experiment.md',parent_manifest_sha256=digest(old_path),bounds=bounds,cases=cases)
    data=(json.dumps(result,indent=2)+'\n').replace('\n','\r\n').encode()
    path=FOLDER/'avl-band-aligned-inputs.json'
    if check:
        assert path.read_bytes()==data
        audit=dict(source_sha256=digest(Path(__file__)),fixture_sha256=digest(path),checked_decks=checked,byte_identical=True)
        (FOLDER/'avl-band-aligned-generation-check.json').write_text(json.dumps(audit,indent=2)+'\n')
    else:
        if path.exists():raise FileExistsError('preserve existing fixture manifest')
        path.write_bytes(data)
    print(json.dumps(dict(checked_decks=checked,byte_identical=check)))

if __name__=='__main__':main()
