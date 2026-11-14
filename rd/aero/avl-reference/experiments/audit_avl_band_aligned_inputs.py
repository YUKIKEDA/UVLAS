import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib
import json
import math

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    path=FOLDER/'avl-band-aligned-inputs.json';fixture=json.loads(path.read_text())
    old_path=FOLDER/'finite-wing-reference-inputs.json';old=json.loads(old_path.read_text())
    assert fixture['parent_manifest_sha256']==digest(old_path)
    expected_meshes=sorted(set([(m,160) for m in [3,6,12,24]]+[(24,n) for n in [20,40,80,160]]+list(zip([3,6,12,24],[20,40,80,160]))))
    for case in fixture['cases']:
        original=next(c for c in old['cases'] if c['id']==case['id'])
        assert [(d['m'],d['n']) for d in case['decks']]==expected_meshes
        assert [s['y'] for s in case['sections']]==fixture['bounds']
        assert len(case['sections'])==21
        for s in case['sections']:
            y=s['y'];eta=abs(y)/4
            assert s['x_le']==abs(y)*math.tan(math.radians(35))
            assert s['chord']==1.5-.9*eta and s['incidence_deg']==0
            assert s['z']==abs(y)*math.tan(math.radians(case['dihedral_deg']))
        for a in original['sections']:
            b=next(s for s in case['sections'] if s['y']==a['y'])
            assert all(a[k]==b[k] for k in ['x_le','chord','z','incidence_deg'])
        area=sum((b['y']-a['y'])*(a['chord']+b['chord'])/2 for a,b in zip(case['sections'][:-1],case['sections'][1:]))
        assert abs(area-8.4)<=1e-13
        for deck in case['decks']:
            p=FOLDER/deck['path'];assert digest(p)==deck['sha256']
            lines=p.read_text().splitlines()
            sections=[list(map(float,lines[i+1].split())) for i,line in enumerate(lines) if line=='SECTION']
            assert len(sections)==21
            for source,encoded in zip(case['sections'],sections):
                assert encoded[:5]==[source[k] for k in ['x_le','y','z','chord','incidence_deg']]
                assert encoded[5:]==[deck['n']//20,1.]
            assert max(deck['m']*deck['n'] for deck in case['decks'])==3840
    result=dict(protocol='rd/aero/avl-reference/records/avl-band-aligned-grid-experiment.md',source_sha256=digest(Path(__file__)),
        fixture_sha256=digest(path),parent_fixture_sha256=digest(old_path),physical_shape_unchanged=True,
        cases=2,decks=20,sections_per_deck=21,max_unknowns=3840)
    (FOLDER/'avl-band-aligned-input-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(physical_shape_unchanged=True,decks=20,max_unknowns=3840)))

if __name__=='__main__':main()
