"""Prepare design fixtures and reference provenance; does not run either solver."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib
import json
import math
import tarfile

ROOT = rd_paths.ROOT
CACHE = ROOT/'build/reference-avl'
OUT = rd_paths.RECORDS
EXPECTED = '0b588ecea9222f5b625d0af0c87ae31daf3cdba1532cf0bbb36f93d6e854849b'

def sha(data): return hashlib.sha256(data).hexdigest()

def main():
    archive = CACHE/'avl3.52.tgz'
    assert sha(archive.read_bytes()) == EXPECTED
    selected = ['version_notes.txt','src/avl.f','src/aero.f','src/aic.f',
        'src/aoper.f','src/aoutput.f','avl_doc-AVL350-04242025.txt']
    provenance = []
    with tarfile.open(archive) as bundle:
        for suffix in selected:
            names = [n for n in bundle.getnames() if n == 'AVL3.52rel09032025/'+suffix]
            assert len(names) == 1
            data = bundle.extractfile(names[0]).read()
            assert (CACHE/Path(suffix).name).read_bytes() == data
            provenance.append(dict(archive_member=names[0],sha256=sha(data),bytes=len(data)))
    meshes = sorted(set([(m,96) for m in [4,8,16,32]]+[(32,n) for n in [12,24,48,96]]+list(zip([4,8,16,32],[12,24,48,96]))))
    fixtures = OUT/'fixtures/avl-fw'
    fixtures.mkdir(parents=True,exist_ok=True)
    cases = []
    for ident,dihedral in [('FW-V01',0.),('FW-V02',5.)]:
        sections = [dict(x_le=abs(y)*math.tan(math.radians(35)), y=y,
            z=abs(y)*math.tan(math.radians(dihedral)), chord=1.5-.9*abs(y)/4, incidence_deg=0.) for y in [-4.,0.,4.]]
        decks = []
        for m,n in meshes:
            name = f'{ident}-M{m}-N{n}.avl'
            text = f'{ident} thin swept tapered wing\n0.0\n0 0 0.0\n8.4 1.1142857142857143 8.0\n0.375 0.0 0.0\n0.0\nSURFACE\nMainWing\n{m} 1.0\nCOMPONENT\n1\n'
            for s in sections:
                text += 'SECTION\n'+' '.join(format(s[k],'.17g') for k in ['x_le','y','z','chord','incidence_deg'])+f' {n//2} 1.0\n'
            path = fixtures/name
            path.write_text(text,encoding='ascii')
            decks.append(dict(m=m,n=n,path=path.relative_to(OUT).as_posix(),sha256=sha(path.read_bytes())))
        cases.append(dict(id=ident,leading_edge_sweep_deg=35.,dihedral_deg=dihedral,
            sections=sections,reference_area_xy=8.4,reference_span_y=8.,
            reference_chord=2/3*1.5*(1+.4+.4**2)/(1+.4),moment_origin=[.375,0.,0.],decks=decks))
    result = dict(kind='design_reference_fixture_not_product_schema',
        specification='finite-wing-load-reference-specification.md',solver_runs_performed=False,
        generator_sha256=sha(Path(__file__).read_bytes()),reference=dict(
            version='AVL 3.52',url='https://web.mit.edu/drela/Public/web/avl/avl3.52.tgz',
            archive_sha256=EXPECTED,archive_bytes=archive.stat().st_size,audited_files=provenance,
            executable_status='not acquired or built',source_core_audit_copy_status='not created',
            required_trailing_leg_forces=True,required_geometric_moment_axes=True,
            wake_direction=[1.,0.,0.]),
        density=1.2,speed=20.,mach=0.,sideslip_deg=0.,angles_deg=[-2.,-1.,0.,1.,2.],
        strip_boundaries_y=[-4+.4*i for i in range(21)],wake_lengths_over_reference_chord=[20,40,80,160],
        cases=cases,limitations=['decks not parsed by AVL yet','no reference numerical outputs',
            'reference option session, residual and conservative strip extraction still required'])
    (OUT/'finite-wing-reference-inputs.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(cases=len(cases),decks=sum(len(c['decks']) for c in cases),audited_files=len(provenance),solver_runs=0)))

if __name__=='__main__': main()
