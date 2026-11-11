"""Prepare isolated native symmetry measurement and positive-side input decks."""
from pathlib import Path
import hashlib,json,shlex,shutil,difflib

ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification';BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parent=BASE/'source-core0-grid8000';target=BASE/'source-symmetry8000';shutil.copytree(parent,target)
    p=target/'src/aoper.f';old=p.read_text();marker="      ELSE IF(COMAND.EQ.'UG  ') THEN";assert old.count(marker)==1
    new=old.replace(marker,"      ELSE IF(COMAND.EQ.'UH  ') THEN\n        IF(LSOL) CALL UVL_SYMMETRY_DUMP(COMARG)\nC\n"+marker);p.write_text(new)
    adapter_names=['reference_avl_symmetry_export.f','reference_avl_symmetry_guard.f']
    objects=['uvl_symmetry_export','uvl_symmetry_guard']
    for name,obj in zip(adapter_names,objects):
        src=Path(__file__).with_name(name)
        assert all(len(line)<=72 for line in src.read_text().splitlines() if line and line[0] not in 'Cc*!')
        shutil.copyfile(src,target/'src'/(obj+'.f'))
    patch=FOLDER/'fixtures/avl-measurement/symmetry-export.patch'
    patch.write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='grid8000/src/aoper.f',tofile='symmetry/src/aoper.f')))
    linux='<repo>';prefix=linux+'/build/reference-avl/toolchain';directory=linux+'/build/reference-avl/source-symmetry8000'
    compiler=prefix+'/bin/x86_64-conda-linux-gnu-gfortran'
    script='#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+shlex.quote(prefix+'/lib')+'\ncd '+shlex.quote(directory+'/bin')+'\n'
    for name in ['aoper']+objects:
        args=[compiler,'-c','-O','-fomit-frame-pointer','-fallow-argument-mismatch','-fdefault-real-8',
            '-I'+directory+'/src','-I'+directory+'/bin',directory+'/src/'+name+'.f','-o',name+'.o']
        script+=shlex.join(args)+'\n'
    link=shlex.split([line for line in (BASE/'build-grid8000.sh').read_text().splitlines() if ' -o avl ' in line][-1])
    pos=link.index('uvl_grid_export.o')+1;link[pos:pos]=[name+'.o' for name in objects]
    script+=shlex.join(link)+'\n';build=BASE/'build-symmetry.sh';build.write_text(script,newline='\n')
    aligned_path=FOLDER/'avl-band-aligned-inputs.json';aligned=json.loads(aligned_path.read_text())
    oldruns=json.loads((FOLDER/'avl-band-aligned-grid-runs.json').read_text())['runs']
    larger=json.loads((FOLDER/'avl-span320-runs.json').read_text())['runs']
    folder=FOLDER/'fixtures/avl-fw-symmetric';folder.mkdir();decks=[]
    for case in aligned['cases']:
        assert len(case['sections'])==21 and case['sections'][10]['y']==0
        for m,n in [(6,40),(24,160),(24,320),(24,640)]:
            title=f'{case["id"]} thin swept tapered wing'
            text=title+'\n0.0\n1 0 0.0\n8.4 1.1142857142857143 8.0\n0.375 0.0 0.0\n0.0\nSURFACE\nMainWing\n'+f'{m} 1.0\nCOMPONENT\n1\n'
            sections=case['sections'][10:]
            for s in sections:text+='SECTION\n'+' '.join(format(s[k],'.17g') for k in ['x_le','y','z','chord','incidence_deg'])+f' {n//20} 1.0\n'
            p=folder/f'{case["id"]}-M{m}-N{n}-half.avl';p.write_bytes(text.replace('\n','\r\n').encode('ascii'))
            record=dict(case=case['id'],m=m,physical_n=n,native_n=n//2,path=p.relative_to(ROOT).as_posix(),sha256=digest(p),sections=sections,
                half_area=sum((b['y']-a['y'])*(a['chord']+b['chord'])/2 for a,b in zip(sections[:-1],sections[1:])))
            assert abs(record['half_area']-4.2)<1e-13
            if n<=160:
                r=next(r for r in oldruns if (r['case'],r['m'],r['n'],r['alpha'])==(case['id'],m,n,1))
                record['comparison_directory']=r['directory'];record['full_input']=(ROOT/r['directory']).parent.joinpath('wing.avl').relative_to(ROOT).as_posix()
            elif n==320:
                r=next(r for r in larger if r['kind']=='new' and r['case']==case['id'])
                record['comparison_directory']=r['directory']+'/out';record['full_input']=r['directory']+'/wing.avl'
            decks.append(record)
    sources={}
    for source in (parent/'src').iterdir():
        if source.suffix.lower() in ['.f','.f90','.inc']:sources[source.name]=dict(parent=digest(source),measured=digest(target/'src'/source.name))
    assert [name for name,v in sources.items() if v['parent']!=v['measured']]==['aoper.f']
    result=dict(protocol='docs/verification/avl-symmetry-reduction-experiment.md',source_sha256=digest(Path(__file__)),
        parent_input_sha256=digest(aligned_path),parent_executable_sha256=digest(parent/'bin/avl'),source_hashes=sources,
        adapter_sha256={name:digest(Path(__file__).with_name(name)) for name in adapter_names},patch_sha256=digest(patch),build_script_sha256=digest(build),decks=decks)
    (FOLDER/'avl-symmetry-preparation.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
