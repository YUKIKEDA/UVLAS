"""Separate capacity-only external reference and frozen N320 input preparation."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import difflib,hashlib,json,shlex,shutil

ROOT=rd_paths.ROOT;BASE=ROOT/'build/reference-avl';FOLDER=rd_paths.RECORDS
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    original=BASE/'source-core0-grid';target=BASE/'source-core0-grid8000'
    shutil.copytree(original,target);changes=[]
    for name,before,after in [('avl_heap_inc.f90','NVX=5000','NVX=8000'),('ADIMEN.INC','NVMAX=5000','NVMAX=8000')]:
        p=target/'src'/name;text=p.read_text();assert text.count(before)==1
        new=text.replace(before,after);p.write_text(new)
        changes.extend(difflib.unified_diff(text.splitlines(True),new.splitlines(True),fromfile='grid/src/'+name,tofile='grid8000/src/'+name))
    patch=FOLDER/'fixtures/avl-measurement/capacity8000.patch';patch.write_text(''.join(changes))
    source_hashes={};changed=[]
    for p in (original/'src').iterdir():
        if p.suffix.lower() in ['.f','.f90','.inc']:
            source_hashes[p.name]=dict(original=digest(p),capacity8000=digest(target/'src'/p.name))
            if source_hashes[p.name]['original']!=source_hashes[p.name]['capacity8000']:changed.append(p.name)
    assert sorted(changed)==['ADIMEN.INC','avl_heap_inc.f90']
    linux='<repo>';prefix=linux+'/build/reference-avl/toolchain';directory=linux+'/build/reference-avl/source-core0-grid8000'
    link=[line for line in (BASE/'build-grid-export.sh').read_text().splitlines() if ' -o avl ' in line][-1]
    objects=[p for p in shlex.split(link) if p.endswith('.o')]
    compiler=prefix+'/bin/x86_64-conda-linux-gnu-gfortran'
    script='#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+shlex.quote(prefix+'/lib')+'\ncd '+shlex.quote(directory+'/bin')+'\n'
    for obj in objects:
        stem=obj[:-2];extension='.f90' if stem in ['avl_heap_inc','avl_heap'] else '.f'
        assert (target/'src'/(stem+extension)).exists(),stem
        args=[compiler,'-c','-O','-fomit-frame-pointer','-fallow-argument-mismatch']
        if stem!='matrix-lapacksubs':args+=['-fdefault-real-8']
        args+=['-I'+directory+'/src','-I'+directory+'/bin',directory+'/src/'+stem+extension,'-o',obj]
        script+=shlex.join(args)+'\n'
    script+=link+'\n';build=BASE/'build-grid8000.sh';build.write_text(script,newline='\n')
    parent_path=FOLDER/'avl-band-aligned-inputs.json';parent=json.loads(parent_path.read_text());decks=[]
    folder=FOLDER/'fixtures/avl-fw-span320';folder.mkdir()
    for case in parent['cases']:
        text=f'{case["id"]} thin swept tapered wing\n0.0\n0 0 0.0\n8.4 1.1142857142857143 8.0\n0.375 0.0 0.0\n0.0\nSURFACE\nMainWing\n24 1.0\nCOMPONENT\n1\n'
        for s in case['sections']:text+='SECTION\n'+' '.join(format(s[k],'.17g') for k in ['x_le','y','z','chord','incidence_deg'])+' 16 1.0\n'
        p=folder/(case['id']+'-M24-N320.avl');p.write_bytes(text.replace('\n','\r\n').encode('ascii'))
        area=sum((b['y']-a['y'])*(a['chord']+b['chord'])/2 for a,b in zip(case['sections'][:-1],case['sections'][1:]));assert abs(area-8.4)<1e-13
        decks.append(dict(case=case['id'],m=24,n=320,path=p.relative_to(ROOT).as_posix(),sha256=digest(p),area_xy=area,sections=case['sections']))
    result=dict(protocol='rd/aero/avl-reference/records/avl-span320-reference-experiment.md',source_sha256=digest(Path(__file__)),
        parent_input_sha256=digest(parent_path),parent_executable_sha256=digest(original/'bin/avl'),
        changed_source_files=changed,source_hashes=source_hashes,patch_sha256=digest(patch),build_script_sha256=digest(build),
        objects_recompiled=objects,decks=decks,capacity=8000,compiler_precision_note='matrix-lapacksubs explicit double compiled without default-real-8')
    (FOLDER/'avl-span320-preparation.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
