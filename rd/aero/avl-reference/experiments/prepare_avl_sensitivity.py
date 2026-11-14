"""Prepare a separate external-reference copy; never change the accepted UG build."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import difflib
import hashlib
import json
import shlex
import shutil

ROOT=rd_paths.ROOT
BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    original=BASE/'source-core0-grid';target=BASE/'source-core0-sensitivity-corrected'
    shutil.copytree(original,target)
    p=target/'src/aoper.f';old=p.read_text()
    marker="      ELSE IF(COMAND.EQ.'UG  ') THEN"
    assert old.count(marker)==1
    new=old.replace(marker,"      ELSE IF(COMAND.EQ.'US  ') THEN\n        IF(LSOL) CALL UVL_SENSITIVITY(COMARG)\nC\n"+marker)
    p.write_text(new)
    adapter=ROOT/'build/third-party-gpl/avl/reference_avl_sensitivity.f'
    shutil.copyfile(adapter,target/'src/uvl_sensitivity.f')
    linux='<repo>';prefix=linux+'/build/reference-avl/toolchain'
    directory=linux+'/build/reference-avl/source-core0-sensitivity-corrected'
    compiler=prefix+'/bin/x86_64-conda-linux-gnu-gfortran'
    script='#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+shlex.quote(prefix+'/lib')+'\n'
    for name in ['aoper','uvl_sensitivity']:
        args=[compiler,'-c','-O','-fomit-frame-pointer','-fdefault-real-8','-fallow-argument-mismatch',
              '-I'+directory+'/src','-I'+directory+'/bin',directory+'/src/'+name+'.f','-o',directory+'/bin/'+name+'.o']
        script+=shlex.join(args)+'\n'
    for name in ['dlacn2','dasum','dcopy']:
        args=[compiler,'-c','-O',linux+'/build/reference-avl/norm-estimator/'+name+'.f','-o',directory+'/bin/'+name+'.o']
        script+=shlex.join(args)+'\n'
    line=[s for s in (BASE/'build-grid-export.sh').read_text().splitlines() if ' -o avl ' in s][-1]
    parts=shlex.split(line)
    pos=parts.index('uvl_grid_export.o')+1
    parts[pos:pos]=['uvl_sensitivity.o','dlacn2.o','dasum.o','dcopy.o']
    script+='cd '+shlex.quote(directory+'/bin')+'\n'+shlex.join(parts)+'\n'
    build=BASE/'build-sensitivity-corrected.sh';build.write_text(script,newline='\n')
    patch=ROOT/'build/third-party-gpl/avl/solve-sensitivity-corrected.patch'
    patch.write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='grid/src/aoper.f',tofile='sensitivity/src/aoper.f')))
    result=dict(protocol='rd/aero/avl-reference/records/avl-solve-sensitivity-experiment.md',source_sha256=digest(Path(__file__)),
        adapter_sha256=digest(adapter),parent_source_sha256=digest(original/'src/aoper.f'),
        modified_source_sha256=digest(p),unchanged_kernel_sha256=digest(target/'src/aic.f'),
        build_script_sha256=digest(build),patch_sha256=digest(patch),parent_executable_sha256=digest(original/'bin/avl'),
        lapack_manifest_sha256=digest(BASE/'norm-estimator/source-manifest.json'))
    (ROOT/'rd/aero/avl-reference/records/avl-sensitivity-preparation-corrected.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
