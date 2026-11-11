import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import difflib
import hashlib
import json
import shutil
import shlex

ROOT=rd_paths.ROOT;BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    original=BASE/'source-core0-stable';target=BASE/'source-core0-grid'
    shutil.copytree(original,target)
    p=target/'src/aoper.f';old=p.read_text()
    marker="      ELSE IF(COMAND.EQ.'UD  ') THEN"
    assert old.count(marker)==1
    new=old.replace(marker,"      ELSE IF(COMAND.EQ.'UG  ') THEN\n        IF(LSOL) CALL UVL_GRID_DUMP(COMARG)\nC\n"+marker)
    p.write_text(new)
    exporter=ROOT/'build/third-party-gpl/avl/reference_avl_grid_export.f'
    shutil.copyfile(exporter,target/'src/uvl_grid_export.f')
    linux='<repo>';prefix=linux+'/build/reference-avl/toolchain'
    directory=linux+'/build/reference-avl/source-core0-grid'
    script='#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+shlex.quote(prefix+'/lib')+'\n'
    for name in ['aoper','uvl_grid_export']:
        args=[prefix+'/bin/x86_64-conda-linux-gnu-gfortran','-c','-O','-fomit-frame-pointer','-fdefault-real-8',
            '-fallow-argument-mismatch','-I'+directory+'/src','-I'+directory+'/bin',directory+'/src/'+name+'.f','-o',directory+'/bin/'+name+'.o']
        script+=shlex.join(args)+'\n'
    line=(BASE/'build-pairprobe.sh').read_text().splitlines()[-1]
    parts=shlex.split(line);parts.insert(parts.index('uvl_pair_probe.o')+1,'uvl_grid_export.o')
    script+='cd '+shlex.quote(directory+'/bin')+'\n'+shlex.join(parts)+'\n'
    build=BASE/'build-grid-export.sh';build.write_text(script,newline='\n')
    patch=ROOT/'build/third-party-gpl/avl/grid-export.patch'
    patch.write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='stable/src/aoper.f',tofile='grid/src/aoper.f')))
    result=dict(protocol='rd/aero/avl-reference/records/avl-fixed-band-convergence-experiment.md',source_sha256=digest(Path(__file__)),
        exporter_sha256=digest(exporter),original_source_sha256=hashlib.sha256(old.encode()).hexdigest(),
        modified_source_sha256=digest(p),unchanged_kernel_sha256=digest(target/'src/aic.f'),
        build_script_sha256=digest(build),patch_sha256=digest(patch),original_executable_sha256=digest(original/'bin/avl'))
    (ROOT/'rd/aero/avl-reference/records/avl-grid-export-preparation.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
