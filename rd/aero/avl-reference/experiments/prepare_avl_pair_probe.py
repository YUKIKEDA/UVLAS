import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib
import json
import shutil
import shlex

ROOT=rd_paths.ROOT
BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    original=BASE/'source-core0-corrected';target=BASE/'source-core0-pairprobe'
    shutil.copytree(original,target)
    adapter=target/'src/uvl_measurement.f'
    old=adapter.read_text();marker="      WRITE(*,*) 'UVL: reference rows, matrix, forces written'"
    assert old.count(marker)==1
    adapter.write_text(old.replace(marker,'      CALL UVL_PAIR_DUMP\n'+marker))
    probe=ROOT/'build/third-party-gpl/avl/reference_avl_pair_probe.f'
    shutil.copyfile(probe,target/'src/uvl_pair_probe.f')
    link=next(line for line in (BASE/'build-core0-corrected.log').read_text().splitlines() if ' -o avl ' in line)
    parts=shlex.split(link);parts.insert(parts.index('uvl_measurement.o')+1,'uvl_pair_probe.o')
    # The reference root and dependency library paths stay within the same workspace.
    linux='<repo>'
    prefix=linux+'/build/reference-avl/toolchain'
    directory=linux+'/build/reference-avl/source-core0-pairprobe'
    script='#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+shlex.quote(prefix+'/lib')+'\n'
    for name in ['uvl_measurement','uvl_pair_probe']:
        args=[prefix+'/bin/x86_64-conda-linux-gnu-gfortran','-c','-O','-fdefault-real-8','-fallow-argument-mismatch',
            '-I'+directory+'/src','-I'+directory+'/bin',directory+'/src/'+name+'.f','-o',directory+'/bin/'+name+'.o']
        script+=shlex.join(args)+'\n'
    script+='cd '+shlex.quote(directory+'/bin')+'\n'+shlex.join(parts)+'\n'
    build=BASE/'build-pairprobe.sh';build.write_text(script,newline='\n')
    result=dict(protocol='rd/aero/avl-reference/records/avl-collinear-kernel-experiment.md',source_sha256=digest(Path(__file__)),
        original_adapter_sha256=hashlib.sha256(old.encode()).hexdigest(),modified_adapter_sha256=digest(adapter),
        probe_sha256=digest(probe),build_script_sha256=digest(build),
        original_executable_sha256=digest(original/'bin/avl'))
    (ROOT/'rd/aero/avl-reference/records/avl-pair-probe-preparation.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
