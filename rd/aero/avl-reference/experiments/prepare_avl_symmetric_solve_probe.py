"""Separate read-only UQ adapter; keep the parent and all accepted artifacts intact."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import difflib,hashlib,json,shlex,shutil
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=F/'avl-symmetric-solve-preparation.json'
    if output.exists():raise FileExistsError(output)
    weights=F/'avl-symmetric-solve-weight-result.json';w=json.loads(weights.read_text());assert w['all_criteria_pass']
    parent=BASE/'source-symmetry16000-ks1000-spacing-pi64';target=BASE/'source-symmetry16000-solve-probe';shutil.copytree(parent,target)
    p=target/'src/aoper.f';before=p.read_bytes();marker=b"      ELSE IF(COMAND.EQ.'UH  ') THEN";assert before.count(marker)==1
    nl=b'\r\n' if b'\r\n' in before else b'\n';insert=nl.join([b"      ELSE IF(COMAND.EQ.'UQ  ') THEN",b'        IF(LSOL) CALL UVL_SYM_SOLVE_PROBE(COMARG)',b'C'])+nl
    after=before.replace(marker,insert+marker);p.write_bytes(after)
    adapter=rd_paths.sibling('reference_avl_sym_solve_probe.f');shutil.copyfile(adapter,target/'src/uvl_sym_solve_probe.f')
    patch=F/'fixtures/avl-measurement/symmetric-solve-probe.patch';patch.write_text(''.join(difflib.unified_diff(before.decode().splitlines(True),after.decode().splitlines(True),fromfile='parent/src/aoper.f',tofile='solve-probe/src/aoper.f')))
    sources={p.name:dict(parent=digest(p),probe=digest(target/'src'/p.name)) for p in (parent/'src').iterdir() if p.suffix.lower() in ['.f','.f90','.inc']}
    assert [n for n,h in sources.items() if h['parent']!=h['probe']]==['aoper.f']
    linux='<repo>';directory=linux+'/build/reference-avl/'+target.name;prefix=linux+'/build/reference-avl/toolchain';compiler=prefix+'/bin/x86_64-conda-linux-gnu-gfortran'
    link=[s for s in (BASE/'build-symmetry16000-ks1000.sh').read_text().splitlines() if ' -o avl ' in s][-1];parts=shlex.split(link);objects=[p for p in parts if p.endswith('.o')]
    pos=parts.index('uvl_symmetry_guard.o')+1;parts[pos:pos]=['uvl_sym_solve_probe.o']
    script='#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+prefix+'/lib\ncd '+directory+'/bin\n'
    for name in ['aoper','uvl_sym_solve_probe']:
        script+=shlex.join([compiler,'-c','-O','-fomit-frame-pointer','-fallow-argument-mismatch','-fdefault-real-8','-I'+directory+'/src','-I'+directory+'/bin',directory+'/src/'+name+'.f','-o',name+'.o'])+'\n'
    script+=shlex.join(parts)+'\n';build=BASE/'build-symmetric-solve-probe.sh';build.write_text(script,newline='\n')
    result=dict(source_sha256=digest(__file__),protocol=w['protocol'],protocol_sha256=w['protocol_sha256'],weights_sha256=digest(weights),
        parent=parent.relative_to(ROOT).as_posix(),target=target.relative_to(ROOT).as_posix(),parent_executable_sha256=digest(parent/'bin/avl'),
        source_files=sources,adapter_sha256=digest(adapter),new_adapter='uvl_sym_solve_probe.f',build_script_sha256=digest(build),patch_sha256=digest(patch),
        recompiled_objects=['aoper.o','uvl_sym_solve_probe.o'],reused_object_sha256={o:digest(parent/'bin'/o) for o in objects if o!='aoper.o'},
        compiler_sha256=digest(Path(compiler.replace(linux,str(ROOT)))),link_inputs={p:digest(ROOT/p[len(linux)+1:]) for p in parts if p.endswith('.a')},cases=w['records'])
    for r in result['cases']:r['input_sha256']=digest(ROOT/r['input'])
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(existing_changed=1,new_adapter=1,recompiled=2,reused=len(result['reused_object_sha256']),cases=8)))
if __name__=='__main__':main()
