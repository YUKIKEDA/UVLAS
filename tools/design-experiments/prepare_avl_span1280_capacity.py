"""Separate capacity-only reference, with complete recompilation of all linked objects."""
from pathlib import Path
import difflib,hashlib,json,shlex,shutil

ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification';BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    output=FOLDER/'avl-span1280-capacity-preparation.json'
    if output.exists():raise FileExistsError(output)
    parent=BASE/'source-symmetry8000-spacing-pi64';target=BASE/'source-symmetry16000-spacing-pi64'
    shutil.copytree(parent,target);diff=[]
    replacements={'avl_heap_inc.f90':[(b'NVX=8000',b'NVX=16000')],
        'ADIMEN.INC':[(b'NVMAX=8000',b'NVMAX=16000'),(b'NSMAX=500',b'NSMAX=1000')]}
    for name,changes in replacements.items():
        p=target/'src'/name;before=p.read_bytes();after=before
        for old,new in changes:assert after.count(old)==1;after=after.replace(old,new)
        p.write_bytes(after)
        diff.extend(difflib.unified_diff(before.decode().splitlines(True),after.decode().splitlines(True),fromfile='parent/src/'+name,tofile='capacity16000/src/'+name))
    patch=FOLDER/'fixtures/avl-measurement/capacity16000.patch';patch.write_text(''.join(diff))
    sources={p.name:dict(parent=digest(p),capacity=digest(target/'src'/p.name)) for p in (parent/'src').iterdir() if p.suffix.lower() in ['.f','.f90','.inc']}
    changed=sorted(n for n,v in sources.items() if v['parent']!=v['capacity']);assert changed==sorted(replacements)
    linux='<repo>';prefix=linux+'/build/reference-avl/toolchain';directory=linux+'/build/reference-avl/'+target.name
    link=[s for s in (BASE/'build-symmetry-pi64.sh').read_text().splitlines() if ' -o avl ' in s][-1]
    objects=[p for p in shlex.split(link) if p.endswith('.o')];compiler=prefix+'/bin/x86_64-conda-linux-gnu-gfortran'
    script='#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+shlex.quote(prefix+'/lib')+'\ncd '+shlex.quote(directory+'/bin')+'\n'
    for obj in objects:
        stem=obj[:-2];ext='.f90' if stem in ['avl_heap_inc','avl_heap'] else '.f';assert (target/'src'/(stem+ext)).exists()
        args=[compiler,'-c','-O','-fomit-frame-pointer','-fallow-argument-mismatch']
        if stem!='matrix-lapacksubs':args+=['-fdefault-real-8']
        args+=['-I'+directory+'/src','-I'+directory+'/bin',directory+'/src/'+stem+ext,'-o',obj]
        script+=shlex.join(args)+'\n'
    script+=link+'\n';build=BASE/'build-symmetry16000.sh';build.write_text(script,newline='\n')
    gate_path=FOLDER/'avl-symmetry-pi-equivalence-result.json';gate=json.loads(gate_path.read_text());assert gate['all_criteria_pass']
    fine_path=FOLDER/'avl-symmetry-pi-fine-runs.json';fine=json.loads(fine_path.read_text());assert fine['completed']
    inputs=[]
    for case in ['FW-V01','FW-V02']:
        for m,n in [(6,40),(24,160),(24,320),(24,640)]:
            if n==640:
                record=next(r for r in fine['runs'] if (r['case'],r['m'],r['physical_n'])==(case,m,n));folder=ROOT/record['directory'];raw=folder/'out'
            else:
                record=next(r for r in gate['records'] if (r['case'],r['m'],r['physical_n'])==(case,m,n));raw=ROOT/record['raw_directory'];folder=raw.parent
            p=folder/'wing.avl';inputs.append(dict(case=case,m=m,physical_n=n,native_unknowns=m*n//2,input=p.relative_to(ROOT).as_posix(),input_sha256=digest(p),parent_directory=raw.relative_to(ROOT).as_posix()))
    result=dict(protocol='docs/verification/avl-span1280-capacity-experiment.md',protocol_sha256=digest(FOLDER/'avl-span1280-capacity-experiment.md'),
        source_sha256=digest(Path(__file__)),parent=parent.relative_to(ROOT).as_posix(),target=target.relative_to(ROOT).as_posix(),
        parent_executable_sha256=digest(parent/'bin/avl'),source_files=sources,changed_source_files=changed,changed_constants=3,
        build_script_sha256=digest(build),patch_sha256=digest(patch),compiler_sha256=digest(Path(compiler.replace(linux,str(ROOT)))),
        objects=objects,link_inputs={p:digest(ROOT/p[len(linux)+1:]) for p in shlex.split(link) if p.endswith('.a')},
        gate_sha256=digest(gate_path),parent_runs_sha256=digest(fine_path),inputs=inputs,capacity=16000,strips=1000)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(changed=changed,objects=len(objects),inputs=len(inputs))))

if __name__=='__main__':main()
