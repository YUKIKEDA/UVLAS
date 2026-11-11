"""Separate local geometry capacity correction; preserve failed first fine attempt."""
from pathlib import Path
import difflib,hashlib,json,shlex,shutil
ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification';BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=FOLDER/'avl-span1280-capacity-ks-preparation.json'
    if output.exists():raise FileExistsError(output)
    old=json.loads((FOLDER/'avl-span1280-capacity-preparation.json').read_text())
    failed=FOLDER/'avl-span1280-fine-runs.json';failure=json.loads(failed.read_text());assert not failure['completed']
    parent=ROOT/old['target'];target=BASE/'source-symmetry16000-ks1000-spacing-pi64';shutil.copytree(parent,target)
    p=target/'src/amake.f';before=p.read_bytes();assert before.count(b'KSMAX=500')==1;p.write_bytes(before.replace(b'KSMAX=500',b'KSMAX=1000'))
    patch=FOLDER/'fixtures/avl-measurement/capacity-ks1000.patch';patch.write_text(''.join(difflib.unified_diff(before.decode().splitlines(True),p.read_text().splitlines(True),fromfile='parent/src/amake.f',tofile='ks1000/src/amake.f')))
    sources={p.name:dict(parent=digest(p),capacity=digest(target/'src'/p.name)) for p in (parent/'src').iterdir() if p.suffix.lower() in ['.f','.f90','.inc']}
    assert [n for n,h in sources.items() if h['parent']!=h['capacity']]==['amake.f']
    linux='<repo>';directory=linux+'/build/reference-avl/'+target.name
    compiler=linux+'/build/reference-avl/toolchain/bin/x86_64-conda-linux-gnu-gfortran'
    link=[s for s in (BASE/'build-symmetry16000.sh').read_text().splitlines() if ' -o avl ' in s][-1]
    script='#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+linux+'/build/reference-avl/toolchain/lib\ncd '+directory+'/bin\n'
    script+=shlex.join([compiler,'-c','-O','-fomit-frame-pointer','-fallow-argument-mismatch','-fdefault-real-8','-I'+directory+'/src',directory+'/src/amake.f','-o','amake.o'])+'\n'+link+'\n'
    build=BASE/'build-symmetry16000-ks1000.sh';build.write_text(script,newline='\n')
    parentruns=FOLDER/'avl-span1280-capacity-runs.json';runs=json.loads(parentruns.read_text());assert runs['completed']
    inputs=[dict(case=r['case'],m=r['m'],physical_n=r['physical_n'],native_unknowns=r['native_unknowns'],
        input=r['directory']+'/wing.avl',input_sha256=digest(ROOT/r['directory']/'wing.avl'),parent_directory=r['directory']+'/out') for r in runs['runs']]
    result=dict(old,protocol='docs/verification/avl-span1280-ksmax-addendum.md',protocol_sha256=digest(FOLDER/'avl-span1280-ksmax-addendum.md'),
        source_sha256=digest(__file__),parent=parent.relative_to(ROOT).as_posix(),target=target.relative_to(ROOT).as_posix(),source_files=sources,
        changed_source_files=['amake.f'],changed_constants=1,parent_executable_sha256=digest(parent/'bin/avl'),
        build_script_sha256=digest(build),patch_sha256=digest(patch),inputs=inputs,parent_runs_sha256=digest(parentruns),failed_runs_sha256=digest(failed),
        reused_object_sha256={n:digest(parent/'bin'/n) for n in old['objects'] if n!='amake.o'})
    output.write_text(json.dumps(result,indent=2)+'\n')
    # Generate separate runners, keeping all executed sources and old manifests unchanged.
    runner=Path(__file__).with_name('run_avl_span1280_capacity.py').read_text()
    for a,b in [('avl-span1280-capacity-runs','avl-span1280-capacity-ks-runs'),('avl-span1280-capacity-preparation','avl-span1280-capacity-ks-preparation'),('build-symmetry16000.','build-symmetry16000-ks1000.'),('avl-capacity16000-','avl-capacity16000-ks1000-')]:runner=runner.replace(a,b)
    Path(__file__).with_name('run_avl_span1280_capacity_ks.py').write_text(runner)
    fine=json.loads((FOLDER/'avl-span1280-fine-preparation.json').read_text());fine.update(protocol='docs/verification/avl-span1280-fine-experiment.md',
        correction_protocol_sha256=digest(FOLDER/'avl-span1280-ksmax-addendum.md'),source_sha256=digest(__file__),
        capacity_preparation_sha256=digest(output),failed_runs_sha256=digest(failed))
    # The fine preparation is finalized after the new regression audit and executable exist.
    (FOLDER/'avl-span1280-fine-ks-pending.json').write_text(json.dumps(fine,indent=2)+'\n')
    print(json.dumps(dict(changed='amake.f KSMAX',recompiled=1,reused=40)))
if __name__=='__main__':main()
