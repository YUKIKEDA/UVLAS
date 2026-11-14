"""Isolated SPACER pi64 correction, applied consistently to full and half references."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import difflib,hashlib,json,shlex,shutil

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS;BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    compiler='<repo>/build/reference-avl/toolchain/bin/x86_64-conda-linux-gnu-gfortran'
    linux='<repo>';prefix=linux+'/build/reference-avl/toolchain'
    script='#!/bin/sh\nset -eu\nexport LD_LIBRARY_PATH='+shlex.quote(prefix+'/lib')+'\n';variants=[]
    for source,build in [('source-core0-grid8000','build-grid8000.sh'),('source-symmetry8000','build-symmetry.sh')]:
        parent=BASE/source;target=BASE/(source+'-spacing-pi64');shutil.copytree(parent,target)
        p=target/'src/sgutil.f';old=p.read_text();marker='      DATA  PI / 3.1415926535 /';assert old.count(marker)==1
        new=old.replace(marker,'      DATA  PI / 3.1415926535897932384626433832795 /');p.write_text(new)
        patch=FOLDER/'fixtures/avl-measurement'/(source+'-spacing-pi64.patch')
        patch.write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=source+'/src/sgutil.f',tofile=target.name+'/src/sgutil.f')))
        directory=linux+'/build/reference-avl/'+target.name
        script+='cd '+shlex.quote(directory+'/bin')+'\n'
        script+=shlex.join([compiler,'-c','-O','-fomit-frame-pointer','-fallow-argument-mismatch','-fdefault-real-8',directory+'/src/sgutil.f','-o','sgutil.o'])+'\n'
        link=[line for line in (BASE/build).read_text().splitlines() if ' -o avl ' in line][-1];script+=link+'\n'
        hashes={s.name:dict(parent=digest(s),pi64=digest(target/'src'/s.name)) for s in (parent/'src').iterdir() if s.suffix.lower() in ['.f','.f90','.inc']}
        assert [k for k,v in hashes.items() if v['parent']!=v['pi64']]==['sgutil.f']
        variants.append(dict(parent=source,target=target.name,sources=hashes,patch_path=patch.relative_to(ROOT).as_posix(),patch_sha256=digest(patch)))
    for label,source in [('old','source-core0-grid8000'),('pi64','source-core0-grid8000-spacing-pi64')]:
        driver=linux+'/tools/design-experiments/reference_avl_spacing_driver.f';executable=linux+'/build/reference-avl/spacing-'+label
        args=[compiler,'-O','-fdefault-real-8','-fallow-argument-mismatch',driver,linux+'/build/reference-avl/'+source+'/src/sgutil.f','-o',executable]
        script+=shlex.join(args)+'\n'+shlex.quote(executable)+' > '+shlex.quote(linux+'/build/reference-avl/spacing-'+label+'.txt')+'\n'
    build=BASE/'build-symmetry-pi64.sh';build.write_text(script,newline='\n')
    parent_path=FOLDER/'avl-symmetry-preparation.json';parent=json.loads(parent_path.read_text());decks=list(parent['decks'])
    aligned=json.loads((FOLDER/'avl-band-aligned-inputs.json').read_text());oldruns=json.loads((FOLDER/'avl-band-aligned-grid-runs.json').read_text())['runs']
    folder=FOLDER/'fixtures/avl-fw-symmetric-pi';folder.mkdir()
    for case in aligned['cases']:
        donor=next(d for d in parent['decks'] if d['case']==case['id'] and d['physical_n']==160)
        content=(ROOT/donor['path']).read_bytes();assert content.count(b' 8 1.0\r\n')==11
        p=folder/(case['id']+'-M24-N80-half.avl');p.write_bytes(content.replace(b' 8 1.0\r\n',b' 4 1.0\r\n'))
        old=next(r for r in oldruns if (r['case'],r['m'],r['n'],r['alpha'])==(case['id'],24,80,1))
        record=dict(donor,m=24,physical_n=80,native_n=40,path=p.relative_to(ROOT).as_posix(),sha256=digest(p),
            comparison_directory=old['directory'],full_input=(ROOT/old['directory']).parent.joinpath('wing.avl').relative_to(ROOT).as_posix())
        decks.append(record)
    result=dict(protocol='rd/aero/avl-reference/records/avl-symmetry-spacing-pi-addendum.md',source_sha256=digest(Path(__file__)),
        parent_preparation_sha256=digest(parent_path),variants=variants,build_script_sha256=digest(build),decks=decks,
        driver_sha256=digest(rd_paths.sibling('reference_avl_spacing_driver.f')),reference_id='AVL3.52/DP/zero-core-exterior-stable/spacing-pi64-v1')
    (FOLDER/'avl-symmetry-pi-preparation.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
