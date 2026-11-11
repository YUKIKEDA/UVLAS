"""Prepare external AVL variants without changing original numerical code."""
from pathlib import Path
import difflib
import hashlib
import json
import shutil

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    archive=BASE/'avl3.52.tgz'
    assert digest(archive)=='0b588ecea9222f5b625d0af0c87ae31daf3cdba1532cf0bbb36f93d6e854849b'
    pristine=BASE/'source-pristine'
    variants=[]
    patch_dir=ROOT/'docs/verification/fixtures/avl-measurement'
    patch_dir.mkdir(parents=True,exist_ok=True)
    adapter=ROOT/'tools/design-experiments/reference_avl_measurement.f'
    for name in ['baseline','measured','core0']:
        target=BASE/('source-'+name)
        if not target.exists():
            shutil.copytree(pristine,target)
        changes=[]
        if name!='baseline':
            p=target/'src/aero.f';old=p.read_text();new=old.replace("      INCLUDE 'AVL.INC'", "      INCLUDE 'AVL.INC'\n      INTEGER IUV\n      COMMON /UVMEAS/ IUV",1)
            bound='          FGAM(3) = 2.0*GAM(I)*F(3)'
            trail='            FGAM(3) = 2.0*GAM(I)*F(3)'
            assert old.splitlines().count(bound)==1 and old.splitlines().count(trail)==1
            new=new.replace(bound+'\n',bound+"\n          IF(IUV.GT.0) CALL UVL_FORCE(I,J,0,RROT,FGAM)\n",1)
            new=new.replace(trail+'\n',trail+"\n            IF(IUV.GT.0) CALL UVL_FORCE(I,J,ILEG,RROT,FGAM)\n",1)
            p.write_text(new);changes.append(('src/aero.f',old,new))
            p=target/'src/aoper.f';old=p.read_text()
            anchor="      ELSE IF(COMAND.EQ.'MRF ') THEN"
            assert old.count(anchor)==1
            new=old.replace(anchor,"      ELSE IF(COMAND.EQ.'UD  ') THEN\n        IF(LSOL) CALL UVL_DUMP\nC\n"+anchor,1)
            p.write_text(new);changes.append(('src/aoper.f',old,new))
            shutil.copyfile(adapter,target/'src/uvl_measurement.f')
        if name=='core0':
            p=target/'src/aic.f';old=p.read_text();anchor='            RCORE = 0.0001*DSYZ'
            assert old.count(anchor)==1
            new=old.replace(anchor,'            RCORE = 0.0*DSYZ',1)
            p.write_text(new);changes.append(('src/aic.f',old,new))
        diff=''.join(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='original/'+file,tofile=name+'/'+file)) for file,old,new in changes)
        patch_path=patch_dir/(name+'.patch');patch_path.write_text(diff)
        variants.append(dict(name=name,path=target.relative_to(ROOT).as_posix(),
            patch=patch_path.relative_to(ROOT).as_posix(),patch_sha256=digest(patch_path),
            changes=[dict(path=file,original_sha256=hashlib.sha256(old.encode()).hexdigest(),
                modified_sha256=digest(target/file)) for file,old,new in changes]))
    manifest=dict(protocol='docs/verification/avl-source-measurement-experiment.md',
        source_sha256=digest(Path(__file__)),adapter_sha256=digest(adapter),archive_sha256=digest(archive),
        variants=variants,license='External AVL and measurement adapter: GPL-2.0-or-later; not linked into UVLAS')
    (ROOT/'docs/verification/avl-source-preparation.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(variants=3,adapter_sha256=manifest['adapter_sha256'])))

if __name__=='__main__':main()
