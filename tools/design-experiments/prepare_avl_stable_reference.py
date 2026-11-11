from pathlib import Path
import difflib
import hashlib
import json
import shutil
import shlex

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'build/reference-avl'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    records=[];script='#!/bin/sh\nset -eu\n'
    linux='<repo>';prefix=linux+'/build/reference-avl/toolchain'
    script+='export LD_LIBRARY_PATH='+shlex.quote(prefix+'/lib')+'\n'
    for name,parent,log in [('measured','source-measured-corrected','build-measured-corrected.log'),
        ('core0','source-core0-pairprobe','build-pairprobe.log')]:
        original=BASE/parent;target=BASE/f'source-{name}-stable';shutil.copytree(original,target)
        path=target/'src/aic.f';old=path.read_text()
        before,body=old.split('      SUBROUTINE VORVELC',1)
        body=body.replace('      LOGICAL LBOUND','      LOGICAL LBOUND,LSTABLE',1)
        body=body.replace('      REAL A(3), B(3), AXB(3)',
            '      REAL A(3), B(3), AXB(3)\n      REAL ED(3),EH(3),RP(3),ELEN,QLEN,SI',1)
        anchor='      IF (LBOUND  .AND.  AMAG*BMAG .NE. 0.0) THEN'
        addition='''
C     UVLAS reference-only stable exterior zero-core finite segment.
        LSTABLE=.FALSE.
        IF(RCORE.EQ.0.0.AND.BETA.EQ.1.0) THEN
          ED(1)=X2-X1
          ED(2)=Y2-Y1
          ED(3)=Z2-Z1
          ELEN=SQRT(ED(1)**2+ED(2)**2+ED(3)**2)
          IF(ELEN.GT.0.0) THEN
            QLEN=0.0
            DO K=1,3
              EH(K)=ED(K)/ELEN
              QLEN=QLEN-A(K)*EH(K)
            ENDDO
            LSTABLE=QLEN.LT.0.0.OR.QLEN.GT.ELEN
          ENDIF
        ENDIF
        IF(LSTABLE) THEN
          DO K=1,3
            RP(K)=-A(K)-QLEN*EH(K)
          ENDDO
          SI=ELEN*(2.0*QLEN-ELEN)/
     &       (AMAG*BMAG*(QLEN*BMAG+(QLEN-ELEN)*AMAG))
          U=(EH(2)*RP(3)-EH(3)*RP(2))*SI
          V=(EH(3)*RP(1)-EH(1)*RP(3))*SI
          W=(EH(1)*RP(2)-EH(2)*RP(1))*SI
        ELSE
'''
        assert body.count(anchor)==1
        body=body.replace(anchor,anchor+addition,1)
        end='         W = AXB(3)*T\n        ENDIF\n      ENDIF'
        assert body.count(end)==1
        body=body.replace(end,'         W = AXB(3)*T\n        ENDIF\n        ENDIF\n      ENDIF',1)
        new=before+'      SUBROUTINE VORVELC'+body;path.write_text(new)
        patch=ROOT/f'docs/verification/fixtures/avl-measurement/{name}-stable.patch'
        patch.write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=parent+'/src/aic.f',tofile=target.name+'/src/aic.f')))
        directory=linux+'/build/reference-avl/'+target.name
        args=[prefix+'/bin/x86_64-conda-linux-gnu-gfortran','-c','-O','-fomit-frame-pointer','-fdefault-real-8',
            '-fallow-argument-mismatch','-I'+directory+'/src','-I'+directory+'/bin',directory+'/src/aic.f','-o',directory+'/bin/aic.o']
        script+=shlex.join(args)+'\ncd '+shlex.quote(directory+'/bin')+'\n'
        # Pair-probe link command is generated and recorded verbatim by preparation.
        if name=='core0':
            line=(BASE/'build-pairprobe.sh').read_text().splitlines()[-1]
        else:line=next(x for x in (BASE/log).read_text().splitlines() if ' -o avl ' in x)
        script+=line+'\n'
        records.append(dict(name=name,parent_path=original.relative_to(ROOT).as_posix(),path=target.relative_to(ROOT).as_posix(),
            original_source_sha256=hashlib.sha256(old.encode()).hexdigest(),modified_source_sha256=digest(path),
            patch=patch.relative_to(ROOT).as_posix(),patch_sha256=digest(patch)))
    build=BASE/'build-stable-reference.sh';build.write_text(script,newline='\n')
    result=dict(protocol='docs/verification/avl-collinear-kernel-experiment.md#6',source_sha256=digest(Path(__file__)),
        variants=records,build_script_sha256=digest(build))
    (ROOT/'docs/verification/avl-stable-reference-preparation.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
