"""Independent rational elimination audit; no finite-field arithmetic."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
import numpy as np

ROOT=rd_paths.ROOT


def exact_null(A):
    M=[row[:] for row in A];pivots=[];r=0;n=len(A[0])
    for col in range(n):
        pivot=next((j for j in range(r,len(M)) if M[j][col]),None)
        if pivot is None:continue
        M[r],M[pivot]=M[pivot],M[r]
        value=M[r][col];M[r]=[v/value for v in M[r]]
        active=[j for j,v in enumerate(M[r]) if v]
        for i in range(len(M)):
            if i==r or not M[i][col]:continue
            factor=M[i][col]
            for j in active:M[i][j]-=factor*M[r][j]
        pivots.append(col);r+=1
    free=[j for j in range(n) if j not in pivots];vectors=[]
    for col in free:
        v=[F(0)]*n;v[col]=F(1)
        for i,pivot in enumerate(pivots):v[pivot]=-M[i][col]
        vectors.append(v)
    return r,free,vectors


def main():
    target=ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-flat-rank-audit.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();parent_path=ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-flat-rank-result.json';parent=json.loads(parent_path.read_text());records=[]
    for item in parent['records']:
        path=ROOT/item['artifact'];assert hashlib.sha256(path.read_bytes()).hexdigest()==item['artifact_sha256']
        saved=np.load(path);num,den=saved['numerator'],saved['denominator'];A=[[F(int(n),int(d)) for n,d in zip(nr,dr)] for nr,dr in zip(num,den)]
        rank,free,vectors=exact_null(A)
        zero=all(sum((a*b for a,b in zip(row,v)),F(0))==0 for v in vectors for row in A)
        independent=all(all(v[free[j]]==int(i==j) for j in range(len(free))) for i,v in enumerate(vectors))
        matches=all(v==rank for v in item['modular_ranks'].values())
        passed=zero and independent and matches and (rank==138 if item['space']=='selected_203' else rank<138)
        witness=path.with_name(path.stem+'-complete-rational-kernel.json')
        witness.write_text(json.dumps(dict(free_columns=free,vectors=[[[v.numerator,v.denominator] for v in row] for row in vectors])),encoding='utf-8')
        records.append(dict(space=item['space'],rational_rank=rank,nullity=len(vectors),exact_zero=zero,independent=independent,matches_modular_lower_bound=matches,passed=passed,witness=witness.relative_to(ROOT).as_posix(),witness_sha256=hashlib.sha256(witness.read_bytes()).hexdigest()))
        print(records[-1],flush=True)
    sources=[Path(__file__),parent_path,ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-flat-rank-audit-protocol.md',ROOT/'rd/structure/main-shell/gsp-q3/records/gsp-q3-flat-rank-audit-preparation.json']
    target.write_text(json.dumps(dict(records=records,all_checks_pass=all(v['passed'] for v in records),scope='affine_positive_thickness_only',sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')


if __name__=='__main__':main()
