"""Exact rational work matrices, modular lower bound and rigid-kernel upper bound."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,itertools,json,time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
COLS=list(itertools.product(range(3),range(3),range(4),range(4)))
PAIRS=[(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]


def stress_rows(reduced=False):
    shapes=[(2,3,2),(3,2,2),(3,3,1),(3,2,2),(2,3,2),(2,2,2)]
    return [(component,i,j,k) for component,bounds in enumerate(shapes)
        for i,j,k in itertools.product(*(range(v+1) for v in ((2,2,2) if reduced else bounds)))]


def integral(power):
    return F(0) if power%2 else F(2,power+1)


def matrix(rows):
    values=[]
    for component,a,b,c in rows:
        row=[];i,j=PAIRS[component]
        for vector,p,s,t in COLS:
            value=F(0)
            directions=[i] if i==j and vector==i else ([j] if vector==i else ([i] if vector==j else []))
            for direction in directions:
                powers=[s,t,p];factor=powers[direction]
                if not factor:continue
                powers[direction]-=1
                value+=F(factor,2 if p==2 else 1)*integral(a+powers[0])*integral(b+powers[1])*integral(c+powers[2])
            row.append(value)
        values.append(row)
    return values


def modular(values,p):
    return np.array([[v.numerator*pow(v.denominator,-1,p)%p for v in row] for row in values],dtype=np.int64)


def rank(A,p):
    A=A.copy();r=0
    for col in range(A.shape[1]):
        candidates=np.flatnonzero(A[r:,col])
        if not len(candidates):continue
        pivot=r+int(candidates[0]);A[[r,pivot]]=A[[pivot,r]]
        A[r]=A[r]*pow(int(A[r,col]),-1,p)%p
        indices=np.flatnonzero(A[:,col]);indices=indices[indices!=r]
        A[indices]=(A[indices]-A[indices,col,None]*A[r,None,:])%p
        r+=1
        if r==min(A.shape):break
    return r


def rigid():
    vectors=[]
    for a in range(3):
        v=[F(0)]*144;v[COLS.index((a,0,0,0))]=F(1);vectors.append(v)
    for entries in [[((1,1,0,0),-1),((2,0,0,1),1)],
                    [((0,1,0,0),1),((2,0,1,0),-1)],
                    [((0,0,0,1),-1),((1,0,1,0),1)]]:
        v=[F(0)]*144
        for c,value in entries:v[COLS.index(c)]=F(value)
        vectors.append(v)
    return vectors


def spurious():
    vectors=[]
    for component,power in itertools.product(range(3),range(3)):
        v=[F(0)]*144
        for s,t in itertools.product([1,3],[1,3]):
            v[COLS.index((component,power,s,t))]=F(-3 if s==1 else 5,2)*F(-3 if t==1 else 5,2)
        vectors.append(v)
    return vectors


def zeros(A,vectors):
    return all(sum((a*b for a,b in zip(row,v)),F(0))==0 for v in vectors for row in A)


def save(path,A,vectors):
    pairs=lambda rows:(np.array([[v.numerator for v in row] for row in rows],dtype=np.int64),np.array([[v.denominator for v in row] for row in rows],dtype=np.int64))
    n,d=pairs(A);vn,vd=pairs(vectors)
    np.savez_compressed(path,numerator=n,denominator=d,kernel_vectors_numerator=vn,kernel_vectors_denominator=vd)


def main():
    target=ROOT/'docs/verification/gsp-q3-flat-rank-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();out=ROOT/'build/analysis'/('gsp-flat-certificate-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    primes=[1000003,1000033,1000037]
    prime_checks={p:all(p%d for d in range(2,int(p**.5)+1)) for p in primes}
    assert all(prime_checks.values())
    records=[]
    for reduced in [False,True]:
        A=matrix(stress_rows(reduced));null=rigid()+(spurious() if reduced else [])
        exact_kernel=zeros(A,null)
        ranks={p:rank(modular(A,p),p) for p in primes}
        independent={p:rank(modular(null,p).T,p) for p in primes}
        upper=144-len(null)
        checks=dict(exact_kernel=exact_kernel,independent=all(v==len(null) for v in independent.values()),rank_matches_upper=all(v==upper for v in ranks.values()))
        artifact=out/('reduced.npz' if reduced else 'selected.npz');save(artifact,A,null)
        records.append(dict(space='negative_control_Q222' if reduced else 'selected_203',rows=len(A),columns=144,modular_ranks=ranks,null_directions=len(null),null_independence_ranks=independent,exact_kernel=exact_kernel,proven_rational_rank=upper if all(checks.values()) else None,checks=checks,artifact=artifact.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(artifact.read_bytes()).hexdigest()))
        print(records[-1],flush=True)
    sources=[Path(__file__),ROOT/'docs/verification/gsp-q3-flat-rank-protocol.md',ROOT/'docs/verification/gsp-q3-flat-rank-preparation.json']
    target.write_text(json.dumps(dict(records=records,prime_checks=prime_checks,all_checks_pass=all(all(r['checks'].values()) for r in records),scope='affine_positive_thickness_only',curved_uniform_accuracy_proven=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')


if __name__=='__main__':main()
