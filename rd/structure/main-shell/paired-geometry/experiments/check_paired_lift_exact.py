"""Rational polynomial identities for both lifting directions and full charts."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time,traceback
import general_paired_polynomials as exact
from paired_reference_admissibility import derivative,combine
ROOT=rd_paths.ROOT
p=exact.p
def dot(a,b):return p.add(*[p.mul(x,y) for x,y in zip(a,b)])

def main():
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-lift-exact-result.json';assert not dest.exists();start=time.monotonic();records=[];stops=[]
    # All monomials through degree four, including terms absent from the default family.
    coefficients={(i,j):F((-1)**(i+j)*(i+2*j+1),100*(i+j+1)) for i in range(5) for j in range(5-i)}
    try:
        for width,skew in [(F(1),F(0)),(F(4,5),F(3,20))]:
            count=0;negative=0
            for cx in range(2):
                for cy in range(2):
                    x,y,f=exact.graph_polynomial(1,cx,cy,coefficients)
                    fx=p.scale(p.derivative(f,0,1),1/width);fy=p.add(p.derivative(f,1,1),p.scale(fx,-skew))
                    one={(0,0):F(1)};zero={};g=p.add(one,p.mul(fx,fx),p.mul(fy,fy))
                    A=[{(0,0):width},zero,p.derivative(f,0,1)];B=[{(0,0):skew},one,p.derivative(f,1,1)];normal=[p.scale(fx,-1),p.scale(fy,-1),one]
                    for axis in [0,1]:
                        for origin in range(2):
                            for k in range(8):
                                v=exact.lifting(1,cx,cy,coefficients,origin,k,axis,width=width,shear=skew)
                                ds=[p.derivative(a,0,1) for a in v];dt=[p.derivative(a,1,1) for a in v]
                                D=exact.director_numerator(v,fx,fy,1,width,skew)
                                identities=[dot(A,ds) if axis==0 else dot(B,dt),p.add(dot(A,dt),dot(B,ds)),dot(normal,D),
                                    p.add(dot(A,D),p.mul(g,dot(normal,ds))),p.add(dot(B,D),p.mul(g,dot(normal,dt)))]
                                assert all(not a for a in identities)
                                # v_{st}=v_{ts}, exactly in the coefficient algebra.
                                assert all(p.derivative(a,1,1)==p.derivative(b,0,1) for a,b in zip(ds,dt))
                                if dot(normal,ds) or dot(normal,dt):negative+=1
                                count+=1
            assert negative>0
            records.append(dict(width=str(width),skew=str(skew),fields=count,identities_per_field=8,
                coefficient_residual_exactly_zero=True,zero_director_negative_controls=negative))
    except Exception as exc:stops.append(dict(reason=repr(exc),traceback=traceback.format_exc()))
    paths=[Path(__file__),Path(exact.__file__),Path(p.__file__),ROOT/'rd/structure/main-shell/paired-geometry/experiments/paired_reference_admissibility.py']
    result=dict(passed=len(records)==2 and not stops,records=records,stops=stops,seconds=time.monotonic()-start,
        arbitrary_degree_derivation='linearity and integration by parts; numeric identities above cover every monomial through total degree four',
        sources={s.relative_to(ROOT).as_posix():hashlib.sha256(s.read_bytes()).hexdigest() for s in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
