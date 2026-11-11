"""Exact rank certificate and Decimal coefficient-space elimination.

FE quadrature/Gram/constraints remain float64. No singular directions are dropped.
"""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from decimal import Decimal,localcontext
from fractions import Fraction as F
import hashlib
import numpy as np
import moment_enriched_exact_isometry as p
DIR=[(i,j) for i in range(6) for j in range(6) if (i,j)!=(5,5)]
def certify(columns,e,x,y,k):
 def combine(terms):
  out={}
  for col,factor in terms:
   for key,value in columns[col].items():
    new=out.get(key,F(0))+factor*value
    if new:out[key]=new
    else:out.pop(key,None)
  return not out
 def mterms(vector,phi):
  # vector must lie in the fixed center tangent plane.
  a,b=e;det=a[0]*b[1]-a[1]*b[0];alpha=(vector[0]*b[1]-vector[1]*b[0])/det;beta=(a[0]*vector[1]-a[1]*vector[0])/det
  assert alpha*a[2]+beta*b[2]==vector[2]
  return [(36*axis+6*i+j,v*factor) for axis,factor in enumerate([alpha,beta]) for (i,j),v in phi.items() if v*factor]
 xyz=[x,y,{key:k*value/2 for key,value in p.mul(x,x).items()}];d=[{key:-k*value for key,value in x.items()},{},{(0,0):F(1)}]
 constants=[]
 for cart,component,factor in [(0,5,F(1,2)),(1,1,F(1)),(2,3,F(1,2))]:
  terms=mterms([F(0),F(1),F(0)],xyz[cart]);terms += [(72+3*DIR.index(key)+1,value) for key,value in d[cart].items()];terms.append((177+component,-factor));constants.append(combine(terms))
 # The three independent fixed-plane rigid fields are in the membrane kernel.
 translation=[combine([(0,F(1))]),combine([(36,F(1))])]
 a,b=e;n=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
 rotation=[p.add({key:n[(axis+1)%3]*v for key,v in xyz[(axis+2)%3].items()},{key:-n[(axis+2)%3]*v for key,v in xyz[(axis+1)%3].items()}) for axis in range(3)]
 aa=sum(v*v for v in a);ab=sum(v*w for v,w in zip(a,b));bb=sum(v*v for v in b);det=aa*bb-ab*ab
 def dot(v):
  out={}
  for z,poly in zip(v,rotation):out=p.add(out,poly,z)
  return out
 da,db=dot(a),dot(b);phi=[{key:v/det for key,v in p.add({key:bb*v for key,v in da.items()},{key:-ab*v for key,v in db.items()}).items()},{key:v/det for key,v in p.add({key:aa*v for key,v in db.items()},{key:-ab*v for key,v in da.items()}).items()}]
 translation.append(combine([(36*axis+6*i+j,value) for axis,poly in enumerate(phi) for (i,j),value in poly.items()]))
 if not k:
  constants.extend([combine(mterms([F(1),F(0),F(0)],x)+[(177,F(-1))]),combine([(72+3*DIR.index((0,0)),F(2)),(181,F(-1))]),combine([(72+3*DIR.index((0,0))+2,F(1)),(179,F(-1))])])
 assert all(constants+translation),'exact_stress_rank_upper_bound_certificate'
 return dict(mem_upper_rank=69,joint_upper_rank=174,upper_rank=174 if not k else 177,three_rigid_identities=translation,constant_stress_identities=constants,pointwise_membrane_director_orthogonality=True)
def modular(columns,prime):
 keys=sorted(set().union(*(c.keys() for c in columns)));index={key:i for i,key in enumerate(keys)};basis=[];pivots=[];ranks=[];selected=[];cache={}
 for j,c in enumerate(columns):
  v=np.zeros(len(keys),dtype=np.int64)
  for key,value in c.items():
   den=value.denominator
   if den not in cache:cache[den]=pow(den%prime,-1,prime)
   v[index[key]]=value.numerator%prime*cache[den]%prime
  for pivot,b in zip(pivots,basis):
   factor=int(v[pivot])
   if factor:v=(v-factor*b)%prime
  nz=np.flatnonzero(v)
  if len(nz):
   pivot=int(nz[0]);v=v*pow(int(v[pivot]),-1,prime)%prime;basis.append(v);pivots.append(pivot);selected.append(j)
  ranks.append(len(basis))
 return ranks,selected
def construct(columns,proof,precision):
 primes=[1000000007,1000000009];mods=[modular(columns,prime) for prime in primes];ranks,selected=mods[0]
 if mods[1][0]!=ranks or (ranks[71],ranks[176],ranks[-1])!=(69,174,proof['upper_rank']):raise ValueError('modular_rank_and_exact_upper_bound_disagree')
 bases=[];pivots=[]
 with localcontext() as ctx:
  ctx.prec=precision
  for j in selected:
   v={key:Decimal(value.numerator)/Decimal(value.denominator) for key,value in columns[j].items()}
   for pivot,b in zip(pivots,bases):
    factor=v.get(pivot,Decimal(0))
    if factor:
     for key,value in b.items():
      new=v.get(key,Decimal(0))-factor*value
      if new:v[key]=new
      else:v.pop(key,None)
     v.pop(pivot,None)
   pivot=max(v,key=lambda key:abs(v[key]));scale=v[pivot];v={key:value/scale for key,value in v.items()};v[pivot]=Decimal(1);bases.append(v);pivots.append(pivot)
 encoded='\n'.join(';'.join(str(key)+':'+str(value) for key,value in sorted(c.items())) for c in columns).encode()
 return bases,ranks,dict(**proof,modular_primes=primes,exact_lower_rank=ranks[-1],rational_space_sha256=hashlib.sha256(encoded).hexdigest(),basis_preprocessing_decimal_digits=precision)
