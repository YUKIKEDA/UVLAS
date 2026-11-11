"""MMS-1 rational numerator stress space for the analytic cylinder benchmark."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib,itertools,time
import numpy as np
from scipy.linalg import qr,svd
import moment_enriched_exact_isometry as p
import moment_enriched_ancf as f
from mixed_mean_stress_basis import certify,construct
PAIRS=[(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]
def sc(a,v):return {k:w*v for k,w in a.items() if w*v}
def sub(a,b):return p.add(a,b,F(-1))
def outer(a,b):return [sc(p.add(p.mul(a[i],b[j]),p.mul(a[j],b[i])),F(1,2)) for i,j in PAIRS]
def vecdot(a,b):
 out={}
 for x,y in zip(a,b):out=p.add(out,p.mul(x,y))
 return out
def mul(a,b):return p.mul(a,b)
def monomial_legendre(field):
 out={}
 for component,poly in enumerate(field):
  for (a,b),value in poly.items():
   for i,vi in f.c.powers(a):
    for j,vj in f.c.powers(b):
     key=(i,j,component);out[key]=out.get(key,F(0))+value*vi*vj
 return {k:v for k,v in out.items() if v}
def exact_basis(columns):
 import math
 bases=[];pivots=[];ranks=[]
 def primitive(v):
  gcd=0
  for value in v.values():gcd=math.gcd(gcd,abs(value))
  return {key:value//gcd for key,value in v.items()} if gcd>1 else v
 for col in columns:
  denominator=1
  for value in col.values():denominator=math.lcm(denominator,value.denominator)
  v=primitive({key:value.numerator*(denominator//value.denominator) for key,value in col.items()})
  for pivot,base in zip(pivots,bases):
   w=v.get(pivot,0)
   if not w:continue
   u=base[pivot];gcd=math.gcd(u,w);left=u//gcd;right=w//gcd
   v={key:left*value for key,value in v.items()}
   for key,value in base.items():
    new=v.get(key,0)-right*value
    if new:v[key]=new
    else:v.pop(key,None)
   v=primitive(v)
  if v:
   pivot=max(v,key=lambda key:abs(v[key]))
   if v[pivot]<0:v={key:-value for key,value in v.items()}
   pivots.append(pivot);bases.append(v)
  ranks.append(len(bases))
 return [{key:F(value,base[pivot]) for key,value in base.items()} for base,pivot in zip(bases,pivots)],ranks
def plane_basis(points,ids):
 va=np.polynomial.legendre.legvander(points[:,0],max(i for i,j in ids));vb=np.polynomial.legendre.legvander(points[:,1],max(j for i,j in ids))
 return np.column_stack([va[:,i]*vb[:,j]*np.sqrt((2*i+1)*(2*j+1)) for i,j in ids])
@lru_cache(maxsize=96)
def coefficients(qbytes,kappa,independent=False):
 start=time.monotonic();q=np.frombuffer(qbytes,dtype=np.float64).reshape(75,3);x=p.bilinear(q[[0,4,20,24],0]);y=p.bilinear(q[[0,4,20,24],1]);k=F(str(kappa));one={(0,0):F(1)};zero={}
 tx=[one,zero,sc(x,k)];d=[sc(x,-k),zero,one];g=p.add(one,sc(mul(x,x),k*k));xx=p.deriv(x,0);xy=p.deriv(x,1);yx=p.deriv(y,0);yy=p.deriv(y,1);jac=sub(mul(xx,yy),mul(xy,yx));den=mul(mul(g,g),jac)
 xc=x.get((0,0),F(0));directions=[[{(0,0):v} if v else {} for v in (xx.get((0,0),F(0)),yx.get((0,0),F(0)),k*xc*xx.get((0,0),F(0)))],[{(0,0):v} if v else {} for v in (xy.get((0,0),F(0)),yy.get((0,0),F(0)),k*xc*xy.get((0,0),F(0)))]]
 columns=[]
 for e in directions:
  pe=[sub(mul(g,a),mul(b,vecdot(d,e))) for a,b in zip(e,d)]
  for i,j in itertools.product(range(6),repeat=2):
   phi={(i,j):F(1)};a=sub(mul(yy,p.deriv(phi,0)),mul(yx,p.deriv(phi,1)));b=sub(mul(xx,p.deriv(phi,1)),mul(xy,p.deriv(phi,0)))
   grad=[mul(t,a) for t in tx];grad[1]=p.add(grad[1],mul(g,b));columns.append(monomial_legendre(outer(pe,grad)))
 for i,j in itertools.product(range(6),repeat=2):
  if (i,j)==(5,5):continue
  phi={(i,j):F(1)}
  for cart in range(3):
   a=[phi if axis==cart else {} for axis in range(3)];columns.append(monomial_legendre([mul(v,mul(g,jac)) for v in outer(a,d)]))
 for component in range(6):columns.append(monomial_legendre([den if a==component else {} for a in range(6)]))
 proof=certify(columns,[[poly.get((0,0),F(0)) for poly in e] for e in directions],x,y,k)
 bases,ranks,basis_certificate=construct(columns,proof,100 if independent else 80);assert ranks[71]==69 and ranks[176]==174,(ranks[71],ranks[176])
 maxi=max(key[0] for b in bases for key in b);maxj=max(key[1] for b in bases for key in b);ids=list(itertools.product(range(maxi+1),range(maxj+1)));lookup={key:i for i,key in enumerate(ids)}
 B=np.zeros((6*len(ids),len(bases)))
 for j,b in enumerate(bases):
  for (a,c,component),v in b.items():B[6*lookup[a,c]+component,j]=float(v)/np.sqrt((2*a+1)*(2*c+1))*(1 if component<3 else np.sqrt(2))
 norms=np.linalg.norm(B,axis=0);B/=norms;condition=float(np.linalg.cond(B));Q,T=qr(B,mode='economic',check_finite=False)
 if condition>1e8:raise ValueError('exact_stress_basis_conversion_condition:'+str(condition))
 raw=np.zeros((6*len(ids),len(columns)))
 for j,col in enumerate(columns):
  for (a,b,component),value in col.items():raw[6*lookup[a,b]+component,j]=float(value)/np.sqrt((2*a+1)*(2*b+1))*(1 if component<3 else np.sqrt(2))
 projection_error=float(np.linalg.norm(raw-Q@(Q.T@raw))/max(1.,np.linalg.norm(raw)))
 if projection_error>1e-12:raise ValueError('coefficient_space_containment:'+str(projection_error))
 encoded='\n'.join(';'.join(str(key)+':'+str(v) for key,v in sorted(b.items())) for b in bases).encode()
 return Q.reshape(len(ids),6,-1),ids,den,dict(mem_rank=ranks[71],normal_shear_rank=105,combined_rank=ranks[176],stress_rank=ranks[-1],constant_completion=ranks[-1]-ranks[176],coefficient_condition=condition,coefficient_basis_sha256=hashlib.sha256(encoded).hexdigest(),exact_space_certificate=basis_certificate,coefficient_projection_error=projection_error,seconds=time.monotonic()-start)
def evaluate(q,kappa,points,independent=False):
 coeff,ids,den,meta=coefficients(np.asarray(q,dtype=np.float64).tobytes(),kappa,independent);D=sum(float(v)*points[:,0]**i*points[:,1]**j for (i,j),v in den.items())
 if np.min(D)<=0:raise ValueError('stress_reference_chart_orientation')
 return np.einsum('pk,kjr->pjr',plane_basis(points,ids),coeff,optimize=True)/D[:,None,None],meta
