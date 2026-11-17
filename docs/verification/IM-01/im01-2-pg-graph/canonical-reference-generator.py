from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
u=lambda n:n.to_bytes(8,'big')
def text(s):b=s.encode();return u(len(b))+b
def rational(x):
 sign=bytes([x<0]);p=abs(x.numerator);q=x.denominator
 a=p.to_bytes((p.bit_length()+7)//8,'big');b=q.to_bytes((q.bit_length()+7)//8,'big')
 return sign+u(len(a))+a+u(len(b))+b
terms=[(0,0,F(1,11)),(0,4,F(-2,5)),(2,1,F(2,3)),(4,0,F(1))]
rows=[]
for kind,coefficients in [('flat',[]),('quartic',terms)]:
 for n in [1,2,4]:
  data=text('uvlas.pg.graph.v1')+u(1)+b''.join(rational(x) for x in [F(1,7),F(4,5),F(3,20)])
  data+=u(len(coefficients))+b''.join(u(a)+u(b)+rational(c) for a,b,c in coefficients)
  data+=u(n)+u((2*n)**2)
  for a in range(2*n):
   for b in range(2*n):data+=b''.join(u(x) for x in [a,b,a*(2*n)+b,a,a+1,b,b+1,2*n,4*n])
  data+=text('anchor-s0-first-cell-v1')+u(0)+rational(F(-1))+text('UnitSquare-v1')+b''.join(u(x) for x in [0,1,0,1])
  rows.append({'kind':kind,'n':n,'sha256':sha256(data).hexdigest(),'canonical_byte_count':len(data),'canonical_bytes_hex':data.hex()})
result={'schema':'uvlas.pg-graph-reference/1','id':'pg-graph-canonical-v1','provenance':'Independent Python standard-library Fraction/hashlib encoder of the specified canonical stream. All coefficient fractions and cell formulas are literal; no product code, product output or R&D rewrite is used.','graph':{'x0':['1','7'],'w':['4','5'],'a':['3','20'],'quartic_terms':[[a,b,str(c.numerator),str(c.denominator)] for a,b,c in terms]},'identities':rows}
Path('tests/structure/fixtures/pg-graph-v1.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
