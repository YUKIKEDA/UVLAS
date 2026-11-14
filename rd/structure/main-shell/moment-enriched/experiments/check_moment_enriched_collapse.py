"""Exact polynomial identities, distorted-edge connectivity, positive L2 gap."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import moment_enriched_exact_isometry as p
import moment_enriched_ancf as f
ROOT=rd_paths.ROOT;k=F(2,5);one={(0,0):F(1)};x={(1,0):F(1)};y={(0,1):F(1)}
def scale(a,s):return {key:v*s for key,v in a.items() if v*s}
def sub(a,b):return p.add(a,b,F(-1))
def power(a,n):
 out=one
 for _ in range(n):out=p.mul(out,a)
 return out
def fields(a,b,c0=F(0)):
 ap=p.deriv(a,0);app=p.deriv(ap,0);appp=p.deriv(app,0);bp=p.deriv(b,0)
 c={(degree+1,0):-k*value/(degree+1) for (degree,j),value in p.mul(x,bp).items()};c=p.add(c,scale(one,c0))
 u=[p.add(p.mul(y,sub(p.mul(x,app),ap)),c),a,sub(b,scale(p.mul(y,app),1/k))]
 dd=[sub(scale(p.mul(y,appp),1/k),bp),p.add(sub(scale(p.mul(power(x,2),app),k),scale(p.mul(x,ap),k)),scale(app,1/k)),sub(p.mul(p.mul(x,y),appp),scale(p.mul(x,bp),k))]
 equations=[p.add(p.deriv(u[0],0),scale(p.mul(x,p.deriv(u[2],0)),k)),p.deriv(u[1],1),p.add(p.add(p.deriv(u[0],1),scale(p.mul(x,p.deriv(u[2],1)),k)),p.deriv(u[1],0)),
 p.add(p.add(dd[0],scale(p.mul(x,dd[2]),k)),sub(p.deriv(u[2],0),scale(p.mul(x,p.deriv(u[0],0)),k))),p.add(dd[1],sub(p.deriv(u[2],1),scale(p.mul(x,p.deriv(u[0],1)),k))),sub(dd[2],scale(p.mul(x,dd[0]),k))]
 return u,dd,all(not eq for eq in equations)
def compose(poly,xx,yy):
 out={}
 for (a,b),v in poly.items():out=p.add(out,p.mul(power(xx,a),power(yy,b)),v)
 return out
def solve(A,b):
 rows=[list(row)+[rhs] for row,rhs in zip(A,b)];n=len(b)
 for j in range(n):
  pivot=next(i for i in range(j,n) if rows[i][j]);rows[j],rows[pivot]=rows[pivot],rows[j];z=rows[j][j];rows[j]=[v/z for v in rows[j]]
  for i in range(n):
   if i!=j:
    z=rows[i][j];rows[i]=[a-z*b for a,b in zip(rows[i],rows[j])]
 return [r[-1] for r in rows]
def main():
 dest=ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-collapse-result.json';assert not dest.exists();identities=[];graphs=[];spaces=[]
 for degree in range(5):
  for family in ['a','b']:
   u,d,ok=fields(power(x,degree) if family=='a' else {},power(x,degree) if family=='b' else {});identities.append(dict(family=family,degree=degree,all_six_identities_zero=ok))
 identities.append(dict(family='c',degree=0,all_six_identities_zero=fields({}, {},F(1))[2]))
 for n in [2,4]:
  grid=f.m.mesh(.4,n,.2);edges=[];cellchecks=[]
  for cell,ids in enumerate(grid['cells']):
   q=grid['ref'][ids].reshape(75,3);xx=p.bilinear(q[[0,4,20,24],0]);yy=p.bilinear(q[[0,4,20,24],1]);d=xx.get((1,1),F(0));within=True
   for degree in range(5):
    for family in ['a','b']:
     u,dd,ok=fields(power(x,degree) if family=='a' else {},power(x,degree) if family=='b' else {})
     for poly in u+dd:within &= all(i<=5 and j<=5 for i,j in compose(poly,xx,yy))
     for poly in dd:
      cp=compose(poly,xx,yy);within &= not cp.get((5,5),F(0))
   high=compose(fields(power(x,5),{})[1][1],xx,yy).get((5,5),F(0));corner=high==15*k*d**5 and high!=0
   cellchecks.append(dict(cell=cell,bilinear_cross_nonzero=d!=0,eleven_fields_in_spaces=bool(within),a5_director_corner_matches=corner))
   j,i=divmod(cell,n)
   for axis,neighbor,endpoints in [(0,cell+1,(4,24)),(1,cell+n,(20,24))]:
    if (axis==0 and i==n-1) or (axis==1 and j==n-1):continue
    slope=(F.from_float(float(q[endpoints[1],0]))-F.from_float(float(q[endpoints[0],0])))/2
    edges.append((cell,neighbor,slope!=0))
  reached={0}
  for _ in range(n*n):
   for a,b,nonconstant in edges:
    if nonconstant and (a in reached or b in reached):reached.update([a,b])
  graphs.append(dict(n=n,edges=len(edges),all_edges_x_nonconstant=all(e[2] for e in edges),reached=len(reached),connected=len(reached)==n*n));spaces.append(dict(n=n,cells=cellchecks))
 G=[[F(1,i+j+1) for j in range(1,5)] for i in range(1,5)];v=[F(1,i+6) for i in range(1,5)];coef=solve(G,v);gap=F(1,11)-sum(a*b for a,b in zip(v,coef));relative2=11*gap
 passed=all(r['all_six_identities_zero'] for r in identities) and all(g['connected'] and g['all_edges_x_nonconstant'] for g in graphs) and all(all(c['bilinear_cross_nonzero'] and c['eleven_fields_in_spaces'] and c['a5_director_corner_matches'] for c in r['cells']) for r in spaces) and gap>0
 paths=[Path(__file__),ROOT/'rd/structure/main-shell/moment-enriched/moment-enriched-thin-collapse.md',ROOT/'rd/structure/main-shell/moment-enriched/experiments/moment_enriched_exact_isometry.py',ROOT/'rd/structure/main-shell/moment-enriched/records/moment-enriched-q0-space-v2-result.json']
 dest.write_text(json.dumps(dict(identities=identities,graphs=graphs,spaces=spaces,xy_L2_w5_gap_squared=str(gap),xy_L2_w5_relative_gap_squared=str(relative2),xy_L2_w5_relative_gap=float(relative2)**.5,passed=bool(passed),geometry_conditional=True,finite_h_lower_bound=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2),encoding='utf-8');print('exact collapse',passed,'graphs',graphs,'L2 gap2',gap,'relative',(float(relative2)**.5),flush=True)
if __name__=='__main__':main()
