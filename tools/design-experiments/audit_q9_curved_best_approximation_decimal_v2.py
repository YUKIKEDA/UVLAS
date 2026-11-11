"""80-digit sparse-row audit of one rounded-data convex certificate."""
from decimal import Decimal,localcontext
from pathlib import Path
import json,time
import numpy as np
import check_q9_curved_best_approximation_bounds as bounds
core=bounds.core;ROOT=core.ROOT;F=core.F
def dec(x):return Decimal.from_float(float(x))
def ldl(matrix):
    n=len(matrix);L=[[Decimal(int(i==j)) for j in range(n)] for i in range(n)];D=[]
    for j in range(n):
        value=matrix[j][j]-sum((L[j][k]**2*D[k] for k in range(j)),Decimal(0));assert value>0;D.append(value)
        for i in range(j+1,n):L[i][j]=(matrix[i][j]-sum((L[i][k]*L[j][k]*D[k] for k in range(j)),Decimal(0)))/value
    return L,D
def solve(factor,b):
    L,D=factor;n=len(D);y=[]
    for i in range(n):y.append(b[i]-sum((L[i][j]*y[j] for j in range(i)),Decimal(0)))
    y=[v/d for v,d in zip(y,D)];x=[Decimal(0)]*n
    for i in range(n-1,-1,-1):x[i]=y[i]-sum((L[j][i]*x[j] for j in range(i+1,n)),Decimal(0))
    return x
def sparse(A,W):
    rows=[]
    for a,w in zip(A,W):
        rows.append([(int(i),dec(a[i])) for i in np.flatnonzero(a)]+[(int(len(a)+i),dec(w[i])) for i in np.flatnonzero(w)])
    return rows
def apply(rows,x,gradient=False):
    y=[];g=[Decimal(0)]*len(x)
    for row in rows:
        value=sum((v*x[i] for i,v in row),Decimal(0));y.append(value)
        if gradient:
            for i,v in row:g[i]+=v*value
    return y,g
def quadratic(G,b,target,r):return (sum((r[i]*sum((G[i][j]*r[j] for j in range(len(r))),Decimal(0)) for i in range(len(r))),Decimal(0))-2*sum((v*w for v,w in zip(r,b)),Decimal(0))+target)/2
def main():
    target=F/'q9-curved-best-approximation-decimal-v2-audit.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-curved-best-approximation-bounds-result.json';br=json.loads(parent.read_text());assert br['checks_pass'];pr=F/'q9-curved-best-approximation-result.json';r=json.loads(pr.read_text());c=next(c for c in r['cases'] if c['geometry']=='SADDLE' and c['n']==2 and c['h']==.001);record=next(c for c in br['cases'] if c['geometry']=='SADDLE' and c['n']==2 and c['relative_tolerance']==1e-6);assert core.digest(ROOT/c['artifact'])==c['artifact_sha256'];a=np.load(ROOT/c['artifact']);index=record['parent_record'];fit=c['records'][index];p0,r0=a['p'][index],a['r'][index];start=time.monotonic()
    with localcontext() as ctx:
        ctx.prec=80;rows=sparse(a['A'],a['W']);x=[dec(v) for v in np.r_[p0,r0]];y,g=apply(rows,x,True);lam=dec(fit['lambda_shear']);eps=dec(fit['epsilon']);B=dec(fit['bound']);np_=len(p0);gp=[lam*v for v in g[:np_]];G=[[dec(v) for v in row] for row in a['G']];M=[[dec(v) for v in row] for row in a['M']];b=[dec(v) for v in a['b'][:,0]];norm=dec(a['norm'][0]);rd=x[np_:];pd=x[:np_];gr=[sum((G[i][j]*rd[j] for j in range(len(rd))),Decimal(0))-b[i]+lam*g[np_+i] for i in range(len(rd))]
        yh,ye=bounds.shear(a['A'],a['W'],p0,r0);forward_ratio=max(float(abs(value-dec(v))/dec(e)) if e else float(value!=dec(v)) for value,v,e in zip(y,yh,ye));assert forward_ratio<=1
        ap,ape=bounds.dot_rows(a['A'].T,yh);ae,aee=bounds.dot_rows(abs(a['A'].T),ye);gph=fit['lambda_shear']*ap;gpe=fit['lambda_shear']*(ape+ae+aee)+bounds.gamma(3)*abs(gph);gradient_ratio=max(float(abs(value-dec(v))/dec(e)) if e else float(value!=dec(v)) for value,v,e in zip(gp,gph,gpe));assert gradient_ratio<=1
        gr0,gr0e=bounds.dot_rows(a['G'],r0);wr,wre=bounds.dot_rows(a['W'].T,yh);we,wee=bounds.dot_rows(abs(a['W'].T),ye);grh=gr0-a['b'][:,0]+fit['lambda_shear']*wr;gre=gr0e+fit['lambda_shear']*(wre+we+wee)+bounds.gamma(6)*(abs(gr0)+abs(a['b'][:,0])+abs(fit['lambda_shear']*wr));position_ratio=max(float(abs(value-dec(v))/dec(e)) if e else float(value!=dec(v)) for value,v,e in zip(gr,grh,gre));assert position_ratio<=1
        gm=dec(record['position_gram_lower']);hm=dec(record['h1_gram_lower']);ldl([[v-(gm if i==j else 0) for j,v in enumerate(row)] for i,row in enumerate(G)]);ldl([[v-(hm if i==j else 0) for j,v in enumerate(row)] for i,row in enumerate(M)]);factor=ldl(M);dual_square=Decimal(0)
        for k in range(6):
            rhs=gp[k::6];z=solve(factor,rhs);dual_square+=sum((v*w for v,w in zip(rhs,z)),Decimal(0))
        objective=quadratic(G,b,norm,rd);lag=objective+lam/2*(sum((v*v for v in y),Decimal(0))-eps*eps);lower=lag-sum((v*v for v in gr),Decimal(0))/(2*gm)-B*dual_square.sqrt()-sum((v*w for v,w in zip(gp,pd)),Decimal(0));lower_check=dec(record['lower_objective'])<=lower
        rr=np.array(r0*record['scale']);pp=np.array(p0*record['scale']);xx=[dec(v) for v in np.r_[pp,rr]];yy,_=apply(rows,xx);other,_=apply(sparse(a['A11'],a['W11']),xx);h1_square=sum((sum((xx[6*i+k]*M[i][j]*xx[6*j+k] for i in range(len(M)) for j in range(len(M))),Decimal(0)) for k in range(6)),Decimal(0));feasible=sum((v*v for v in yy),Decimal(0))<=eps**2 and sum((v*v for v in other),Decimal(0))<=eps**2 and h1_square<=B**2;upper=quadratic(G,b,norm,xx[np_:]);upper_check=upper<=dec(record['upper_objective'])
        result=dict(source_sha256=core.digest(__file__),parent_sha256=core.digest(parent),artifact_sha256=c['artifact_sha256'],protocol_sha256=core.digest(F/'q9-curved-best-approximation-decimal-protocol.md'),dependency_sha256=core.digest(bounds.__file__),precision=80,forward_interval_maximum_ratio=forward_ratio,gradient_interval_maximum_ratio=gradient_ratio,position_interval_maximum_ratio=position_ratio,decimal_lower_expression=str(lower),reported_lower_objective=record['lower_objective'],decimal_feasible_objective=str(upper),reported_upper_objective=record['upper_objective'],checks=dict(forward_interval=forward_ratio<=1,gradient_interval=gradient_ratio<=1,position_interval=position_ratio<=1,positive_gram=True,lower=bool(lower_check),feasible=bool(feasible),upper=bool(upper_check)),seconds=time.monotonic()-start,peak_bytes=core.memory.memory_bytes())
    result['all_criteria_pass']=bool(all(result['checks'].values()) and result['seconds']<=300 and result['peak_bytes']<=1024**3);target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
