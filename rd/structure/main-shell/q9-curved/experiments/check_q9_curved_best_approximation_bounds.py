"""Convex lower bounds and feasible upper bounds, with rounding envelopes."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,time
import numpy as np
import check_q9_curved_bounded_lift as core
ROOT=core.ROOT;F=core.F;UNIT=np.finfo(float).eps/2
def gamma(n):return n*UNIT/(1-n*UNIT)
def tree(x):
    while x.shape[1]>1:
        length=x.shape[1];even=length//2*2;y=x[:,:even:2]+x[:,1:even:2]
        if length%2:y=np.column_stack((y,x[:,-1]))
        x=y
    return x[:,0]
def dot_rows(A,x):
    values=[];errors=[];depth=int(np.ceil(np.log2(max(1,len(x)))))+2
    for start in range(0,len(A),32):
        products=A[start:start+32]*x[None,:];value=tree(products.copy());absolute=tree(abs(products));error=gamma(depth)*absolute/(1-gamma(depth));values.extend(value);errors.extend(np.nextafter(error,np.inf))
    return np.array(values),np.array(errors)
def upper_sum(x):return float(np.sum(abs(x))/(1-gamma(len(x)+2)))
def norms(x,error=None):
    value=float(np.linalg.norm(x));lo=value*(1-gamma(len(x)+3));hi=value/(1-gamma(len(x)+3))
    if error is not None:_,e=norms(error);lo=max(0.,lo-e);hi+=e
    return float(np.nextafter(lo,0.)),float(np.nextafter(hi,np.inf))
def dot(a,b):
    v,e=dot_rows(np.asarray(a)[None,:],np.asarray(b));return float(v[0]),float(e[0])
def positive_lower(G):
    _,Q=np.linalg.eigh(G);n=len(G);g=gamma(n+2);T=G@Q;errT=g*(abs(G)@abs(Q))/(1-g);B=Q.T@T;errB=(g*(abs(Q).T@abs(T))+abs(Q).T@errT)/(1-g)**2;O=Q.T@Q;errO=g*(abs(Q).T@abs(Q))/(1-g)
    radii=np.sum(abs(B),axis=1)-abs(B.diagonal());low=np.min(B.diagonal()-radii-np.sum(errB,axis=1)-gamma(n+5)*np.sum(abs(B),axis=1));orthlow=np.min(O.diagonal()-(np.sum(abs(O),axis=1)-abs(O.diagonal()))-np.sum(errO,axis=1));orthup=np.max(np.sum(abs(O),axis=1)+np.sum(errO,axis=1))/(1-gamma(n+3))
    if low<=0 or orthlow<=0:raise ValueError('uncertified_positive_gram')
    return float(np.nextafter(.5*low/orthup,0.))
def shear(A,W,p,r):
    ap,ae=dot_rows(A,p);wr,we=dot_rows(W,r);value=ap+wr;error=ae+we+gamma(2)*(abs(ap)+abs(wr));return value,error
def metric_square(M,p):
    P=p.reshape(len(M),-1);value=0.;error=0.
    for k in range(P.shape[1]):
        mp,me=dot_rows(M,P[:,k]);v,e=dot(P[:,k],mp);value+=v;error+=e+upper_sum(abs(P[:,k])*me)
    error+=gamma(P.shape[1]+3)*(abs(value)+error)
    return float(value),float(np.nextafter(error,np.inf))
def objective(G,b,target,r):
    gr,ge=dot_rows(G,r);quad,qe=dot(r,gr);linear,le=dot(r,b);value=(quad-2*linear+target)/2;error=(qe+upper_sum(abs(r)*ge)+2*le)/2+gamma(10)*(abs(quad)+2*abs(linear)+abs(target))
    return float(np.nextafter(value-error,-np.inf)),float(np.nextafter(value+error,np.inf))
def certify(a,r,p,lam,eps,bound,k):
    G,M,A,W,b,target=a['G'],a['M'],a['A'],a['W'],a['b'][:,k],float(a['norm'][k]);gmin=positive_lower(G);hmin=positive_lower(M);hroot=float(np.nextafter(np.sqrt(hmin),0.));res,re=shear(A,W,p,r);rlow,rup=norms(res,re);gr0,gr0e=dot_rows(G,r);wr,wre=dot_rows(W.T,res);we,wee=dot_rows(abs(W.T),re);gr=gr0-b+lam*wr;gre=gr0e+lam*(wre+we+wee)+gamma(6)*(abs(gr0)+abs(b)+abs(lam*wr));ap,ape=dot_rows(A.T,res);ae,aee=dot_rows(abs(A.T),re);gp=lam*ap;gpe=lam*(ape+ae+aee)+gamma(3)*abs(gp)
    lower_f,_=objective(G,b,target,r);lagterm=lam/2*(rlow**2-eps**2);lagmargin=gamma(10)*abs(lam)*(rup**2+eps**2);laglower=lower_f+lagterm-lagmargin
    z=np.linalg.solve(M,gp.reshape(len(M),-1));q=np.zeros_like(z);qe=np.zeros_like(z)
    for j in range(z.shape[1]):
        mz,mze=dot_rows(M,z[:,j]);q[:,j]=gp.reshape(z.shape)[:,j]-mz;qe[:,j]=gpe.reshape(z.shape)[:,j]+mze+gamma(3)*(abs(gp.reshape(z.shape)[:,j])+abs(mz))
    zn,ze=metric_square(M,z.reshape(-1));dual_upper=np.sqrt(max(0.,zn+ze))+norms(q.reshape(-1),qe.reshape(-1))[1]/hroot;dual_upper=float(np.nextafter(dual_upper*(1+gamma(10)),np.inf));gp_dot,e=dot(gp,p);gp_dot_up=gp_dot+e+upper_sum(abs(p)*gpe);gnorm=norms(gr,gre)[1];penalty=gnorm**2/(2*gmin);terms=[laglower,-penalty,-bound*dual_upper,-gp_dot_up];lb=sum(terms)-gamma(32)*sum(abs(t) for t in terms);lb=float(np.nextafter(lb,-np.inf))
    rr=r.copy();pp=p.copy();scale=1.
    for _ in range(6):
        v9,e9=shear(A,W,pp,rr);v11,e11=shear(a.get('A11',A),a.get('W11',W),pp,rr);s9=norms(v9,e9)[1];s11=norms(v11,e11)[1];mn,me=metric_square(M,pp);pn=np.sqrt(max(0.,mn+me));factor=min(1.,eps/max(s9,s11,1e-300),bound/max(pn,1e-300))
        if s9<=eps and s11<=eps and pn<=bound:break
        factor*=1-1e-8;rr*=factor;pp*=factor;scale*=factor
    feasible=bool(s9<=eps and s11<=eps and pn<=bound)
    _,upper=objective(G,b,target,rr);_,upper11=objective(a.get('G11',G),a.get('b11',a['b'])[:,k],float(a.get('norm11',a['norm'])[k]),rr)
    result=dict(lower_objective=lb,upper_objective=upper,lower_relative_error=float(np.sqrt(max(0.,2*lb/target))),upper_relative_error=float(np.sqrt(max(0.,2*upper/target))),upper_relative_error11=float(np.sqrt(max(0.,2*upper11/float(a.get('norm11',a['norm'])[k])))),position_gram_lower=gmin,h1_gram_lower=hmin,scale=scale,shear_upper9=s9,shear_upper11=s11,h1_upper=pn,feasible=feasible,objective_order=bool(lb<=upper),gradient_position_upper=gnorm,gradient_direction_uncertainty=norms(gpe)[1],dual_norm_upper=dual_upper)
    return result
def synthetic():
    cases=[]
    for title,Aval,B,r,p,lam,opt in [('shear_only',0.,1.,.2,0.,4.,.8),('two_active',1.,.1,.3,-.1,3.5,.7),('nonoptimal_point',1.,.1,.1,0.,3.5,.7),('inactive',1.,2.,1.,-1.,0.,0.)]:
        a=dict(A=np.array([[Aval]]),W=np.ones((1,1)),M=np.ones((1,1)),G=np.ones((1,1)),b=np.ones((1,1)),norm=np.ones(1));record=certify(a,np.array([r]),np.array([p]),lam,.2,B,0);record.update(name=title,known_error=opt,contains=bool(record['lower_relative_error']<=opt+1e-12 and record['upper_relative_error']>=opt-1e-12));cases.append(record)
    return cases
def main():
    target=F/'q9-curved-best-approximation-bounds-result.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();tests=synthetic();assert all(t['contains'] and t['feasible'] and t['objective_order'] for t in tests)
    parent=F/'q9-curved-best-approximation-result.json';data=json.loads(parent.read_text());assert len(data['cases'])==8 and data['stop_reason'] is None;records=[]
    for c in data['cases']:
        if c['geometry'] not in ['SADDLE','NORMAL-C'] or c['h']!=.001:continue
        assert core.digest(ROOT/c['artifact'])==c['artifact_sha256'];a=np.load(ROOT/c['artifact'])
        for i,record in enumerate(c['records']):
            if record['wave_index']!=0 or record['bound_factor']!=2 or record['relative_tolerance'] not in [1e-4,1e-6]:continue
            result=certify(a,a['r'][i],a['p'][i],record['lambda_shear'],record['epsilon'],record['bound'],0);result.update(geometry=c['geometry'],n=c['n'],h=c['h'],relative_tolerance=record['relative_tolerance'],epsilon=record['epsilon'],bound=record['bound'],parent_record=i);records.append(result);print(json.dumps(result),flush=True)
        if time.monotonic()-start>300 or core.memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
    result=dict(source_sha256=core.digest(__file__),parent_sha256=core.digest(parent),protocol_sha256=core.digest(F/'q9-curved-best-approximation-bounds-protocol.md'),dependency_sha256=core.digest(core.__file__),synthetic=tests,cases=records,seconds=time.monotonic()-start,peak_bytes=core.memory.memory_bytes(),checks_pass=bool(len(records)==8 and all(r['feasible'] and r['objective_order'] for r in records)),kkt_gate_replaced=False,product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['cases','synthetic']}))
    if not result['checks_pass']:raise SystemExit(1)
if __name__=='__main__':main()
