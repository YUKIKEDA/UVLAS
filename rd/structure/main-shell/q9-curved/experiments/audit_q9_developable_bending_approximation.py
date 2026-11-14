"""Independent gradients/isotropic tensor energies and a Decimal constrained LS solve."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from decimal import Decimal,localcontext
import json,time
import numpy as np
import q9_developable_bending_approximation as method
import audit_q9_curved_strain_screen as independent
import check_q9_developable_bending_reference as ref
from audit_q9_curved_best_approximation_decimal import dec,ldl,solve
core=method.core;ROOT=core.ROOT;F=core.F
def decimal_optimum(R,nr,M,B):
    n=len(R)-nr-1;D=[[dec(v) for v in row] for row in R[nr:,nr:-1]];e=[dec(v) for v in R[nr:,-1]];H=[[dec(M[i//6,j//6]) if i%6==j%6 else Decimal(0) for j in range(n)] for i in range(n)];K=[[sum((row[i]*row[j] for row in D),Decimal(0)) for j in range(n)] for i in range(n)];b=[sum((row[i]*v for row,v in zip(D,e)),Decimal(0)) for i in range(n)]
    def point(tau):
        p=solve(ldl([[v+tau*H[i][j] for j,v in enumerate(row)] for i,row in enumerate(K)]),b);norm=sum((p[i]*sum((H[i][j]*p[j] for j in range(n)),Decimal(0)) for i in range(n)),Decimal(0));return p,norm
    p,norm=point(Decimal(0));lo=Decimal(0);hi=max(K[i][i]/H[i][i] for i in range(n));iterations=0
    if norm>B*B:
        while point(hi)[1]>B*B:hi*=4
        for iterations in range(160):
            mid=(lo+hi)/2;q,nq=point(mid)
            if nq>B*B:lo=mid
            else:hi=mid
        p,norm=point(hi)
    residual=[sum((v*w for v,w in zip(row,p)),Decimal(0))-t for row,t in zip(D,e)];objective=sum((v*v for v in residual),Decimal(0))/2;reference=sum((dec(v)**2 for v in R[:,-1]),Decimal(0));return dict(error=(2*objective/reference).sqrt(),objective=objective,h1=norm.sqrt(),tau=hi if iterations else Decimal(0),iterations=iterations+1 if iterations else 0)
def main():
    target=F/'q9-developable-bending-approximation-audit.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-developable-bending-approximation-result.json';r=json.loads(parent.read_text());assert len(r['cases'])==38 and r['stop_reason'] is None;start=time.monotonic()
    assert r['source_sha256']==core.digest(ROOT/'rd/structure/main-shell/q9-curved/experiments/check_q9_developable_bending_approximation.py')
    for name,sha in r['dependencies'].items():assert core.digest(rd_paths.EXPERIMENTS/name)==sha
    for c in r['cases']:assert core.digest(ROOT/c['artifact'])==c['artifact_sha256']
    duplicates=[]
    for c in r['cases']:
        if c['n']!=1 or c['amplitude']!=0:continue
        other=next(v for v in r['cases'] if v['n']==1 and v['kappa']==c['kappa'] and v['h']==c['h'] and v['amplitude']==.2)
        with np.load(ROOT/c['artifact']) as a,np.load(ROOT/other['artifact']) as b:identical=all(np.array_equal(a[key],b[key]) for key in ['coords','xy','Rraw','Rcg','M'])
        duplicates.append(dict(kappa=c['kappa'],h=c['h'],identical=bool(identical)))
    records=[];MB=independent.E.copy();MB[3:]*=np.sqrt(2.)
    for c in r['cases']:
        if c['kappa']!=.4 or c['amplitude']!=.2 or c['h']!=.0001 or c['n'] not in [2,4,8]:continue
        a=np.load(ROOT/c['artifact']);xy,cells,coords,_,free=method.mesh(c['kappa'],c['n'],c['amplitude']);points,w,_,_=method.table(11);energy=0.;values=[0.]*len(c['solutions']);states={}
        for cell,(_,_,ids,_) in enumerate(cells):
            q=coords[ids].reshape(81);phi,J,grad=independent.data(q,c['h'],points);x=phi[:,0::3]@xy[ids,0];star=np.array([ref.tables(c['kappa'],xx,c['h']*pt[2]/2)[1] for xx,pt in zip(x,points)]);weight=w*np.linalg.det(J);star_norm=core.shapes.lam*np.trace(star,axis1=1,axis2=2)**2+2*core.shapes.mu*np.einsum('pij,pij->p',star,star);energy+=float(weight@star_norm)
            if any(v['route']=='cg' for v in c['solutions']):states[cell]=core.method.prepare(q,c['h'],0.,9);_,_,BR=core.method.evaluate(states[cell],points)
            for i,v in enumerate(c['solutions']):
                u=np.zeros((len(xy),9));u[free,:3]=a['r'][i].reshape(-1,3);u[free,3:6]=a['p'][i].reshape(-1,6)[:,:3];u[free,6:]=a['p'][i].reshape(-1,6)[:,3:]/core.LENGTH
                if v['route']=='raw':G=np.einsum('ni,pnj->pij',u[ids].reshape(27,3),grad);E=(G+G.transpose(0,2,1))/2
                else:E=np.einsum('aij,pa->pij',MB,BR[:,:,:81]@u[ids].reshape(81))
                difference=E-star;norm=core.shapes.lam*np.trace(difference,axis1=1,axis2=2)**2+2*core.shapes.mu*np.einsum('pij,pij->p',difference,difference);values[i]+=float(weight@norm)
        expected=next(v for v in json.loads((F/'q9-developable-bending-reference-result.json').read_text())['cases'] if v['kappa']==c['kappa'] and v['h']==c['h']);reference_error=abs(energy/(2*float(expected['decimal_energy']))-1)
        for v,value in zip(c['solutions'],values):
            error=float(np.sqrt(value/energy));diff=abs(error-v['observation11']['error']);records.append(dict(n=c['n'],route=v['route'],factor=v['factor'],energy_error=error,difference=diff,reference_error=reference_error,checks_pass=bool(diff<=1e-7 and reference_error<=1e-10)))
        del a,states
    chosen=next(c for c in r['cases'] if c['kappa']==.4 and c['amplitude']==.2 and c['n']==1 and c['h']==.0001);source=np.load(ROOT/chosen['artifact']);index=next(i for i,v in enumerate(chosen['solutions']) if v['route']=='cg' and v['factor']==2);v=chosen['solutions'][index];R=source['Rcg'];M=source['M'];nr=len(source['free'])*3;assert np.all(R[nr:,:nr]==0) and np.all(np.diag(R[:nr,:nr])!=0)
    with localcontext() as ctx:
        ctx.prec=80;opt=decimal_optimum(R,nr,M,dec(v['reference_h1']*2));low=dec(v['bounds']['qr_lower_relative']);high=dec(v['bounds']['qr_upper_relative']);inside=low<=opt['error']<=high;hm=dec(v['bounds']['h1_gram_lower']);ldl([[dec(value)-(hm if i==j else 0) for j,value in enumerate(row)] for i,row in enumerate(M)]);W=[[dec(value) for value in row] for row in R[:nr,:nr]];sigma=dec(v['bounds']['position_singular_lower']);ldl([[sum((W[k][i]*W[k][j] for k in range(nr)),Decimal(0))-(sigma*sigma if i==j else 0) for j in range(nr)] for i in range(nr)])
        decimal=dict(error=str(opt['error']),objective=str(opt['objective']),h1=str(opt['h1']),tau=str(opt['tau']),iterations=opt['iterations'],inside_bounds=bool(inside),positive_gram_bounds=True,qr_lower=v['bounds']['qr_lower_relative'],qr_upper=v['bounds']['qr_upper_relative'])
    result=dict(source_sha256=core.digest(__file__),parent_sha256=core.digest(parent),protocol_sha256=core.digest(F/'q9-developable-bending-approximation-audit-protocol.md'),dependencies={Path(m.__file__).name:core.digest(m.__file__) for m in [method,independent,ref,core,core.method]},duplicate_cases=duplicates,energy_cases=records,decimal=decimal,seconds=time.monotonic()-start,peak_bytes=core.memory.memory_bytes(),all_criteria_pass=bool(len(duplicates)==6 and all(v['identical'] for v in duplicates) and len(records)==8 and all(v['checks_pass'] for v in records) and decimal['inside_bounds'] and time.monotonic()-start<=300 and core.memory.memory_bytes()<=1024**3),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
