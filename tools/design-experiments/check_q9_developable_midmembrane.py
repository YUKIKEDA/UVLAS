"""Leading mid-surface membrane observation in the interpolated reference field."""
from pathlib import Path
import json,time
import numpy as np
import q9_developable_bending_approximation as method
import audit_q9_curved_strain_screen as independent
core=method.core;ROOT=core.ROOT;F=core.F
def main():
    target=F/'q9-developable-midmembrane-result.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-developable-bending-approximation-result.json';r=json.loads(parent.read_text());assert r['numeric_checks_pass'];start=time.monotonic();records=[];N,g,w=method.surface_table(9);p=np.polynomial.legendre.leggauss(9)[0];points=np.array([[s,t,0.] for s in p for t in p]);CI=method.material.CI;MB=independent.E.copy();MB[3:]*=np.sqrt(2.)
    for k,n in [(0.,4),(.4,4),(.4,8)]:
        h=.0001;xy,cells,coords,u,free=method.mesh(k,n,.2);values=dict(raw=0.,cg=0.);positive=True
        for _,_,ids,_ in cells:
            q=coords[ids].reshape(81);_,J,grad=independent.data(q,h,points);A=J[:,:,:2];L=np.linalg.cholesky(A.transpose(0,2,1)@A);U=A@np.linalg.inv(L).transpose(0,2,1);H=np.stack((np.einsum('pi,kij,pj->pk',U[:,:,0],MB,U[:,:,0]),np.einsum('pi,kij,pj->pk',U[:,:,1],MB,U[:,:,1]),np.sqrt(2.)*np.einsum('pi,kij,pj->pk',U[:,:,0],MB,U[:,:,1])),axis=1);Cmem=np.linalg.inv(H@CI@H.transpose(0,2,1));positive=positive and np.min(np.linalg.eigvalsh(Cmem))>0
            du=np.einsum('ni,pnj->pij',u[ids].reshape(27,3),grad);E=(du+du.transpose(0,2,1))/2;raw=np.einsum('kij,pij->pk',MB,E);state=core.method.prepare(q,h,0.,9);_,_,B=core.method.evaluate(state,points);cg=B[:,:,:81]@u[ids].reshape(81);weight=w*np.linalg.det(J)*2/h
            for name,field in [('raw',raw),('cg',cg)]:
                observed=np.einsum('pki,pi->pk',H,field);values[name]+=float(h*np.einsum('p,pi,pij,pj->',weight,observed,Cmem,observed)/2)
        case=next(c for c in r['cases'] if c['kappa']==k and c['amplitude']==.2 and c['n']==n and c['h']==h);relative={name:float(np.sqrt(value/case['reference_energy'])) for name,value in values.items()};checks=dict(metric=bool(positive),finite=bool(all(np.isfinite(v) for v in relative.values())),flat=k!=0 or max(relative.values())<=1e-8);records.append(dict(kappa=k,n=n,h=h,membrane_scale=values,normalized_membrane_scale=relative,full_interpolation_error=case['reference_interpolation'],interpolation_position_error=case['interpolation_position_error'],checks=checks));print(json.dumps(records[-1]),flush=True)
    result=dict(source_sha256=core.digest(__file__),parent_sha256=core.digest(parent),protocol_sha256=core.digest(F/'q9-developable-midmembrane-protocol.md'),dependencies={Path(m.__file__).name:core.digest(m.__file__) for m in [method,core,core.method,independent]},cases=records,seconds=time.monotonic()-start,peak_bytes=core.memory.memory_bytes(),all_criteria_pass=bool(len(records)==3 and all(all(c['checks'].values()) for c in records) and time.monotonic()-start<=300 and core.memory.memory_bytes()<=1024**3),product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
