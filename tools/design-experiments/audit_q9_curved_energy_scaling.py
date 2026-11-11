"""Independent unscaled-covariant KKT construction of energy observation metric."""
from pathlib import Path
import json,time
import numpy as np
import q9_curved_energy_scaling as method
import audit_q9_curved_strain_screen as independent
core=method.core;ROOT=core.ROOT;F=core.F
def main():
    target=F/'q9-curved-energy-scaling-audit.json'
    if target.exists():raise FileExistsError(target)
    parent=F/'q9-curved-energy-scaling-result.json';r=json.loads(parent.read_text());assert r['bounds_checks_pass']
    assert r['source_sha256']==core.digest(ROOT/'tools/design-experiments/check_q9_curved_energy_scaling.py')
    assert r['contract_sha256']==core.digest(ROOT/'docs/q9-curved-energy-scaling-contract.md');assert r['protocol_sha256']==core.digest(F/'q9-curved-energy-scaling-experiment.md')
    for name,sha in r['dependencies'].items():assert core.digest(ROOT/'tools/design-experiments'/name)==sha
    start=time.monotonic();records=[];fixed={}
    for c in r['cases']:
        assert core.digest(ROOT/c['artifact'])==c['artifact_sha256'];xy,cells,coords,_,_,_=core.geometry.make_mesh(c['geometry'],c['n']);points=[[-.71,.22,-.31],[.18,-.64,.73],[.62,.43,-.27]];maximum=0.
        for cell in [0,len(cells)-1]:
            q=coords[cells[cell][2]].reshape(81);_,J,_=independent.data(q,c['h'],points);A=J[:,:,:2];d=J[:,:,2]*2/c['h'];H=2*np.einsum('pia,kij,pj->pak',A,core.method.MB,d)
            Jp,_=core.method.base.raw.kinematics(q,c['h'],0.,points);Ap=Jp[:,:,:2];dp=Jp[:,:,2]*2/c['h'];Gp=Ap.transpose(0,2,1)@Ap;Hp=2*np.einsum('pia,kij,pj->pak',Ap,core.method.MB,dp);Hp=np.linalg.solve(np.linalg.cholesky(Gp),Hp)/np.linalg.norm(dp,axis=1)[:,None,None];weighted,_,_,_=method.energy_rows(Hp)
            for h,w in zip(H,weighted):
                saddle=np.block([[method.CM,h.T],[h,np.zeros((2,2))]]);P=np.linalg.solve(saddle,np.vstack((np.zeros((6,6)),h)))[:6];energy=P.T@method.CM@P/method.MU;primary=w.T@w;maximum=max(maximum,float(np.linalg.norm(energy-primary)/np.linalg.norm(energy)))
        budget_error=0.;energy_excess=0.
        for x in c['records']:
            if x['family']=='fixed':
                if c['geometry'] not in fixed:fixed[c['geometry']]=x['requested_epsilon']
                budget_error=max(budget_error,abs(x['requested_epsilon']/fixed[c['geometry']]-1));budget=.5*method.MU*c['volume9']*x['requested_epsilon']**2
            else:
                scale=c['membrane_scale'] if x['family']=='membrane' else c['bending_scale'];budget=x['eta']*scale;budget_error=max(budget_error,abs(.5*method.MU*c['volume9']*x['requested_epsilon']**2/budget-1))
            energy_excess=max(energy_excess,max(0.,x['certified_component_energy_upper9']/budget-1,x['certified_component_energy_upper11']/budget-1))
        checks=dict(metric=maximum<=1e-11,budget=budget_error<=1e-12,upper_energy=energy_excess<=1e-12);records.append(dict(geometry=c['geometry'],n=c['n'],h=c['h'],metric_difference=maximum,budget_difference=budget_error,upper_energy_excess=energy_excess,checks={k:bool(v) for k,v in checks.items()}))
        if time.monotonic()-start>300 or core.memory.memory_bytes()>1024**3:raise RuntimeError('resource_limit')
    result=dict(source_sha256=core.digest(__file__),parent_sha256=core.digest(parent),protocol_sha256=core.digest(F/'q9-curved-energy-scaling-audit-protocol.md'),dependencies={Path(m.__file__).name:core.digest(m.__file__) for m in [method,core,core.method,independent]},cases=records,seconds=time.monotonic()-start,peak_bytes=core.memory.memory_bytes(),all_criteria_pass=bool(len(records)==16 and all(all(c['checks'].values()) for c in records)),kkt_gate_replaced=False,product_shell_adopted=False)
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(all_criteria_pass=result['all_criteria_pass'],cases=len(records),metric=max(c['metric_difference'] for c in records),budget=max(c['budget_difference'] for c in records),energy_excess=max(c['upper_energy_excess'] for c in records),seconds=result['seconds'])))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
