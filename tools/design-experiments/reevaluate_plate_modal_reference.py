"""Re-evaluate frozen fields with independently rounded MR-R1 references."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,math,time
import numpy as np
import check_q9_distorted_plate as primary
import audit_q9_shear_projection_plate as independent
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=F/'plate-modal-reference-reevaluation.json'
    if output.exists():raise FileExistsError(output)
    start=time.monotonic();reference_path=F/'plate-modal-exact-reference-result.json';raw=json.loads(reference_path.read_text());assert raw['all_criteria_pass']
    assert raw['source_sha256']==digest(ROOT/'tools/design-experiments/check_plate_modal_exact_reference.py')
    names={'full':'q9-distorted-plate-result.json','P-S36':'q9-shear-projection-plate-result.json'};hashes={n:digest(F/n) for n in names.values()};old_audit=F/'q9-shear-projection-plate-artifact-audit.json';old_audit_hash=digest(old_audit);datasets=[]
    for label,name in names.items():
        parent=json.loads((F/name).read_text());records=[]
        for c in parent['cases']:
            assert c['eligible'] and digest(ROOT/c['artifact'])==c['artifact_sha256'];arrays=np.load(ROOT/c['artifact']);refraw=next(r for r in raw['cases'] if r['h']==c['h'])
            r=refraw['rounded_reference'];r2=dict(coefficients=[float(v) for v in refraw['decimal_paths']['100']['coefficients']],center=float(refraw['decimal_paths']['100']['coefficients'][2]),compliance=float(refraw['decimal_paths']['100']['compliance']),pressure=r['pressure'])
            assert r==r2
            cells=primary.plate.connectivity(c['n']);observations=[]
            for order in [7,9]:
                for digits in [70,90]:
                    q=arrays[f'q{order}_{digits}'];f=arrays[f'f{order}']
                    a=primary.observables(c['n'],arrays['vertices'],arrays['nodes'],cells,q,f,r)
                    b=independent.evaluate(c['n'],arrays['vertices'],arrays['nodes'],q,f,r2)
                    center=abs((a['center']-b['center'])/r['center']);work=abs((a['compliance']-b['compliance'])/r['compliance']);delta=abs(a['l2_relative']-b['l2_relative']);limit=1e-8*max(a['l2_relative'],1e-12)
                    exact=Q(refraw['coefficients_rational'][2]);rounding=abs(Q.from_float(r['center'])/exact-1);bound=float(rounding/(1-rounding))*(1+a['l2_relative'])
                    checks=dict(center=center<=1e-10,work=work<=1e-10,l2=delta<=limit,reference_rounding_budget=bound<=limit)
                    observations.append(dict(order=order,digits=digits,primary=a,independent=b,defects=dict(center=center,work=work,l2=delta),reference_rounding_l2_bound=bound,checks=checks))
            chosen=next(o for o in observations if o['order']==9 and o['digits']==90)
            records.append(dict(n=c['n'],h=c['h'],artifact=c['artifact'],artifact_sha256=c['artifact_sha256'],observations=observations,primary=chosen['primary'],independent=chosen['independent']))
            if time.monotonic()-start>120 or memory.memory_bytes()>512*1024**2:raise RuntimeError('reevaluation_resource_limit')
        checks=[]
        for h in [.01,.001]:
            a=next(r for r in records if r['n']==8 and r['h']==h);b=next(r for r in records if r['n']==16 and r['h']==h)
            rate=math.log2(a['primary']['l2_relative']/b['primary']['l2_relative']);ri=math.log2(a['independent']['l2_relative']/b['independent']['l2_relative']);accuracy=all(b['primary'][k]<=.01 for k in ['center_relative','l2_relative','compliance_relative']);order=b['primary']['l2_relative']<=1e-8 or rate>=1.5
            assert accuracy==all(b['independent'][k]<=.01 for k in ['center_relative','l2_relative','compliance_relative'])
            assert order==(b['independent']['l2_relative']<=1e-8 or ri>=1.5)
            checks.append(dict(h=h,accuracy_pass=accuracy,order_pass=order,last_l2_order=rate,independent_order=ri,order_audit_pass=abs(rate-ri)<=1e-8))
        mid=next(r for r in records if r['n']==8 and r['h']==.01);thin=next(r for r in records if r['n']==8 and r['h']==.001)
        uniform=thin['primary']['l2_relative']<=2*mid['primary']['l2_relative']+1e-8
        assert uniform==(thin['independent']['l2_relative']<=2*mid['independent']['l2_relative']+1e-8)
        datasets.append(dict(name=label,cases=records,thickness_checks=checks,thin_uniformity_pass=uniform,thin_to_mid_l2_ratio=thin['primary']['l2_relative']/mid['primary']['l2_relative'],
            independent_observation_audit_pass=all(all(o['checks'].values()) for r in records for o in r['observations']) and all(c['order_audit_pass'] for c in checks),experiment_acceptance_pass=uniform and all(c['accuracy_pass'] and c['order_pass'] for c in checks)))
    assert all(digest(F/n)==v for n,v in hashes.items()) and digest(old_audit)==old_audit_hash
    result=dict(source_sha256=digest(__file__),reference_sha256=digest(reference_path),protocol_sha256=digest(ROOT/'docs/plate-modal-reference-accuracy-contract.md'),
        parent_sha256=hashes,old_audit_sha256=old_audit_hash,dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [primary,independent,memory]},datasets=datasets,
        seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),old_files_unchanged=True,all_criteria_pass=all(d['independent_observation_audit_pass'] for d in datasets),product_shell_adopted=False)
    output.write_text(json.dumps(result,indent=2)+'\n');assert output.stat().st_size<=16*1024**2
    print(json.dumps(dict(all_criteria_pass=result['all_criteria_pass'],seconds=result['seconds'],datasets=[{k:d[k] for k in ['name','thin_to_mid_l2_ratio','independent_observation_audit_pass','experiment_acceptance_pass']} for d in datasets]),indent=2))
    if not result['all_criteria_pass']:raise SystemExit(1)
if __name__=='__main__':main()
