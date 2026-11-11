"""Independent MR-R1 observations and provenance for the T-S36 global trial."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,math,time
import numpy as np
import audit_q9_shear_projection_plate as observer
import check_q9_plate_convergence as memory
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=F/'q9-trace-shear-plate-audit.json'
    if output.exists():raise FileExistsError(output)
    start=time.monotonic();parent=F/'q9-trace-shear-plate-result.json';data=json.loads(parent.read_text());refs=json.loads((F/'plate-modal-exact-reference-result.json').read_text());cases=[];observations=[]
    assert data['source_sha256']==digest(ROOT/'rd/structure/main-shell/q9-edge-trace/experiments/check_q9_trace_shear_plate.py') and data['protocol_sha256']==digest(ROOT/data['protocol'])
    assert data['candidate_sha256']==digest(ROOT/'rd/structure/main-shell/q9-edge-trace/q9-trace-shear-candidate.md') and data['stop_reason'] is None
    for n,h in data['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/n)==h
    for n,h in data['gate_sha256'].items():
        saved=json.loads((F/n).read_text());assert digest(F/n)==h and saved['all_criteria_pass']
        for dep,value in saved.get('dependency_sha256',{}).items():assert digest(rd_paths.EXPERIMENTS/dep)==value
    ancestor=F/'q9-distorted-plate-result.json';assert refs['parent_sha256']==digest(ancestor)
    for dep,value in json.loads(ancestor.read_text())['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/dep)==value
    assert len(data['cases'])==8 and {(c['h'],c['n']) for c in data['cases']}=={(h,n) for h in [.01,.001] for n in [2,4,8,16]}
    assert data['seconds']<=1200 and data['peak_bytes']<=2*1024**3
    for c in data['cases']:
        assert c['eligible'] and c['assembly_pass'] and c['stiffness_7_vs_9']<=1e-8 and c['load_7_vs_9']<=1e-8 and c['symmetry']<=1e-12
        assert c['seconds']<=240 and c['peak_bytes']<=2*1024**3
        assert all(v<=1e-8 for group in c['precision_comparisons'].values() for v in group.values()) and all(v<=1e-5 for v in c['quadrature_response'].values())
        for group in c['solves'].values():
            for r in group.values():
                last=r['history'][-1];assert r['eligible'] and len(r['history'])<=6 and last['relative_residual']<=1e-9 and last['backward_normwise']<=1e-14 and last['backward_componentwise']<=1e-14
        file=ROOT/c['artifact'];assert digest(file)==c['artifact_sha256'];arrays=np.load(file);ref=next(r for r in refs['cases'] if r['h']==c['h']);path=ref['decimal_paths']['100'];independent=dict(coefficients=[float(v) for v in path['coefficients']],center=float(path['coefficients'][2]),compliance=float(path['compliance']),pressure=ref['rounded_reference']['pressure']);assert independent==c['reference']
        n=c['n']
        for j in range(n+1):
            for i in range(n+1):
                shift=.2/n*(-1)**(i+j)*math.sin(math.pi*i/n)*math.sin(math.pi*j/n) if 0<i<n and 0<j<n else 0.
                assert np.max(abs(arrays['vertices'][j,i]-[i/n+shift,j/n]))<=1e-13
        chosen=None;maximum=0.
        for order in [7,9]:
            for digits in [70,90]:
                r=observer.evaluate(n,arrays['vertices'],arrays['nodes'],arrays[f'q{order}_{digits}'],arrays[f'f{order}'],independent);old=c['solves'][str(order)][str(digits)]['observables']
                dc=abs((r['center']-old['center'])/independent['center']);dw=abs((r['compliance']-old['compliance'])/independent['compliance']);de=abs(r['l2_relative']-old['l2_relative']);maximum=max(maximum,de)
                assert dc<=1e-10 and dw<=1e-10 and de<=1e-8*max(old['l2_relative'],1e-12)
                observations.append(dict(n=n,h=c['h'],order=order,digits=digits,center_difference=dc,work_difference=dw,l2_difference=de))
                if order==9 and digits==90:chosen=r
        cases.append(dict(n=n,h=c['h'],independent=chosen,maximum_l2_difference=maximum))
        assert time.monotonic()-start<=120 and memory.memory_bytes()<=512*1024**2
    checks=[]
    for h in [.01,.001]:
        a=next(r for r in cases if r['n']==8 and r['h']==h)['independent'];b=next(r for r in cases if r['n']==16 and r['h']==h)['independent'];rate=math.log2(a['l2_relative']/b['l2_relative']);saved=next(r for r in data['thickness_checks'] if r['h']==h);accuracy=all(b[k]<=.01 for k in ['center_relative','l2_relative','compliance_relative']);order=b['l2_relative']<=1e-8 or rate>=1.5
        assert accuracy==saved['accuracy_pass'] and order==saved['order_pass'] and abs(rate-saved['last_l2_order'])<=1e-8
        checks.append(dict(h=h,accuracy_pass=accuracy,order_pass=order,last_l2_order=rate))
    mid=next(r for r in cases if r['n']==8 and r['h']==.01)['independent'];thin=next(r for r in cases if r['n']==8 and r['h']==.001)['independent'];uniform=thin['l2_relative']<=2*mid['l2_relative']+1e-8;assert uniform==data['thin_uniformity_pass']
    acceptance=uniform and all(c['accuracy_pass'] and c['order_pass'] for c in checks);assert acceptance==data['all_criteria_pass']
    assert sum((ROOT/c['artifact']).stat().st_size for c in data['cases'])<=64*1024**2
    result=dict(source_sha256=digest(__file__),parent_sha256=digest(parent),reference_sha256=digest(F/'plate-modal-exact-reference-result.json'),protocol_sha256=digest(F/'q9-trace-shear-plate-experiment.md'),
        dependency_sha256={Path(m.__file__).name:digest(m.__file__) for m in [observer,memory]},cases=cases,observations=observations,thickness_checks=checks,thin_uniformity_pass=uniform,thin_to_mid_l2_ratio=thin['l2_relative']/mid['l2_relative'],
        seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes(),all_criteria_pass=True,experiment_acceptance_pass=acceptance,product_shell_adopted=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_criteria_pass=True,experiment_acceptance_pass=acceptance,thin_to_mid_l2_ratio=result['thin_to_mid_l2_ratio'],seconds=result['seconds'])))
if __name__=='__main__':main()
