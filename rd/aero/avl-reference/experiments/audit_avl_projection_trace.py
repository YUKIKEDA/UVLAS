"""Audit saved projection/trace evidence without rerunning numerical solves."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,math
import numpy as np

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def finite(x):
    if isinstance(x,float):assert math.isfinite(x)
    elif isinstance(x,dict):
        for v in x.values():finite(v)
    elif isinstance(x,list):
        for v in x:finite(v)

def main():
    path=FOLDER/'avl-projection-trace-result.json';r=json.loads(path.read_text())
    assert r['source_sha256']==digest(rd_paths.sibling('check_avl_projection_trace.py'))
    assert r['protocol_sha256']==digest(FOLDER/'avl-projection-trace-experiment.md')
    for name,value in r['dependency_sha256'].items():assert digest(rd_paths.sibling(name))==value
    for name,value in r['parent_sha256'].items():assert digest(FOLDER/name)==value
    for name,value in r['input_sha256'].items():assert digest(ROOT/name)==value
    assert r['methods_criteria_pass'] and len(r['checks'])==8 and all(r['checks'].values()) and r['stop_reason'] is None
    assert not r['reference_convergence_accepted'] and not r['manufactured_projection_is_physical_error_estimate']
    assert len(r['independent']['uniform'])==8 and len(r['independent']['quadrature'])==12
    assert len(r['independent']['checks'])==4 and all(r['independent']['checks'].values())
    assert all(v['maximum_difference']<=1e-11 and v['polynomials']==10 for v in r['independent']['quadrature'])
    assert len(r['projections'])==80 and len(r['projection_series'])==24 and len(r['actual_trace'])==36
    keys=set()
    for p in r['projections']:
        key=(p['family'],p['case'],p['m'],p['n'],p['kind']);assert key not in keys;keys.add(key)
        assert len(p['coefficients'])==p['n'] and all(len(row)==p['m'] for row in p['coefficients'])
        assert p['orthogonality_max']<=1e-12 and p['integral_difference']<=1e-10
        assert abs(p['relative_error_l2']-p['error_l2']/p['reference_l2'])<=1e-12
        o=p['observations'];assert len(p['cell_means'])==len(o['cell_areas'])==320
        assert abs(math.fsum(o['cell_areas'])-8.4)<=1e-10 and abs(math.fsum(o['cell_integrals'][0])-p['projected_integral'])<=1e-10
    for s in r['projection_series']:
        e=np.array(s['known_errors_l2']);assert s['known_errors_decrease']==bool(np.all(np.diff(e)<0))
        assert np.max(abs(np.log2(e[:-1]/e[1:])-s['known_error_orders']))<=1e-12
        assert len(s['pairs'])==3 and len(s['field_diagnostic']['difference_cosines'])==2
    for t in r['actual_trace']:
        assert max(t['preservation'],t['lift_change_difference'],t['previous_trace_difference'])<=1e-10
        rows=t['intervals'];sums={}
        sums['signed']=math.fsum((x['high']-x['low'])*x['difference'] for x in rows)
        sums['l1']=math.fsum((x['high']-x['low'])*abs(x['difference']) for x in rows)
        sums['l2_squared']=math.fsum((x['high']-x['low'])*x['difference']**2 for x in rows)
        for k,v in sums.items():assert abs(v-t['total'][k])<=1e-11
        for group in [t['bands'],list(t['zones'].values())]:
            assert all(abs(math.fsum(x[k] for x in group)-t['total'][k])<=1e-10 for k in ['positive','negative','signed','l1','l2_squared'])
        total=t['total'];assert abs(total['positive']+total['negative']-total['signed'])<=1e-11
        assert abs(total['signed'])<=total['l1']+1e-11
        assert total['l1']<=math.sqrt(8*total['l2_squared'])+1e-11
        assert abs(total['cancellation_ratio']-total['l1']/max(abs(total['signed']),1e-12))<=1e-12
    assert r['seconds']<600 and r['peak_process_working_bytes']<512*1024**2 and path.stat().st_size<64*1024**2
    finite(r)
    report=dict(source_sha256=digest(Path(__file__)),result_sha256=digest(path),method_checks=8,
        independent_gauss_cases=12,uniform_analytic_cases=8,projections=80,projection_series=24,actual_trace_pairs=36,
        input_files=len(r['input_sha256']),known_error_decreasing_series=sum(s['known_errors_decrease'] for s in r['projection_series']),
        field_direction_candidates=sum(s['field_diagnostic']['status']=='asymptotic_candidate' for s in r['projection_series']),
        all_criteria_pass=True)
    (FOLDER/'avl-projection-trace-artifact-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
