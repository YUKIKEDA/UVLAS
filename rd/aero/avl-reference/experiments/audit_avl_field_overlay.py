"""Read-only verification of saved field-overlay evidence and frozen inputs."""
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
    path=FOLDER/'avl-field-overlay-result.json';r=json.loads(path.read_text())
    assert r['source_sha256']==digest(rd_paths.sibling('check_avl_field_overlay.py'))
    assert r['protocol_sha256']==digest(FOLDER/'avl-field-overlay-experiment.md')
    assert r['geometry_addendum_sha256']==digest(FOLDER/'avl-field-overlay-exact-geometry-addendum.md')
    for name,value in r['dependency_sha256'].items():assert digest(rd_paths.sibling(name))==value
    for name,value in r['parent_sha256'].items():assert digest(FOLDER/name)==value
    assert r['methods_criteria_pass'] and len(r['checks'])==6 and all(r['checks'].values()) and r['stop_reason'] is None
    assert len(r['basis'])==2 and all(x['eligible'] and x['inputs_unchanged'] for x in r['basis'])
    assert len(r['manufactured']['checks'])==12 and all(r['manufactured']['checks'].values())
    assert len(r['manufactured']['failures'])==8 and all(x['refused'] and x['input_unchanged'] and x['output_unchanged'] for x in r['manufactured']['failures'])
    assert r['manufactured']['micro']['linf']==1 and abs(r['manufactured']['micro']['l1']-2**-40)<1e-15
    assert not r['reference_convergence_accepted'] and len(r['meshes'])==40 and len(r['series'])==12
    overlays=[];inputs=0
    for mesh in r['meshes']:
        assert mesh['load_ledger_difference']<=1e-10 and mesh['geometry_z_error']<=1e-11
        for info in mesh['inputs'].values():assert digest(ROOT/info['path'])==info['sha256'];inputs+=1
        obs=mesh['observations'];overlays.append(obs)
        assert len(obs['cell_areas'])==320 and len(obs['cell_means'])==320 and min(obs['cell_areas'])>0
        assert abs(math.fsum(obs['cell_areas'])-8.4)<=1e-10
        assert abs(math.fsum(obs['cell_integrals'][0])-obs['integrals'][0])<=1e-10
        assert np.max(abs(np.array(obs['cell_means'])*obs['cell_areas']-obs['cell_integrals'][0]))<=1e-10
        assert abs(mesh['moment_line']+mesh['moment_couple']-mesh['moment_point'])<=1e-10
    for s in r['series']:
        assert len(s['pairs'])==3 and len(s['triples'])==2
        overlays.extend(x['overlay'] for x in s['pairs']+s['triples'])
        for k,d in s['diagnostics'].items():
            assert d['differences']==[p[k]['l2'] for p in s['pairs']]
            for c in d['difference_cosines']:assert c is None or abs(c)<=1+1e-12
            if d['status']=='asymptotic_candidate':
                assert all(c>=.9 for c in d['difference_cosines']) and .5<=d['orders'][-1]<=6 and abs(d['orders'][-1]-d['orders'][-2])<=.5
    assert len(overlays)==100 and all(abs(o['area']-8.4)<=1e-10 and o['integral_preservation']<=1e-10 for o in overlays)
    assert all(o['minimum_region_area']>0 for o in overlays)
    assert r['seconds']<600 and r['peak_process_working_bytes']<512*1024**2 and path.stat().st_size<64*1024**2
    initial=json.loads((FOLDER/'avl-field-overlay-initial-result.json').read_text())
    assert initial['dependency_sha256']['avl_field_overlay.py']==digest(ROOT/'build/analysis/avl_field_overlay_initial.py')
    decimal=json.loads((FOLDER/'avl-field-overlay-decimal-result.json').read_text())
    assert decimal['source_sha256']==digest(ROOT/'build/analysis/check_avl_field_overlay_decimal.py')
    assert decimal['dependency_sha256']['avl_field_overlay.py']==digest(ROOT/'build/analysis/avl_field_overlay_decimal.py')
    changes=[]
    for before,after in zip(decimal['series'],r['series']):
        assert all(before[k]==after[k] for k in ['family','case','direction'])
        for j,(b,a) in enumerate(zip(before['pairs'],after['pairs'])):
            changes.append(dict(family=after['family'],case=after['case'],direction=after['direction'],pair=j,
                l1_difference=abs(b['field']['l1']-a['field']['l1']),l2_difference=abs(b['field']['l2']-a['field']['l2']),
                linf_difference=abs(b['field']['linf']-a['field']['linf'])))
    finite(r)
    report=dict(source_sha256=digest(Path(__file__)),result_sha256=digest(path),parent_input_files=inputs,
        meshes=40,series=12,overlays=100,manufactured_checks=12,rejected_inputs=8,method_checks=6,
        initial_result_preserved=True,decimal_to_exact_comparison=changes,all_criteria_pass=True,
        diagnostic_candidates={k:sum(s['diagnostics'][k]['status']=='asymptotic_candidate' for s in r['series']) for k in ['field','trace','cell']})
    (FOLDER/'avl-field-overlay-artifact-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
