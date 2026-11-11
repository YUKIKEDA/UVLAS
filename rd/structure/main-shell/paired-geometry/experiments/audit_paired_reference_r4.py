"""Independent acceptance arithmetic and immutable evidence graph for R4."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import ast,hashlib,json,platform,time
import numpy as np
import scipy,flint
ROOT=rd_paths.ROOT
NAMES=['mass-paired-quality','mass-paired-coefficient-input','mass-paired-reference','paired-lift-exact',
    'mass-paired-body','mass-paired-interface','mass-paired-mechanics','mass-paired-orientation','mass-paired-capacity']

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def require(condition,message):
    if not condition:raise AssertionError(message)
def maximum(values):return max(float(x) for x in values)

def quality_checks(records,expected,extra=False):
    require(len(records)==expected,'case count')
    for r in records:
        o=r['observation'];require(maximum(o['stationarity'])<=1e-8,'original stationarity')
        require(r['rcond']>=1e-10 and r['mass_min']>1e-10,'conditioning')
        ns=r.get('ns',r['mass_root_ranks'][0]);require(r['stiffness_ranks']==[3*ns]*3 and r['mass_root_ranks']==[ns]*3,'completed roots rank')
        require(max(max(a.values()) for a in r['inclusion'])<=1e-8,'reference coarse inclusion')
        require(all(r['checks'].values()),'reported numerical gates')
        if not extra:
            require(r['trace']<=1e-9 and max(r['equivalence'].values())<=1e-12,'same operator and trace')
            require(r['cold_case_seconds']<=300 and r['cold_preparation']>=r['preparation'],'cold resource')
            require(r['normalization']['verified_congruence_bound']<=1e-20 and r['normalization']['max_radius']<=1e-25,'coefficient arithmetic')
    groups={tuple(r.get('parameters',[])) for r in records}
    for pars in groups:
        subset=[r for r in records if tuple(r.get('parameters',[]))==pars]
        for n in {r['n'] for r in subset}:
            pair=[next(r for r in subset if r['n']==n and r['h']==h) for h in [.01,.0001]]
            for quantity in ['error','position']:
                ratio=np.maximum(pair[1]['observation'][quantity],1e-12)/np.maximum(pair[0]['observation'][quantity],1e-12)
                require(np.all(ratio<=2),'thin amplification '+str((pars,n,quantity,ratio)))
        finest=max(r['n'] for r in subset);coarsest=min(r['n'] for r in subset)
        for h in [.01,.0001]:
            a=next(r for r in subset if r['n']==coarsest and r['h']==h);b=next(r for r in subset if r['n']==finest and r['h']==h)
            for quantity in ['error','position']:
                require(maximum(b['observation'][quantity])<=.01,'finest accuracy')
                require(np.all(np.array(b['observation'][quantity])<=np.maximum(a['observation'][quantity],1e-12)),'spatial decrease')

def artifact_pairs(value):
    if isinstance(value,dict):
        if isinstance(value.get('artifact'),str) and 'artifact_sha256' in value:yield value['artifact'],value['artifact_sha256']
        for v in value.values():yield from artifact_pairs(v)
    elif isinstance(value,list):
        for v in value:yield from artifact_pairs(v)

def main():
    start=time.monotonic();dest=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-reference-r4-audit.json';require(not dest.exists(),'immutable destination')
    reports={name:json.loads((rd_paths.RECORDS/f'{name}-result.json').read_text(encoding='utf-8')) for name in NAMES}
    for name,r in reports.items():require(r['passed'] and not r.get('stops') and not r.get('failures'),'report '+name)
    q=reports['mass-paired-quality'];quality_checks(q['records'],24)
    require(sum(r['family']=='primary' for r in q['records'])==16 and sum(r['family']=='nonquadratic' for r in q['records'])==8,'quality coverage')
    quality_checks(reports['mass-paired-coefficient-input']['records'],4,True)
    reference=reports['mass-paired-reference'];require(len(reference['records'])==21,'reference coverage')
    for r in reference['records']:
        if 'metrics' in r:require(max(r['metrics'].values())<=1e-8,'reference metrics')
        if 'certificate' in r:require(r['certificate']['global_det_J_lower']>0 and r['certificate']['max_half_thickness_curvature']<1,'reference tube')
    lift=reports['paired-lift-exact'];require(sum(r['fields'] for r in lift['records'])==256 and all(r['coefficient_residual_exactly_zero'] and r['zero_director_negative_controls']>0 for r in lift['records']),'exact identities')
    body=reports['mass-paired-body'];require(len(body['records'])==4,'physical body coverage')
    for r in body['records']:
        require(maximum(r['observation']['error'])<=.01 and maximum(r['observation']['position'])<=.01,'physical accuracy')
        require(maximum(r['observation']['stationarity'])<=1e-8 and r['physical_vs_projection']<=1e-8,'physical stationarity')
        require(r['external_work_error']<=.01 and np.linalg.norm(r['loads']['resultant'])<=1e-10,'physical work and force')
    interface=reports['mass-paired-interface'];require([r['elements'] for r in interface['records']]==[[1,1],[1,2]],'physical interface coverage')
    for r in interface['records']:
        require(max(r['equilibrium'],r['interface_work'],r['solution_error'])<=1e-8 and r['rcond']>=1e-10,'physical interface balance')
        require(len(set(r['trace']['ranks']))==1 and max(r['trace']['independent_interface_holdout'],r['trace']['fixed_holdout'])<=1e-9,'all trace')
        require(r['stiffness_ranks']==[3*r['trace']['scalar_dimension']]*3,'interface stiffness rank')
    mechanical=reports['mass-paired-mechanics'];require(len(mechanical['records'])==6,'mechanics coverage')
    for r in mechanical['records']:
        if 'metrics' in r:require(max(r['metrics'].values())<=1e-8,'mechanics metrics')
        if r['case']=='mass_spectrum_inertia':require(r['free_nullities']==[6]*3 and max(r[k] for k in ['eigen_backward','mass_orthogonality','spatial_inertia'])<=1e-8 and r['eigen_min']>0,'mass modes')
        if r['case']=='negative_prestress_retained':require(r['negative_count']>0 and r['minimum']<0 and not r['filtered'],'signed spectrum')
        if r['case']=='same_system_restart_reference':require(r['all_modes_excited'] and r['bitwise_replay'] and r['negative_controls']==16 and r['energy_error']<=1e-8 and not r['product_integrator_adopted'],'restart')
    orientation=reports['mass-paired-orientation'];require(len(orientation['records'])==3,'orientation coverage')
    require(max(orientation['records'][0]['metrics'].values())<=1e-8 and orientation['records'][1]['mass_error']<=1e-10 and orientation['records'][1]['forgotten_layer_flip_defect']>1e-5,'quadrature layers')
    require(max(orientation['records'][2]['errors'])<=1e-10 and min(orientation['records'][2]['wrong_signed_measure_defects'])>1,'all face work')
    cap=reports['mass-paired-capacity']['record'];require(cap['total_dofs']>=100000 and cap['memory_after']['peak_rss']<4*1024**3 and max(cap['energy_error'],cap['backward_error'])<=1e-8 and all(cap['checks'].values()),'actual capacity')
    # Every declared hash is checked; then freeze transitive local Python imports.
    manifest={};artifacts={}
    for name,r in reports.items():
        path=rd_paths.RECORDS/f'{name}-result.json';manifest[path.relative_to(ROOT).as_posix()]=digest(path)
        for path,sha in r['sources'].items():
            require(digest(ROOT/path)==sha,'changed source '+path);manifest[path]=sha
        for path,sha in artifact_pairs(r):require(digest(ROOT/path)==sha,'changed artifact '+path);artifacts[path]=sha
    directory=rd_paths.EXPERIMENTS;queue=[ROOT/p for p in manifest if p.endswith('.py')];seen=set()
    while queue:
        path=queue.pop()
        if path in seen:continue
        seen.add(path);manifest[path.relative_to(ROOT).as_posix()]=digest(path)
        tree=ast.parse(path.read_text(encoding='utf-8-sig'))
        modules=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):modules += [a.name.split('.')[0] for a in node.names]
            elif isinstance(node,ast.ImportFrom) and node.module:modules.append(node.module.split('.')[0])
        for module in modules:
            other=directory/(module+'.py')
            if other.exists():queue.append(other)
    for path in ['rd/structure/main-shell/paired-geometry/paired-reference-coefficient-definition.md','rd/structure/main-shell/paired-geometry/paired-reference-local-global-contract.md','rd/structure/main-shell/paired-geometry/records/paired-reference-r4-protocol.md']:
        manifest[path]=digest(ROOT/path)
    manifest[Path(__file__).relative_to(ROOT).as_posix()]=digest(Path(__file__))
    finest=[r for r in q['records'] if r['n']==(4 if r['family']=='primary' else 2)]
    metrics=dict(quality_configurations=28,quality_rhs=84,primary_configurations=16,nonquadratic_configurations=12,
        finest_energy_error=max(max(r['observation']['error']) for r in finest+reports['mass-paired-coefficient-input']['records'][2:]),
        finest_position_error=max(max(r['observation']['position']) for r in finest+reports['mass-paired-coefficient-input']['records'][2:]),
        max_thin_amplification=max(max(c['values']) for c in q['comparisons']+reports['mass-paired-coefficient-input']['comparisons'] if c['kind']=='thin_amplification'),
        max_original_stationarity=max(max(r['observation']['stationarity']) for r in q['records']),
        max_cold_quality_case_seconds=max(r['cold_case_seconds'] for r in q['records']),actual_capacity_dofs=cap['total_dofs'],capacity_peak_rss=cap['memory_after']['peak_rss'])
    result=dict(passed=True,R4_evidence_gates_passed=True,main_shell_adopted=False,product_implementation_start=False,
        metrics=metrics,environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,python_flint=flint.__version__),
        sources=manifest,artifacts=artifacts,counts=dict(reports=len(reports),sources=len(manifest),artifacts=len(artifacts)),seconds=time.monotonic()-start)
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(metrics=metrics,counts=result['counts']),indent=2))

if __name__=='__main__':main()
