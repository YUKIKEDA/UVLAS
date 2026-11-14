"""Archive/provenance audit and independent Decimal reconstruction of effects."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from decimal import Decimal,localcontext
import bisect,hashlib,json,zipfile
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;GIB=1024**3
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(n):return json.loads((F/n).read_text())
def main():
    raw=read('avl-allocation-factorial-runs.json');obs=read('avl-allocation-factorial-result.json');prep=read('avl-allocation-factorial-preparation.json')
    assert raw['completed'] and raw['stop_reason'] is None and len(raw['runs'])==4 and raw['seconds']<=1600
    assert raw['preparation_sha256']==digest(F/'avl-allocation-factorial-preparation.json')
    assert raw['executable_sha256']==prep['executable_sha256']==digest(ROOT/'build/reference-avl/source-symmetry16000-ks1000-spacing-pi64/bin/avl')
    assert raw['capacity_preparation_sha256']==digest(F/'avl-span1280-capacity-ks-preparation.json')
    for n,h in raw['gate_sha256'].items():assert digest(F/n)==h and read(n)['all_criteria_pass']
    for data,script in [(raw,'run_avl_allocation_factorial.py'),(obs,'check_avl_allocation_factorial.py')]:
        assert data['source_sha256']==digest(rd_paths.EXPERIMENTS/script)
        assert data['protocol_sha256']==digest(ROOT/data['protocol'])
        for n,h in data['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/n)==h
    watch=read('avl-factorial-watchdog-result.json')
    assert watch['watchdog_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/avl_factorial_watchdog.py')
    assert watch['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/check_avl_factorial_watchdog.py')
    assert watch['real_termination'] and len(watch['cases'])==10 and watch['start_cases']==3
    for n,h in obs['parent_sha256'].items():assert digest(F/n)==h
    members=0;total=0
    for r in raw['runs']:
        assert r['exit_code']==0 and r['stopped'] is None and r['seconds']<=400
        assert max(r['peak_resident_bytes'],r['maximum_group_live_bytes'])<=6*GIB
        assert r['memory_before']['cgroup_known'] and r['memory_before']['effective_available']>=10*GIB
        assert r['minimum_effective_available']>=2*GIB and r['memory_after']['cgroup_known']
        assert r['maximum_swap_used']-min(r['memory_before']['swap_used'],raw['stage_memory_before']['swap_used'])<=256*1024**2
        folder=ROOT/r['directory'];assert digest(folder/'wing.avl')==r['sha256']==digest(ROOT/r['path'])
        assert len(r['before_after_equal'])==8 and all(r['before_after_equal'].values())
        with zipfile.ZipFile(folder/'after.zip') as archive:
            assert set(archive.namelist())==set(r['archive_members']) and len(archive.namelist())==8
            for n,h in r['archive_members'].items():
                content=archive.read(n);assert hashlib.sha256(content).hexdigest()==h['sha256']==digest(folder/'out'/n) and len(content)==h['bytes'];members+=1
        for n,h in r['files'].items():assert digest(folder/n)==h['sha256'] and (folder/n).stat().st_size==h['bytes'];total+=h['bytes']
    assert total==raw['total_bytes'] and total<=128*1024**2
    assert obs['methods_criteria_pass'] and all(obs['checks'].values()) and obs['stop_reason'] is None
    assert len(obs['meshes'])==8 and len(obs['responses'])==2 and obs['seconds']<180 and obs['peak_bytes']<512*1024**2
    assert obs['richardson_order'] is None and obs['continuum_tail_estimate'] is None and obs['m24_error_budget'] is None and not obs['reference_acceptance']
    records={(r['case'],r['allocation']):r for r in obs['meshes']};maximum=0.;norm_error=0.
    with localcontext() as ctx:
        ctx.prec=60
        D=lambda v:Decimal(str(v))
        for key,r in records.items():
            assert r['eligible'] and r['residual_relative']<=1e-11
            folder=ROOT/r['raw_directory'];lines=(folder.parent/'wing.avl').read_text().splitlines();boundaries=[float(lines[j+1].split()[1]) for j,line in enumerate(lines) if line=='SECTION'];counts=[0]*10
            for line in (folder/'uv-strips.txt').read_text().splitlines()[2:]:
                v=line.split();counts[bisect.bisect_right(boundaries[1:],(float(v[1])+float(v[2]))/2)]+=1
            assert counts==r['actual_interval_counts']
            for n,h in r['files'].items():assert digest(ROOT/r['derived_directory']/n)==h['sha256']
        def vector(r,c):
            bands=r['bands'];root=[D(b['width']).sqrt() for b in bands]
            if c=='tau_p1':return vector(r,'tau0')+[(D(b['width'])/3).sqrt()*D(b['tau1']) for b in bands]
            field={'tau0':'tau0','k0':'kbar','j_line':'jbar','j_point':'jcbar'}[c]
            return [w*D(b[field]) for w,b in zip(root,bands)]
        for response in obs['responses']:
            case=response['case']
            for channel,data in response['channels'].items():
                A,B,C,E=[vector(records[(case,a)],channel) for a in ['uniform','inner-only','tip-only','tip-rich']]
                effects={'inner':[b-a for a,b in zip(A,B)],'tip':[c-a for a,c in zip(A,C)],'interaction':[e-b-c+a for a,b,c,e in zip(A,B,C,E)],'combined':[e-a for a,e in zip(A,E)]}
                for name,v in effects.items():
                    saved=data['effects'][name];err=max(abs(float(x)-y) for x,y in zip(v,saved['difference']));maximum=max(maximum,err)
                    assert err<=1e-12
                    norm=sum((x*x for x in v),Decimal(0)).sqrt();err=abs(float(norm)-saved['norm']);norm_error=max(norm_error,err);assert err<=1e-12
    result=dict(source_sha256=digest(__file__),runs_sha256=digest(F/'avl-allocation-factorial-runs.json'),result_sha256=digest(F/'avl-allocation-factorial-result.json'),
        archive_members=members,raw_bytes=total,decimal_effect_maximum_difference=maximum,decimal_norm_maximum_difference=norm_error,all_criteria_pass=True,reference_acceptance=False)
    (F/'avl-allocation-factorial-artifact-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
