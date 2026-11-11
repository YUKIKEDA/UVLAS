"""Independent archive/provenance/resource/allocation and sensitivity audit."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from decimal import Decimal,localcontext
import bisect,hashlib,json,math,zipfile
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;GIB=1024**3
NAMES={'total.txt','strips.txt','elements.txt','uv-rows.txt','uv-forces.txt','uv-strips.txt','uh-meta.txt','uh-guard-tests.txt'}
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(n):return json.loads((F/n).read_text())
def main():
    raw=read('avl-tip-redistribution-runs.json');obs=read('avl-tip-redistribution-result.json');prep=read('avl-tip-redistribution-preparation.json');watch=read('avl-tip-watchdog-result.json')
    assert raw['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/run_avl_tip_redistribution.py')
    assert obs['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/check_avl_tip_redistribution.py')
    assert prep['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/prepare_avl_tip_redistribution.py')
    assert prep['protocol_sha256']==raw['protocol_sha256']==obs['protocol_sha256']==digest(ROOT/prep['protocol'])
    assert raw['execution_protocol_sha256']==obs['execution_protocol_sha256']==digest(ROOT/raw['execution_protocol'])
    assert raw['preparation_sha256']==digest(F/'avl-tip-redistribution-preparation.json')
    for n,h in raw['gate_sha256'].items():assert digest(F/n)==h and read(n)['all_criteria_pass']
    assert watch['watchdog_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/avl_tip_watchdog.py')
    assert watch['termination_source_sha256']==raw['dependency_sha256']['avl_span1280_watchdog.py']
    assert watch['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/check_avl_tip_watchdog.py') and len(watch['cases'])==10 and watch['start_cases']==3 and watch['real_termination']
    for data in [raw,obs]:
        for n,h in data['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/n)==h
    for n,h in obs['parent_sha256'].items():assert digest(F/n)==h
    assert digest(ROOT/'build/reference-avl/source-symmetry16000-ks1000-spacing-pi64/bin/avl')==raw['executable_sha256']==prep['executable_sha256']
    assert raw['capacity_preparation_sha256']==digest(F/'avl-span1280-capacity-ks-preparation.json')
    capacity=read('avl-span1280-capacity-ks-preparation.json')
    for n,h in capacity['source_files'].items():assert digest(ROOT/capacity['target']/'src'/n)==h['capacity']
    assert raw['completed'] and raw['stop_reason'] is None and len(raw['runs'])==4 and raw['seconds']<=1200
    expected={(case,allocation) for case in ['FW-V01','FW-V02'] for allocation in ['tip-rich','tip-poor']}
    assert {(r['case'],r['allocation']) for r in raw['runs']}==expected
    total=0;members=0
    for r in raw['runs']:
        assert r['exit_code']==0 and r['stopped'] is None and r['seconds']<=300
        assert max(r['peak_resident_bytes'],r['maximum_group_live_bytes'])<=5*GIB
        assert r['memory_before']['cgroup_known'] and r['memory_before']['effective_available']>=8*GIB
        assert r['minimum_effective_available']>=2*GIB and r['memory_after']['cgroup_known']
        assert r['maximum_swap_used']-min(r['memory_before']['swap_used'],raw['stage_memory_before']['swap_used'])<=256*1024**2
        folder=ROOT/r['directory'];assert digest(folder/'wing.avl')==r['sha256']==digest(ROOT/r['path'])
        assert set(r['before_after_equal'])==NAMES and all(r['before_after_equal'].values()) and set(r['archive_members'])==NAMES
        with zipfile.ZipFile(folder/'after.zip') as archive:
            assert set(archive.namelist())==NAMES
            for n,h in r['archive_members'].items():
                data=archive.read(n);assert hashlib.sha256(data).hexdigest()==h['sha256']==digest(folder/'out'/n) and len(data)==h['bytes'];members+=1
        assert not any((folder/'after').iterdir())
        for n,h in r['files'].items():assert digest(folder/n)==h['sha256'] and (folder/n).stat().st_size==h['bytes'];total+=h['bytes']
    assert total==raw['total_bytes'] and total<=128*1024**2
    assert obs['methods_criteria_pass'] and all(obs['checks'].values()) and obs['stop_reason'] is None
    assert len(obs['meshes'])==6 and len(obs['pairs'])==6 and len(obs['responses'])==2
    assert obs['seconds']<120 and obs['peak_bytes']<512*1024**2
    assert obs['richardson_order'] is None and obs['continuum_tail_estimate'] is None and obs['m24_error_budget'] is None
    records={(r['case'],r['allocation']):r for r in obs['meshes']};max_norm_error=0.
    for key,r in records.items():
        assert r['eligible'] and r['residual_relative']<=1e-11 and r['native_unknowns']==7680 and r['physical_panels']==15360
        expected=[64]*10 if key[1]=='uniform' else ([48]*8+[128]*2 if key[1]=='tip-rich' else [72]*8+[32]*2)
        assert r['actual_interval_counts']==expected
        folder=ROOT/r['raw_directory'];lines=(folder.parent/'wing.avl').read_text().splitlines();sections=[list(map(float,lines[j+1].split()[:5])) for j,s in enumerate(lines) if s=='SECTION'];bounds=[s[1] for s in sections];counts=[0]*10
        for line in (folder/'uv-strips.txt').read_text().splitlines()[2:]:
            values=line.split();index=bisect.bisect_right(bounds[1:],(float(values[1])+float(values[2]))/2);counts[index]+=1
        assert counts==expected
        for n,h in r['files'].items():p=ROOT/r['derived_directory']/n;assert digest(p)==h['sha256'] and p.stat().st_size==h['bytes']
    with localcontext() as ctx:
        ctx.prec=70
        def vector(bands,channel):
            if channel=='tau_p1':return [math.sqrt(b['width'])*b['tau0'] for b in bands]+[math.sqrt(b['width']/3)*b['tau1'] for b in bands]
            name={'tau0':'tau0','k0':'kbar','j_line':'jbar','j_point':'jcbar'}[channel]
            return [math.sqrt(b['width'])*b[name] for b in bands]
        def norm(v):return sum((Decimal.from_float(float(x))**2 for x in v),Decimal(0)).sqrt()
        for p in obs['pairs']:
            a=records[(p['case'],p['from_allocation'])]['bands'];b=records[(p['case'],p['to_allocation'])]['bands'];base=records[(p['case'],'uniform')]['bands']
            assert min(p['margins'].values())>=-1e-10 and p['contraction_excess']<=1e-10
            for c,s in p['channels'].items():
                difference=[y-x for x,y in zip(vector(a,c),vector(b,c))];n=float(norm(difference));max_norm_error=max(max_norm_error,abs(n-s['norm']))
                assert max(abs(x-y) for x,y in zip(difference,s['difference']))<=1e-12
                assert abs(n/float(norm(vector(base,c)))-s['relative_to_uniform'])<=1e-12
                assert abs(sum(s['regional_energy'].values())-n*n)<=1e-12
                if n>0:assert abs(sum(s['regional_fractions'].values())-1)<=1e-12
            for c,s in p['scalar_changes'].items():
                delta=[y[c]-x[c] for x,y in zip(a,b)];total_delta=math.fsum(delta);assert abs(total_delta-s['signed']['total'])<=1e-12
                den=max(abs(math.fsum(v[c] for v in base)),.1 if c.startswith('moment') else 0.)
                assert abs(abs(total_delta)/den-s['relative_to_uniform'])<=1e-12
        for response in obs['responses']:
            case=response['case'];base=records[(case,'uniform')]['bands'];rich=records[(case,'tip-rich')]['bands'];poor=records[(case,'tip-poor')]['bands']
            for c,s in response['channels'].items():
                dr=[y-x for x,y in zip(vector(base,c),vector(rich,c))];dp=[y-x for x,y in zip(vector(base,c),vector(poor,c))]
                den=float(norm(dr)*norm(dp))
                expected=math.fsum(x*y for x,y in zip(dr,dp))/den if den else None
                assert (expected is None and s['direction_cosine'] is None) or abs(expected-s['direction_cosine'])<=1e-12
    assert max_norm_error<=1e-12
    components=read('avl-tip-couple-component-result.json')
    assert components['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/check_avl_tip_couple_components.py')
    assert components['parent_sha256']==digest(F/'avl-tip-redistribution-result.json') and components['protocol_sha256']==digest(F/'avl-tip-couple-diagnostic-addendum.md')
    assert components['all_criteria_pass'] and len(components['pairs'])==6 and components['maximum_reconstruction_error']<=1e-12
    assert components['seconds']<60 and components['peak_bytes']<128*1024**2 and (F/'avl-tip-couple-component-result.json').stat().st_size<1024**2
    for n,h in components['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/n)==h
    for c in components['pairs']:
        pair=next(p for p in obs['pairs'] if (p['case'],p['from_allocation'],p['to_allocation'])==(c['case'],c['from_allocation'],c['to_allocation']))
        a=records[(c['case'],c['from_allocation'])]['bands'];b=records[(c['case'],c['to_allocation'])]['bands']
        couple=[-(y['Cy']-x['Cy'])/(2*math.sqrt(x['width'])) for x,y in zip(a,b)]
        assert max(abs(x-y) for x,y in zip(couple,c['couple_difference']))<=1e-12
        assert max(abs(x+y-z) for x,y,z in zip(c['line_difference'],couple,c['point_difference']))<=1e-12
        line2=math.fsum(v*v for v in c['line_difference']);couple2=math.fsum(v*v for v in couple);cross=2*math.fsum(x*y for x,y in zip(c['line_difference'],couple))
        assert abs(line2+couple2+cross-c['squared_terms']['total'])<=1e-12
        assert abs(c['scalar_point']-pair['scalar_changes']['moment_point']['signed']['total'])<=1e-12
        assert abs(c['scalar_line']+c['scalar_couple']-c['scalar_point'])<=1e-12
    actual=sum(p.stat().st_size for p in (ROOT/raw['root']).rglob('*') if p.is_file());assert actual==obs['combined_output_bytes'] and actual<=128*1024**2
    assert not raw['reference_acceptance'] and not obs['reference_acceptance']
    report=dict(source_sha256=digest(__file__),runs_sha256=digest(F/'avl-tip-redistribution-runs.json'),observations_sha256=digest(F/'avl-tip-redistribution-result.json'),
        completed_runs=4,fields=6,stored_field_pairs=6,archive_members_verified=members,decimal_norm_maximum_difference=max_norm_error,
        run_seconds=raw['seconds'],peak_resident_bytes=max(r['peak_resident_bytes'] for r in raw['runs']),raw_archived_bytes=raw['total_bytes'],combined_bytes=actual,
        couple_component_sha256=digest(F/'avl-tip-couple-component-result.json'),reference_acceptance=False,all_criteria_pass=True)
    (F/'avl-tip-redistribution-artifact-audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
