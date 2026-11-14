"""Independent geometry-column, source, state, archive and sensitivity audit."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from decimal import Decimal,localcontext
import hashlib,json,math,re,shlex,zipfile
import numpy as np
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS;GIB=1024**3
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(n):return json.loads((F/n).read_text())
def main():
    prep=read('avl-symmetric-solve-preparation.json');w=read('avl-symmetric-solve-weight-result.json');small=read('avl-symmetric-solve-small-check.json');out=read('avl-symmetric-solve-sensitivity-result.json')
    assert prep['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/prepare_avl_symmetric_solve_probe.py')
    assert w['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/prepare_avl_symmetric_solve_weights.py')
    assert small['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/check_avl_symmetric_solve_small.py')
    assert out['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/check_avl_symmetric_solve_sensitivity.py')
    assert prep['weights_sha256']==digest(F/'avl-symmetric-solve-weight-result.json')
    assert prep['protocol_sha256']==w['protocol_sha256']==out['protocol_sha256']==digest(ROOT/prep['protocol'])
    for d in [w,out]:
        for n,h in d['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/n)==h
    for n,h in w['parent_sha256'].items():assert digest(F/n)==h
    for n,h in out['parent_sha256'].items():assert digest(F/n)==h
    initial=read('avl-symmetric-solve-weight-initial-result.json')
    for n,h in initial['sources'].items():assert digest(ROOT/n)==h
    assert initial['manufacturing_error'] and initial['reference_unchanged']
    assert w['all_criteria_pass'] and w['manufacturing']['all_criteria_pass'] and len(w['records'])==8 and w['seconds']<600 and w['peak_bytes']<512*1024**2 and w['output_bytes']<32*1024**2
    parent=ROOT/prep['parent'];target=ROOT/prep['target'];changed=[]
    for n,h in prep['source_files'].items():
        assert digest(parent/'src'/n)==h['parent'] and digest(target/'src'/n)==h['probe']
        if h['parent']!=h['probe']:changed.append(n)
    assert changed==['aoper.f'] and digest(parent/'bin/avl')==prep['parent_executable_sha256']
    assert digest(target/'src'/prep['new_adapter'])==prep['adapter_sha256']==digest(ROOT/'build/third-party-gpl/avl/reference_avl_sym_solve_probe.f')
    assert len(prep['reused_object_sha256'])==40
    for n,h in prep['reused_object_sha256'].items():assert digest(parent/'bin'/n)==h==digest(target/'bin'/n)
    linux='<repo>';build=ROOT/'build/reference-avl/build-symmetric-solve-probe.sh'
    assert digest(build)==prep['build_script_sha256'] and digest(F/'fixtures/avl-measurement/symmetric-solve-probe.patch')==prep['patch_sha256']
    commands=[shlex.split(s) for s in build.read_text().splitlines() if ' -c ' in s];assert len(commands)==2
    for c in commands:assert '-fdefault-real-8' in c and digest(ROOT/c[0][len(linux)+1:])==prep['compiler_sha256']
    for n,h in prep['link_inputs'].items():assert digest(ROOT/n[len(linux)+1:])==h
    generation=read('avl-symmetric-solve-runner-preparation.json')
    assert generation['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/prepare_avl_symmetric_solve_runner.py')
    assert generation['template_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/run_avl_tip_redistribution.py')
    assert generation['runner_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/run_avl_symmetric_solve.py')
    column_error=0.;sample_columns=0
    # Closed-form independent columns using Decimal80, not the weight generator.
    with localcontext() as ctx:
        ctx.prec=80
        def D(x):return Decimal.from_float(float(x))
        bounds=[D(-4+.4*j) for j in range(21)];slope=D(math.tan(math.radians(35)))-D(.225)
        for rec in w['records']:
            assert rec['all_criteria_pass'] and max(rec['original_reproduction_error'],*(max(t[n] for n in ['direct_error','fold_error','symmetry_error','transpose_error']) for t in rec['trials']))<=1e-12
            geometry=ROOT/rec['geometry_directory']
            for n,h in rec['geometry_files'].items():assert digest(geometry/n)==h
            assert digest(ROOT/rec['weight_path'])==rec['weight_sha256'];cache=np.load(ROOT/rec['weight_path']);matrix=cache['matrix']
            n=rec['native_unknowns'];m=rec['m'];assert matrix.shape==(100,n)
            rows=np.loadtxt(geometry/'uv-rows.txt',skiprows=2);forces=np.loadtxt(geometry/'uv-forces.txt',skiprows=2);points=forces[forces[:,2]==0,3:6]
            def column(index):
                aa=rows[index,7:10];bb=rows[index,10:13];p=points[index];lo=D(aa[1]);hi=D(bb[1]);dy=hi-lo
                fx=(D(bb[0])-D(aa[0]))/dy;intercept=D(aa[0])-fx*lo
                cp=-float(p[0]-(aa[0]+bb[0])/2)*(2*float(bb[1]-aa[1]));answer=[Decimal(0)]*100
                for j in range(20):
                    left=max(lo,bounds[j]);right=min(hi,bounds[j+1])
                    if right<=left:continue
                    width=bounds[j+1]-bounds[j];center=(bounds[j]+bounds[j+1])/2;te_slope=slope if center>0 else -slope
                    length=right-left;t1=((right-center)**2-(left-center)**2)/2
                    kval=(D(1.5)-intercept)*length+(te_slope-fx)*(right*right-left*left)/2
                    line=(intercept-D(.375))*length+fx*(right*right-left*left)/2
                    values=[length/width,6*t1/(width*width),kval/width,line/width,(line-D(cp)*length/(2*dy))/width]
                    for k,v in enumerate(values):answer[20*k+j]=v
                return answer
            indices=sorted({0,1,n//4,n//2,3*n//4,n-2,n-1})
            for i in indices:
                negative=(n//m-1-i//m)*m+i%m;positive=n+i;a=column(negative);b=column(positive)
                error=max(abs(float(x+y)-v)/max(1.,abs(float(x+y))) for x,y,v in zip(a,b,matrix[:,i]));column_error=max(column_error,error);sample_columns+=1
    assert column_error<=1e-12
    stage_results=[];raw_bytes=0;archive_members=0
    for stage,count in [('small',2),('main',6)]:
        raw=read(f'avl-symmetric-solve-{stage}-runs.json');assert raw['completed'] and raw['stop_reason'] is None and raw['stage']==stage and len(raw['runs'])==count
        assert raw['source_sha256']==generation['runner_sha256'] and raw['preparation_sha256']==digest(F/'avl-symmetric-solve-preparation.json')
        assert raw['executable_sha256']==digest(target/'bin/avl')
        for n,h in raw['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/n)==h
        for n,h in raw['gate_sha256'].items():assert digest(F/n)==h and read(n)['all_criteria_pass']
        assert raw['build_script_sha256']==digest(build)
        log=ROOT/'build/reference-avl/build-symmetric-solve-probe.log';timing=ROOT/'build/reference-avl/build-symmetric-solve-probe-resource.txt'
        assert raw['build_log_sha256']==digest(log) and raw['build_resource_sha256']==digest(timing)
        text=timing.read_text();assert re.search(r'Exit status:\s*0',text) and int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',text)[1])*1024<2*GIB
        clock=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)',text)[1]
        assert sum(float(v)*60**j for j,v in enumerate(reversed(clock.split(':'))))<300
        per,total,mem=(60,180,2*GIB) if stage=='small' else (300,2400,5*GIB)
        assert raw['seconds']<=total
        expected={(c,a) for c in ['FW-V01','FW-V02'] for a in (['small'] if stage=='small' else ['uniform','tip-rich','tip-poor'])}
        assert {(r['case'],r['allocation']) for r in raw['runs']}==expected
        size=0
        for run in raw['runs']:
            assert run['seconds']<=per and max(run['peak_resident_bytes'],run['maximum_group_live_bytes'])<=mem and run['exit_code']==0 and run['stopped'] is None
            assert run['memory_before']['cgroup_known'] and run['memory_before']['effective_available']>=8*GIB and run['minimum_effective_available']>=2*GIB
            assert run['maximum_swap_used']-min(run['memory_before']['swap_used'],raw['stage_memory_before']['swap_used'])<=256*1024**2
            folder=ROOT/run['directory'];assert digest(folder/'wing.avl')==run['input_sha256']==digest(ROOT/run['input'])
            assert all(run['before_after_equal'].values()) and all(run['parent_outputs_equal'].values()) and len(run['before_after_equal'])==(9 if stage=='small' else 8)
            with zipfile.ZipFile(folder/'after.zip') as archive:
                assert set(archive.namelist())==set(run['archive_members'])
                for n,h in run['archive_members'].items():
                    data=archive.read(n);assert hashlib.sha256(data).hexdigest()==h['sha256']==digest(folder/'out'/n)==digest(ROOT/run['raw_directory']/n) and len(data)==h['bytes'];archive_members+=1
            for n,h in run['files'].items():assert digest(folder/n)==h['sha256'] and (folder/n).stat().st_size==h['bytes'];size+=h['bytes']
            meta=run['probe_meta'];assert meta['N']==run['native_unknowns'] and meta['REAL_BITS']==64 and meta['QUAD_BITS']==128 and meta['QUAD_DIGITS']>=30 and meta['STATE_UNCHANGED']==1
            tests=np.loadtxt(folder/'uq/uq-guard-tests.txt',skiprows=1,dtype=int);assert tests.shape==(8,3) and np.array_equal(tests[:,1]!=0,tests[:,2]!=0)
            cross=np.loadtxt(folder/'uq/uq-residual-cross.txt',skiprows=1);v=np.loadtxt(folder/'uq/uq-vectors.txt',skiprows=1)
            assert v.shape==(run['native_unknowns'],10) and np.all(np.isfinite(v))
            assert cross.shape==((4*run['native_unknowns'] if stage=='small' else 32),5)
            assert np.max(abs(cross[:,2]-cross[:,3])/np.maximum(1.,cross[:,4]))<=1e-12
            if stage=='small':assert digest(folder/'uq/uq-matrix.bin')==digest(folder/'out/uh-matrix.bin')
        assert size==raw['total_bytes'] and size<=128*1024**2;raw_bytes+=size;stage_results.append(dict(stage=stage,cases=count,seconds=raw['seconds'],peak_bytes=max(r['peak_resident_bytes'] for r in raw['runs'])))
    assert raw_bytes<=128*1024**2
    assert small['all_criteria_pass'] and len(small['records'])==2 and small['small_runs_sha256']==digest(F/'avl-symmetric-solve-small-runs.json') and small['weights_module_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/avl_symmetric_observation_weights.py')
    assert out['all_criteria_pass'] and all(out['checks'].values()) and len(out['systems'])==6 and len(out['comparisons'])==4 and out['seconds']<300 and out['peak_bytes']<512*1024**2
    assert out['forward_error_bound'] is None and out['coefficient_generation_error_bound'] is None and out['m24_error_budget'] is None and not out['reference_acceptance']
    systems={(r['case'],r['allocation']):r for r in out['systems']}
    delta_channels={};observation_delta_error=0.
    main_raw=read('avl-symmetric-solve-main-runs.json')
    with localcontext() as ctx:
        ctx.prec=80
        for run in main_raw['runs']:
            key=(run['case'],run['allocation']);system=systems[key];cache=np.load(ROOT/run['weight_path']);matrix=cache['matrix'];widths=cache['widths']
            vectors=np.loadtxt(ROOT/run['directory']/'uq/uq-vectors.txt',skiprows=1);q=vectors[:,2:6]
            sparse=[[(int(i),Decimal.from_float(float(row[i]))) for i in np.flatnonzero(row)] for row in matrix]
            for j,state in enumerate(system['states'],1):
                difference=q[:,j]-q[:,0];d=[Decimal.from_float(float(x)) for x in difference]
                actual=np.array([float(sum((a*d[i] for i,a in row),Decimal(0))) for row in sparse]);saved=np.array(state['observation_change'])
                observation_delta_error=max(observation_delta_error,float(np.max(abs(actual-saved))))
                assert float(np.max(abs(actual-saved)))<=1e-12
                root=np.sqrt(widths)
                channel=dict(tau0=root*actual[:20],tau_p1=np.concatenate([root*actual[:20],np.sqrt(widths/3)*actual[20:40]]),
                    k0=root*actual[40:60],j_line=root*actual[60:80],j_point=root*actual[80:100])
                delta_channels[(key,state['state'])]=channel
                for c,v in channel.items():assert abs(float(np.linalg.norm(v))-state['channels'][c]['norm'])<=1e-12
    for comparison in out['comparisons']:
        case=comparison['case'];allocation=comparison['allocation'];base=systems[(case,'uniform')];target_sys=systems[(case,allocation)]
        for state in comparison['states']:
            for c,v in state['channels'].items():
                a=next(s for s in base['states'] if s['state']==state['state'])['channels'][c]['norm'];b=next(s for s in target_sys['states'] if s['state']==state['state'])['channels'][c]['norm']
                assert abs((a+b)-v['individual_change_norm_sum'])<=1e-12
                assert abs((a+b)/v['original_response_norm']-v['individual_sum_ratio'])<=1e-12
                da=delta_channels[((case,'uniform'),state['state'])][c];db=delta_channels[((case,allocation),state['state'])][c]
                assert abs(float(np.linalg.norm(db-da))-v['pair_response_change_norm'])<=1e-12
                assert abs(float(np.linalg.norm(db-da))/v['original_response_norm']-v['response_change_ratio'])<=1e-12
    result=dict(source_sha256=digest(__file__),sensitivity_sha256=digest(F/'avl-symmetric-solve-sensitivity-result.json'),weight_sha256=digest(F/'avl-symmetric-solve-weight-result.json'),
        small_check_sha256=digest(F/'avl-symmetric-solve-small-check.json'),initial_manufacturing_failure_preserved=True,sampled_native_columns=sample_columns,
        independent_decimal_column_error=column_error,independent_decimal_observation_change_error=observation_delta_error,stages=stage_results,archive_members_verified=archive_members,raw_bytes=raw_bytes,
        all_criteria_pass=True,reference_acceptance=False)
    (F/'avl-symmetric-solve-artifact-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
