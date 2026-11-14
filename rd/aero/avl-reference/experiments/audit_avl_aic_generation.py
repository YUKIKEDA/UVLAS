"""Audit all saved coefficient records and provenance without running AVL."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from decimal import Decimal, localcontext
import gzip
import hashlib
import json
import math
import numpy as np

ROOT=rd_paths.ROOT;FOLDER=rd_paths.RECORDS
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(name):return json.loads((FOLDER/name).read_text())
def finite(v):
    if isinstance(v,float):assert math.isfinite(v)
    elif isinstance(v,dict):
        for x in v.values():finite(x)
    elif isinstance(v,list):
        for x in v:finite(x)

def main():
    result=read('avl-aic-generation-result.json');parent=read('avl-solve-sensitivity-runs.json')
    assert result['source_sha256']==digest(rd_paths.sibling('check_avl_aic_generation.py'))
    assert result['parent_runs_sha256']==digest(FOLDER/'avl-solve-sensitivity-runs.json')
    assert result['protocol_sha256']==digest(FOLDER/'avl-aic-generation-experiment.md')
    for name,value in result['dependency_sha256'].items():assert digest(rd_paths.sibling(name))==value
    for key,name in [('kernel_sha256','aic.f'),('setup_sha256','asetup.f')]:
        assert result[key]==digest(ROOT/'build/reference-avl/source-core0-grid/src'/name)
    assert len(result['records'])==8 and len(result['checks'])==7 and all(result['checks'].values())
    assert result['all_criteria_pass'] and result['stop_reason'] is None
    assert not result['reference_convergence_accepted'] and not result['large_matrix_global_error_bounded']
    assert result['seconds']<300 and result['peak_process_working_bytes']<512*1024**2 and result['total_bytes']<128*1024**2
    count=0;selfaudits=0;quadratures=0;files=0
    for rec,run in zip(result['records'],parent['runs']):
        assert all(rec[k]==run[k] for k in ['family','case','m','n'])
        folder=ROOT/rec['directory'];n=rec['unknowns']
        expected=n*n if n<=192 else 16*n-64
        assert rec['coefficient_count']==expected
        assert rec['diagonal_includes_bound'] and not rec['partial_is_error_bound']
        assert rec['full_matrix']==rec['partial_is_full']==(n<=192)
        for name,value in rec['files'].items():
            p=folder/name;assert digest(p)==value['sha256'] and p.stat().st_size==value['bytes'];files+=1
        with np.load(folder/'coefficient-diagnostics.npz') as archive:
            arrays={name:archive[name] for name in archive.files}
        for a in arrays.values():assert np.all(np.isfinite(a))
        ids=set();seen=0
        with localcontext() as ctx,gzip.open(folder/'coefficients.tsv.gz','rt',encoding='utf-8') as stream:
            ctx.prec=80
            assert stream.readline().startswith('i\tj\tstored_float64')
            for k,line in enumerate(stream):
                i,j,s,v,d=line.split();i=int(i);j=int(j)
                assert (i,j) not in ids;ids.add((i,j))
                assert i==arrays['i'][k] and j==arrays['j'][k]
                assert float(s)==arrays['stored'][k]
                assert Decimal(v)-Decimal.from_float(float(s))==Decimal(d)
                assert float(v)==arrays['reference_projected'][k] and float(d)==arrays['decimal_difference'][k]
                seen+=1
        assert seen==expected and (1,1) in ids
        count+=seen;selfaudits+=len(rec['decimal_self_audits']);quadratures+=len(rec['independent_line_integrals'])
        assert all(v['eligible'] for v in rec['decimal_self_audits']+rec['independent_line_integrals'])
    assert count==296832 and quadratures==16
    initial=read('avl-aic-generation-initial-status.json')
    assert initial['status']=='postprocessing_serialization_failed' and initial['progress']['completed']==8
    assert initial['source_sha256']==digest(ROOT/'build/analysis/check_avl_aic_generation_initial.py')
    for name,value in initial['files'].items():
        p=ROOT/initial['directory']/name;assert digest(p)==value['sha256'] and p.stat().st_size==value['bytes']
    finite(result)
    report=dict(source_sha256=digest(Path(__file__)),result_sha256=digest(FOLDER/'avl-aic-generation-result.json'),
        coefficients=count,decimal160_cases=selfaudits,line_integral_cases=quadratures,
        files_verified=files,preserved_initial_files=len(initial['files']),criteria=7,all_criteria_pass=True)
    (FOLDER/'avl-aic-generation-artifact-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
