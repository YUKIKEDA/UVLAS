"""Independent recurrence, partition, provenance and next-input audit; no new solves."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
from decimal import Decimal,localcontext
from fractions import Fraction as Q
import hashlib,json,math,xml.etree.ElementTree as ET
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(n):return json.loads((F/n).read_text())
def main():
    r=read('avl-span-order-diagnostic-result.json');parent=read('avl-span1280-observation-result.json');fig=read('avl-span-order-figures.json');prep=read('avl-tip-redistribution-preparation.json')
    assert r['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/check_avl_span_order_diagnostic.py')
    assert r['protocol_sha256']==digest(ROOT/r['protocol']) and r['parent_sha256']==digest(F/'avl-span1280-observation-result.json')
    for n,h in r['dependency_sha256'].items():assert digest(rd_paths.EXPERIMENTS/n)==h
    assert r['all_criteria_pass'] and len(r['checks'])==8 and all(r['checks'].values())
    assert len(r['manufacturing']['cases'])==9 and all(c['passed'] for c in r['manufacturing']['cases'])
    assert len(r['manufacturing']['refusals'])==4 and all(c['refused'] and c['input_unchanged'] for c in r['manufacturing']['refusals'])
    assert r['seconds']<120 and r['peak_bytes']<512*1024**2 and r['new_avl_runs']==0 and r['continuum_tail_estimate'] is None and not r['reference_acceptance']
    recurrences={'h_h2':(Q(-1,8),Q(3,4)),'h2_h3':(Q(-1,32),Q(3,8)),'h2_h4':(Q(-1,64),Q(5,16)),'h2_h2log':(Q(-1,16),Q(1,2))}
    max_error=0.;tipshares=[];cancellations=[]
    with localcontext() as ctx:
        ctx.prec=70
        def D(x):return Decimal.from_float(float(x))
        def norm(values):return sum((D(x)**2 for x in values),Decimal(0)).sqrt()
        for c in r['channels']:
            assert c['region_counts']==dict(root=4,tip=4,interior=12)
            v=c['values'];diff=[[D(b)-D(a) for a,b in zip(v[j],v[j+1])] for j in range(3)]
            saved=next(s['channels'][c['channel']] for s in parent['series'] if s['case']==c['case'] and s['label']=='span1280')
            for j,stage in enumerate(c['stages']):
                exact_norm=sum((x*x for x in diff[j]),Decimal(0)).sqrt();max_error=max(max_error,abs(float(exact_norm)-saved['differences'][j]))
                assert abs(sum(stage['regional_energy'].values())-stage['norm']**2)<1e-12
                assert abs(sum(stage['regional_fractions'].values())-1)<1e-12
            tipshares.append(c['stages'][-1]['regional_fractions']['tip'])
            for m in c['models']:
                w=recurrences[m['model']];assert [Q(x) for x in m['rational_weights']]==list(w)
                pred=[Decimal(w[0].numerator)/w[0].denominator*a+Decimal(w[1].numerator)/w[1].denominator*b for a,b in zip(diff[0],diff[1])]
                for x,y in zip(pred,m['predicted_last_difference']):max_error=max(max_error,abs(float(x)-y))
                error=sum(((x-y)**2 for x,y in zip(pred,diff[2])),Decimal(0)).sqrt()/sum((x*x for x in diff[2]),Decimal(0)).sqrt()
                assert abs(float(error)-m['relative_residual'])<1e-10
        for s in r['scalars']:
            values=[D(x) for x in s['values']];diff=[values[j+1]-values[j] for j in range(3)]
            for m in s['models']:
                w=recurrences[m['model']]
                prediction=Decimal(w[0].numerator)/w[0].denominator*diff[0]+Decimal(w[1].numerator)/w[1].denominator*diff[1]
                assert abs(float(prediction)-m['predicted_last_difference'][0])<1e-12
                assert abs(float(abs(prediction-diff[2])/abs(diff[2]))-m['relative_residual'])<1e-10
            for stage in s['stages']:
                values=stage['bands'];total=math.fsum(values);absolute=math.fsum(abs(v) for v in values)
                assert abs(total-stage['signed']['total'])<1e-12
                assert abs(absolute-stage['signed']['absolute_sum'])<1e-12
                cancellations.append(stage['signed']['cancellation_ratio'])
                if 'physical_parts' in stage:assert abs(math.fsum(v['total'] for v in stage['physical_parts'].values())-total)<1e-12
    assert max_error<1e-12
    assert fig['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/plot_avl_span_order_diagnostic.py') and fig['parent_sha256']==digest(F/'avl-span-order-diagnostic-result.json')
    for n,h in fig['files'].items():assert digest(F/'figures'/n)==h;ET.parse(F/'figures'/n)
    assert all(0<=y<=.5 and -1<=x<=1 for series in fig['data']['normalized_tau0'] for x,y in series)
    for j in range(10):assert abs(sum(row[j] for row in fig['data']['regional_percentages'])-100)<1e-10
    assert prep['source_sha256']==digest(ROOT/'rd/aero/avl-reference/experiments/prepare_avl_tip_redistribution.py')
    assert prep['protocol_sha256']==digest(ROOT/prep['protocol']) and prep['parent_sha256']==digest(F/'avl-span1280-fine-ks-runs.json')
    assert prep['criteria_pass'] and len(prep['decks'])==4 and prep['new_avl_runs']==0 and not prep['reference_acceptance']
    assert digest(ROOT/'build/reference-avl/source-symmetry16000-ks1000-spacing-pi64/bin/avl')==prep['executable_sha256']
    for d in prep['decks']:
        source=ROOT/d['donor'];target=ROOT/d['path'];assert digest(source)==d['donor_sha256'] and digest(target)==d['sha256']
        a=source.read_bytes().splitlines();b=target.read_bytes().splitlines();indices=[j+1 for j,s in enumerate(a) if s==b'SECTION']
        expected=[48]*8+[128]*2 if d['allocation']=='tip-rich' else [72]*8+[32]*2
        assert expected==d['interval_counts'] and sum(expected)==640 and d['native_unknowns']==7680
        assert [j for j,(x,y) in enumerate(zip(a,b)) if x!=y]==indices[:10]
        for j,n in zip(indices[:10],expected):assert a[j].split()[:5]==b[j].split()[:5] and b[j].split()[-2:]==[str(n).encode(),b'1.0']
        assert a[indices[-1]]==b[indices[-1]]
    output_paths=[F/'avl-span-order-diagnostic-result.json',F/'avl-span-order-figures.json',*(F/'figures'/n for n in fig['files'])]
    total_bytes=sum(p.stat().st_size for p in output_paths);assert total_bytes<16*1024**2
    result=dict(source_sha256=digest(__file__),diagnostic_sha256=digest(F/'avl-span-order-diagnostic-result.json'),figures_sha256=digest(F/'avl-span-order-figures.json'),
        preparation_sha256=digest(F/'avl-tip-redistribution-preparation.json'),channels=10,scalar_series=6,manufactured_cases=9,refusals=4,
        independent_decimal_maximum_difference=max_error,tip_energy_range=[min(tipshares),max(tipshares)],
        scalar_spatial_cancellation_range=[min(cancellations),max(cancellations)],next_inputs=4,new_avl_runs=0,output_bytes=total_bytes,
        reference_acceptance=False,all_criteria_pass=True)
    (F/'avl-span-order-diagnostic-artifact-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
