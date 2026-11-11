from pathlib import Path
import json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/finite-wing-2-3'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    records=[];artifacts={};failures=[]
    def get(side,d,key):
        p=OUT/f'{side}-d{d}-{key}.json';artifacts[p.relative_to(ROOT).as_posix()]=sha(p)
        if p.with_suffix('.npz').exists():artifacts[p.with_suffix('.npz').relative_to(ROOT).as_posix()]=sha(p.with_suffix('.npz'))
        return json.loads(p.read_text())
    for d in [0,5]:
        data={}
        for side in ['candidate','reference']:
            a=get(side,d,'band-n192');b=get(side,d,'band-n384');c=get(side,d,'band-n768');chord=get(side,d,'span384');data[side]=c;span=abs(np.array(b['bands'])/c['bands']-1);dx=abs(np.array(chord['bands'])/b['bands']-1);indicator=span+dx;contraction=abs((np.array(b['bands'])-c['bands'])/(np.array(a['bands'])-b['bands']));records.append(dict(case='fixed_band_space',side=side,dihedral=d,span_change=span.tolist(),chord_change=dx.tolist(),sum_indicator=indicator.tolist(),contraction=contraction.tolist(),max_indicator=float(max(indicator)),passed=bool(max(indicator)<.01),rigorous_error_bound=False))
            if max(indicator)>=.01:failures.append(side+' band space')
        a=data['candidate'];b=data['reference'];difference=abs(np.array(a['bands'])/b['bands']-1);core0=get('reference',d,'band768-core1e-07');core1=get('reference',d,'band768-core1e-09');sensitivity=max(float(max(abs(np.array(core0['bands'])/b['bands']-1))),float(max(abs(np.array(core1['bands'])/b['bands']-1))));combined=np.array(records[-2]['sum_indicator'])+np.array(records[-1]['sum_indicator'])+difference+sensitivity;records.append(dict(case='fixed_band_comparison_and_combined_indicator',dihedral=d,max_difference=float(max(difference)),reference_core_max_sensitivity=sensitivity,combined_indicators=combined.tolist(),max_combined_indicator=float(max(combined)),passed=bool(max(combined)<.02)))
        if max(combined)>=.02:failures.append('combined band indicator')
    sources={p.relative_to(ROOT).as_posix():sha(p) for p in [Path(__file__),ROOT/'tools/design-experiments/run_fixed_band_extension.py',ROOT/'tools/design-experiments/probe_fixed_band_reference_core.py',ROOT/'tools/design-experiments/uvlm_frequency_fast.py',ROOT/'docs/verification/fixed-band-final-extension-protocol.md',ROOT/'docs/verification/finite-wing-extended-campaign-result.json']}
    result=dict(passed=not failures,records=records,failures=failures,sources=sources,artifacts=artifacts,uvlm_2_3_complete=False);dest=ROOT/'docs/verification/fixed-band-final-extension-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(passed=result['passed'],summary=[{k:v for k,v in r.items() if k not in ['span_change','chord_change','sum_indicator','contraction','combined_indicators']} for r in records],failures=failures),indent=2))
if __name__=='__main__':main()
