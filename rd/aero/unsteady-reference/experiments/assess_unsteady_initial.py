"""Keep the failed finite campaign and quantify independent frequency validation."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,numpy as np
from extract_unsteady_reference import extract
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/unsteady-2-3'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def value(label):return np.array([complex(*v) for v in json.loads((OUT/(label+'.json')).read_text())['transfer']])
def main():
    artifacts={};sources={};records=[]
    raw0=np.load(OUT/'original.npz');raw1=np.load(OUT/'original-core-exact-parity.npz');assert np.array_equal(raw0['F'],raw1['F']) and np.array_equal(raw0['M'],raw1['M']);records.append(dict(case='unmodified_core_parameter_parity',bitwise_loads=True,uses_actual_computed_original_core=True,rounded_initial_fixture_not_bitwise=True))
    ref=value('ptera-frequency-m2-n8');times=[]
    for s in [64,128,256,512,1024]:
        m,c,fit=extract(f'time{s}');times.append(c);records.append(dict(case='original_time_to_frequency',samples=s,relative_difference=abs((c-ref)/ref).tolist(),fit_residual=fit.tolist()))
    ext0=2*times[-2]-times[-3];ext1=2*times[-1]-times[-2];records.append(dict(case='independent_first_order_time_extrapolation',previous_difference=abs((ext0-ref)/ref).tolist(),last_difference=abs((ext1-ref)/ref).tolist(),extrapolation_change=abs((ext1-ext0)/ref).tolist(),raw_time_qualified=False))
    for label in ['cycles6','amplitude1e-5','nu1e-13','core1e-3','core1e-4']:
        _,c,_=extract(label);records.append(dict(case='reference_sensitivity',label=label,relative_difference=abs((c-times[2])/times[2]).tolist()))
    for m,n in [(2,8),(4,24),(8,48),(16,48),(16,96)]:
        a=value(f'candidate-m{m}-n{n}');b=value(f'ptera-frequency-m{m}-n{n}');records.append(dict(case='frequency_spatial_comparison',m=m,n=n,relative_difference=abs((a-b)/b).tolist()))
    selected=['original','original-core-parity','original-core-exact-parity','time64','time128','time256','time512','time1024','cycles6','core1e-3','core1e-4','amplitude1e-5','nu1e-13','space4x24','space8x48']
    selected += [f'{side}-m{m}-n{n}' for side in ['candidate','ptera-frequency'] for m,n in [(2,8),(4,24),(8,48),(16,48),(16,96)]]
    selected += [side+'-'+label for side in ['candidate','ptera-frequency'] for label in ['quadrature8','length128','length512','core1e-4']]
    for label in selected:
        for suffix in ['.json','.npz']:
            p=OUT/(label+suffix);artifacts[p.relative_to(ROOT).as_posix()]=sha(p)
    for name in ['assess_unsteady_initial.py','probe_ptera_reference.py','check_ptera_exact_core_parity.py','run_unsteady_reference_campaign.py','run_unsteady_candidate_campaign.py','run_ptera_frequency_campaign.py','run_ptera_time_extension.py','ptera_frequency_reference.py','uvlm_harmonic_reference_v2.py','extract_unsteady_reference.py']:
        p=rd_paths.EXPERIMENTS/name;sources[p.relative_to(ROOT).as_posix()]=sha(p)
    for name in ['unsteady-independent-reference-protocol.md','unsteady-frequency-reference-protocol.md','harmonic-age-integral-result.json','frequency-acceleration-result.json']:
        p=rd_paths.RECORDS/name;sources[p.relative_to(ROOT).as_posix()]=sha(p)
    env=ROOT/'build/reference-ptera/environment-manifest.json';sources[env.relative_to(ROOT).as_posix()]=sha(env);sources.update(json.loads(env.read_text())['sources'])
    result=dict(passed=False,failures=['raw time discretization exceeds 0.20%','M16/N96 frequency My difference exceeds 1%','full spatial budget not yet qualified'],records=records,sources=sources,artifacts=artifacts,uvlm_2_3_complete=False,exploratory_moment_pole_sign_corrected_from_raw=True)
    dest=ROOT/'rd/aero/unsteady-reference/records/unsteady-initial-campaign-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(passed=False,records=records),indent=2))
if __name__=='__main__':main()
