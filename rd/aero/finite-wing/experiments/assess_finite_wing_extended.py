"""Bounded reference disposition; global and distribution evidence are distinguished."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,numpy as np
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/finite-wing-2-3'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def vector(a):return np.array([a['F'][2],a['M'][1]])
def main():
    records=[];sources={};artifacts={}
    def get(side,d,key):
        p=OUT/f'{side}-d{d}-{key}.json';artifacts[p.relative_to(ROOT).as_posix()]=sha(p);snap=p.with_suffix('.npz');artifacts[snap.relative_to(ROOT).as_posix()]=sha(snap);return json.loads(p.read_text())
    for d in [0,5]:
        for side in ['candidate','reference']:
            fine=get(side,d,'span384');base=get(side,d,'m16-n192');coarse=get(side,d,'m16-n96');chord=get(side,d,'m32-n192');chord0=get(side,d,'m32-n96')
            ds=abs(vector(base)/vector(fine)-1);dc=np.maximum(abs(vector(chord)/vector(base)-1),abs(vector(chord0)/vector(coarse)-1));indicator=ds+dc;contraction=abs((vector(base)-vector(fine))/(vector(coarse)-vector(base)))
            bands=abs(np.array(base['bands'])/fine['bands']-1);band_chord=abs(np.array(chord['bands'])/base['bands']-1)
            records.append(dict(case='global_space_estimate',side=side,dihedral=d,span_difference=ds.tolist(),chord_difference=dc.tolist(),sum_indicator=indicator.tolist(),contraction=contraction.tolist(),passed=bool(np.max(indicator)<.003),band_max_span_difference=float(max(bands)),band_max_chord_difference=float(max(band_chord)),band_l2_span_difference=float(np.linalg.norm(np.array(base['bands'])-fine['bands'])/np.linalg.norm(fine['bands'])),rigorous_error_bound=False))
        own=get('candidate',d,'span384');ref=get('reference',d,'span384');base=get('candidate',d,'m16-n192');tail0=get('candidate',d,'length256.0');tail1=get('candidate',d,'length512.0');core1=get('candidate',d,'core0.0001');core0=get('candidate',d,'core0.001');rb=get('reference',d,'m16-n192');rc=get('reference',d,'core1e-7')
        diff=abs(vector(own)/vector(ref)-1);bdiff=float(max(abs(np.array(own['bands'])/ref['bands']-1)));tail=abs(vector(tail1)/vector(base)-1);tailratio=abs((vector(tail1)-vector(base))/(vector(tail0)-vector(tail1)));refcore=abs(vector(rc)/vector(rb)-1);owncore=abs(vector(core1)/vector(base)-1)
        # Actual signed angle runs check the analytic first-order circulation reconstruction.
        odd=[]
        for angle in [1.,2.]:
            plus=get('reference',d,'m4-n12') if angle==1 else get('reference',d,'angle2.0');minus=get('reference',d,f'angle{-angle}');rad=np.deg2rad(angle);raw=np.array([(plus['raw_F'][2]-minus['raw_F'][2])/(2*rad),(plus['raw_M'][1]-minus['raw_M'][1])/(2*rad)]);prediction=vector(plus)*np.sin(rad)*np.cos(rad)/rad;err=abs(raw/prediction-1);assert max(err)<1e-10;odd.append(raw)
        angleerror=abs(odd[1]/odd[0]-1)
        assert max(diff)<.01 and bdiff<.02 and max(refcore)<.0005 and max(tail)<.0015 and max(angleerror)<.001
        records.append(dict(case='independent_global_and_20_band_comparison',dihedral=d,global_difference=diff.tolist(),band_max_difference=bdiff,tail_512_to_1024=tail.tolist(),tail_contraction=tailratio.tolist(),reference_core_sensitivity=refcore.tolist(),own_core_1e4_to_1e5_sensitivity=owncore.tolist(),angle_1_to_2_sensitivity=angleerror.tolist(),reconstructed_odd_angles_verified=True,passed=True))
    for p in [Path(__file__),ROOT/'rd/aero/finite-wing/experiments/run_finite_wing_span_extension.py',ROOT/'rd/aero/finite-wing/records/finite-wing-span-extension-protocol.md',ROOT/'rd/aero/finite-wing/records/finite-wing-initial-campaign-result.json',ROOT/'rd/aero/finite-wing/records/finite-wing-specialization-result.json',ROOT/'build/reference-aerosandbox/environment-manifest.json']:
        sources[p.relative_to(ROOT).as_posix()]=sha(p)
    env=json.loads((ROOT/'build/reference-aerosandbox/environment-manifest.json').read_text());sources.update(env['sources'])
    assert all(r['passed'] for r in records)
    result=dict(passed=True,scope='two finite swept tapered wings, zero-angle total lift and My derivatives; conservative 20-band comparison',records=records,sources=sources,artifacts=artifacts,global_spatial_estimates_within_budget=True,individual_band_continuum_uncertainty_qualified=False,unsteady_reference_qualified=False,induced_drag_qualified=False,uvlm_2_3_complete=False)
    dest=ROOT/'rd/aero/finite-wing/records/finite-wing-extended-campaign-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['sources','artifacts']},indent=2))
if __name__=='__main__':main()
