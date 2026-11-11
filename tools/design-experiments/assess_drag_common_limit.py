from pathlib import Path
import json,hashlib,numpy as np
from uvlm_steady_reference import geometry
from uvlm_frequency_fast import bound,segment
from uvlm_drag_force_kernel import matrix
from probe_drag_far_series import trace,circulation
from trefftz_observations import midpoint
ROOT=Path(__file__).resolve().parents[2];OLD=ROOT/'build/analysis/finite-wing-2-3';FAR=ROOT/'build/analysis/drag-common-limit';NEAR=ROOT/'build/analysis/drag-common-limit-v2'
USED=set()
def use(p):USED.add(p);return p
def far(m,n,d,core=1e-10):
    p=use(FAR/f'reference-d{d}-m{m}-n{n}-core{core}.npz');raw=np.load(p);return float(midpoint(trace(n,d),circulation(raw,m,n)))
def main():
    records=[];series=json.loads(use(FAR/'far-reference-series.json').read_text())['records']
    for d in [0,5]:
        values=np.array([far(8,n,d) for n in [96,192,384,768]]);dv=np.diff(values);assert np.all(dv>0);orders=np.log2(dv[:-1]/dv[1:]);assert min(orders)>.5 and abs(np.diff(orders)[0])<.25;limits=values[2:]+dv[1:]/(2**orders-1);chord=far(16,384,d)-values[2];limit=limits[-1]+chord
        def prior(m,n):
            r=next(r for r in series if r['dihedral']==d and r['m']==m and r['n']==n and r['core']==1e-8);use(ROOT/r['source']);return r['far_drag']
        cd=np.array([prior(32,n)-prior(16,n) for n in [96,192]]);parts=dict(width_extrapolation=1.5*abs(np.diff(limits)[0]),chord_remainder=3*max(abs(cd)),mixed_chord_span=abs(np.diff(cd)[0]),core_fine=abs(far(8,768,d,1e-9)-values[-1]),core_chord=abs(prior(16,384)-far(16,384,d)));unc=sum(parts.values())
        path=use(NEAR/f'candidate-d{d}-m128-n192.json');raw=json.loads(path.read_text());gamma=np.load(use(path.with_suffix('.npz')))['gamma'];g=geometry(128,192,d);n=192;dl=g['right']-g['left'];p=(g['left']+g['right'])/2;N=np.column_stack([np.zeros(len(dl)),-dl[:,2],dl[:,1]]);horse=(gamma.reshape(128,n)-np.vstack([np.zeros((1,n)),gamma.reshape(128,n)[:-1]])).ravel();eps=1.114285714285714e-5
        correction=1.2*sum(horse[j]*gamma[j-n]*segment(p[j],g['right'][j],g['left'][j],eps,N[j])[0] for j in range(n,len(gamma)));near=raw['near_drag']+correction;tail=g['X'].copy();tail[-1,:,0]+=1024;chosen=np.linspace(0,len(gamma)-1,32,dtype=int);C=bound(p[chosen],N[chosen],tail,eps)
        for k,j in enumerate(chosen):
            C[k,j]-=segment(p[j],g['left'][j],g['right'][j],eps,N[j])[0]
            if j>=n:C[k,j-n]-=segment(p[j],g['right'][j],g['left'][j],eps,N[j])[0]
        fresh=matrix(p[chosen],N[chosen],tail,eps,chosen);kernel_difference=float(np.linalg.norm(C-fresh)/np.linalg.norm(fresh));assert kernel_difference<1e-12
        comparison=abs(near-limit)/abs(limit);combined=comparison+unc/abs(limit);nearfar=abs(near/raw['trefftz_midpoint_drag']-1);passed=combined<.01 and nearfar<.01 and max(raw['residual'],raw['uncondensed_sample_residual'])<1e-11 and raw['peak_bytes']<24*2**30 and raw['seconds']<1200
        records.append(dict(dihedral=d,reference_finite_width_series=values.tolist(),width_orders=orders.tolist(),neighbor_width_limits=limits.tolist(),chord_correction=chord,reference_limit=limit,reference_absolute_uncertainty_components=parts,reference_relative_uncertainty=unc/abs(limit),candidate_raw_near=raw['near_drag'],material_self_roundoff_correction=float(correction),candidate_near=float(near),candidate_trefftz=raw['trefftz_midpoint_drag'],same_width_near_far=nearfar,candidate_to_reference_limit=comparison,combined_relative_indicator=combined,material_self_sample_matrix_difference=kernel_difference,passed=bool(passed)))
    USED.update([Path(__file__),ROOT/'tools/design-experiments/run_drag_final_chord.py',ROOT/'tools/design-experiments/run_steady_drag_common_limit_v2.py',ROOT/'tools/design-experiments/run_drag_reference_convergence.py',ROOT/'tools/design-experiments/uvlm_drag_force_kernel.py',ROOT/'docs/verification/drag-common-limit-assessment-policy.md',ROOT/'docs/verification/induced-drag-common-limit-protocol.md',ROOT/'docs/verification/drag-reference-continuum-protocol.md',ROOT/'docs/verification/drag-material-self-result.json',ROOT/'docs/verification/trefftz-midpoint-operator-result.json'])
    result=dict(passed=all(r['passed'] for r in records),scope='FW-V01 and FW-V02 small-angle second-order induced drag coefficient D/alpha^2; midpoint near and far common limit',reference_uncertainty_is_rigorous_bound=False,old_log_energy_comparison_remains_failed=True,whole_uvlm_complete=False,records=records,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(USED)},artifacts={})
    dest=ROOT/'docs/verification/drag-common-limit-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['sources','artifacts']},indent=2));assert result['passed']
if __name__=='__main__':main()
