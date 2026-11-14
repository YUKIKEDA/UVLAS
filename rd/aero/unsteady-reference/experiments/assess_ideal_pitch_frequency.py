import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,numpy as np
from uvlm_error_certification import estimate
ROOT=rd_paths.ROOT;REF=ROOT/'build/analysis/frequency-ideal-reference';FINE=ROOT/'build/analysis/frequency-delta';Q24=ROOT/'build/analysis/frequency-ideal-reference-q24';OLD=ROOT/'build/analysis/pitch-frequency-2-3';USED=set()
def use(p):USED.add(p);return p
def get(folder,side,m,n,k):
    p=use(folder/f'{side}-m{m}-n{n}-k{k}.npz');use(p.with_suffix('.json'));return np.load(p)['transfer']
def pairs(x):return [float(np.real(x)),float(np.imag(x))]
def main():
    records=[]
    for k in [.05,.2,.8]:
        widths=np.array([get(REF,'reference',16,n,k) for n in [48,96,192,384]]);chords=np.array([get(REF,'reference',m,192,k) for m in [8,16,32,64]]);candidate=get(FINE,'candidate',128,192,k);mixed=abs((chords[2]-chords[1])-(get(REF,'reference',32,96,k)-widths[1]));quad=np.maximum(abs(get(Q24,'reference',16,96,k)-widths[1]),abs(get(Q24,'reference',64,192,k)-chords[-1]));native=get(FINE,'reference',64,192,k)
        for a in range(2):
            for b in range(2):
                w=estimate(widths[:,a,b]);c=estimate(chords[:,a,b]);r=dict(k=k,motion=['heave','pitch'][a],observable=['Fz','My'][b],width_resolved=w['resolved'],chord_resolved=c['resolved'],width_orders=w.get('orders'),chord_orders=c.get('orders'),candidate=pairs(candidate[a,b]),passed=False,reference_native_core_family_change_at_m64n192=float(abs(native[a,b]-chords[-1,a,b])/abs(chords[-1,a,b])))
                if w['resolved'] and c['resolved']:
                    limit=c['limit']+w['limit']-widths[2,a,b];parts=dict(chord=c['reference_uncertainty_indicator'],width=w['reference_uncertainty_indicator'],mixed=mixed[a,b],quadrature=2*quad[a,b]);unc=sum(parts.values())/abs(limit);diff=abs(candidate[a,b]-limit)/abs(limit);r.update(reference_limit=pairs(limit),absolute_uncertainty_components={name:float(v) for name,v in parts.items()},reference_relative_uncertainty=float(unc),finite_candidate_to_reference_limit=float(diff),combined_relative_indicator=float(diff+unc),passed=bool(unc<.0025 and diff+unc<.01))
                else:r.update(width_reason=w.get('reason'),chord_reason=c.get('reason'))
                records.append(r)
    paths=[Path(__file__),ROOT/'rd/aero/unsteady-reference/records/ideal-pitch-frequency-assessment-policy.md',ROOT/'rd/aero/unsteady-reference/records/ptera-scalar-reference-result.json',ROOT/'rd/aero/unsteady-reference/records/ptera-core-model-diagnostic.md',ROOT/'rd/aero/unsteady-reference/records/ideal-frequency-quadrature-protocol.md',ROOT/'rd/aero/unsteady-reference/experiments/run_ptera_ideal_frequency.py',ROOT/'rd/aero/unsteady-reference/experiments/run_ptera_ideal_quadrature.py',ROOT/'rd/aero/unsteady-reference/experiments/launch_ptera_ideal_frequency.py'];USED.update(paths);result=dict(passed=all(r['passed'] for r in records),scope='FW-V01 infinitesimal heave/pitch k=.05,.2,.8; actual finite RM-core candidate versus ideal-core independent extrapolated reference, both closed wake128m',rigorous_error_bound=False,records=records,whole_uvlm_complete=False,sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(USED)});dest=ROOT/'rd/aero/unsteady-reference/records/ideal-pitch-frequency-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k!='sources'},indent=2))
if __name__=='__main__':main()
