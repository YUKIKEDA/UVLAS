"""Four unchanged finite-h quality cases using primary physical objective evaluation."""
from pathlib import Path
import hashlib,json,time
import check_reference_work_v3_shared_v2 as shared
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    dest=ROOT/'docs/verification/reference-work-v3-refinement-result.json';assert not dest.exists();start=time.monotonic();records=[];stop=None
    parent=json.loads((ROOT/'docs/verification/reference-work-v3-shared-v2-result.json').read_text());assert parent['all_global_checks'] and parent['all_curved_improved']
    out=ROOT/'build/analysis'/('reference-work-v3-refinement-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
    try:
        for n,h in [(4,.01),(4,.0001),(8,.001),(8,.0001)]:
            r=shared.case(.4,n,.2,h,out);records.append(r);(out/'checkpoint.json').write_text(json.dumps(dict(records=records),indent=2),encoding='utf-8');print(n,h,r['observed'],r['checks'],r['seconds'],flush=True)
            if not all(r['checks'].values()) or time.monotonic()-start>900:raise RuntimeError('refinement_necessary_gate')
    except Exception as exc:stop=repr(exc)
    complete=len(records)==4 and stop is None;ratios=[records[1]['observed']['error'][j]/records[0]['observed']['error'][j] for j in range(2)] if len(records)>=2 else None
    accuracy=complete and all(r['observed']['error'][0]<=.01 for r in records if r['n']==8) and ratios[0]<=2
    ps=[Path(__file__),Path(shared.__file__),ROOT/'tools/design-experiments/reference_work_fit_v2.py',ROOT/'docs/verification/reference-work-v3-protocol.md',ROOT/'docs/verification/reference-work-v3-refinement-preparation.json',ROOT/'docs/verification/reference-work-v3-shared-v2-result.json']
    dest.write_text(json.dumps(dict(records=records,stop_reason=stop,all_numerical_checks=complete,quartic_accuracy_gate_pass=bool(accuracy),thin_amplification=ratios,exponential_additional_observation_only=True,physical_adoption=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in ps},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes()),indent=2),encoding='utf-8')
if __name__=='__main__':main()
