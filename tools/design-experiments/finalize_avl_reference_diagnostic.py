"""Freeze diagnostic disposition without accepting the unresolved reference."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    names=['avl-allocation-factorial-runs.json','avl-allocation-factorial-result.json','avl-allocation-factorial-artifact-audit.json',
        'avl-span1280-observation-result.json','avl-span1280-results.md','avl-symmetric-solve-sensitivity-results.md',
        'avl-allocation-factorial-results.md','avl-reference-diagnostic-disposition.md']
    raw,obs,audit=[json.loads((F/n).read_text()) for n in names[:3]]
    assert raw['completed'] and obs['methods_criteria_pass'] and audit['all_criteria_pass']
    assert audit['runs_sha256']==digest(F/names[0]) and audit['result_sha256']==digest(F/names[1])
    assert all(obs['checks'].values()) and not obs['reference_acceptance'] and obs['m24_error_budget'] is None
    result=dict(source_sha256=digest(__file__),evidence_sha256={n:digest(F/n) for n in names},
        current_diagnostic_completed=True,disposition='not-adopted-as-accuracy-acceptance-reference',
        allowed_use=['fixed-artifact-regression','finite-grid-diagnostics'],reference_accuracy_accepted=False,
        m24_combined_error_budget=None,root_cause_fully_resolved=False,e_aero_closed=False,j1_closed=False,
        production_core_start_authorized=False,next_work='structural-model-adoption-review')
    output=F/'avl-reference-diagnostic-disposition.json'
    if output.exists():raise FileExistsError(output)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
