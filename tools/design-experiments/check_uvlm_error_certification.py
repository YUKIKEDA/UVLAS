from pathlib import Path
import json,hashlib
from uvlm_error_certification import estimate,certify
ROOT=Path(__file__).resolve().parents[2]
def main():
    rows=[];a=estimate([108,104,102,101]);assert a['resolved'] and abs(a['limit']-100)<1e-12 and a['finite_error_indicator']>1.49 and a['reference_uncertainty_indicator']<1e-10;rows.append('finite_error_distinct_from_corrected_reference_uncertainty')
    for name,values in [('oscillation',[101,99,100.5,99.75]),('false_plateau',[110,103,101,101]),('complex_rotating',[1+1j,1+.5j,1.25+.5j,1.25+.375j]),('unstable_order',[110,104,102,101.9]),('all_equal',[100]*4),('nonfinite',[100,99,98,float('nan')])]:assert not estimate(values)['resolved'];rows.append(name)
    plan=dict(case_id='wing-heave',scale=100.,relative_tolerance=.01,required=['space','time']);good=[dict(case_id='wing-heave',covers=['space'],resolved=True,absolute_indicator=.7),dict(case_id='wing-heave',covers=['time'],resolved=True,absolute_indicator=.2)];assert certify(plan,good)['qualified'];rows.append('additive_budget_pass')
    assert not certify(plan,good[:1])['qualified'];rows.append('missing_is_not_zero')
    bad=[dict(x) for x in good];bad[1]['absolute_indicator']=.4;assert 'total_budget_exceeded' in certify(plan,bad)['reasons'];rows.append('additive_budget_failure')
    bad=[dict(x) for x in good];bad[1]['case_id']='other-axis';assert 'physical_case_mismatch' in certify(plan,bad)['reasons'];rows.append('different_observation_case_rejected')
    caps=dict(plan,absolute_group_caps=dict(space=.25));assert 'individual_cap_exceeded' in certify(caps,good)['reasons'];rows.append('no_silent_budget_transfer')
    bad=[dict(good[0]),dict(good[0])];assert 'overlapping_error_groups' in certify(plan,bad)['reasons'];rows.append('overlapping_estimates_rejected')
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'tools/design-experiments/uvlm_error_certification.py',ROOT/'docs/verification/uvlm-error-certification-contract.md']};r=dict(passed=True,groups=len(rows),records=rows,sources=sources,not_a_physical_validation_certificate=True);dest=ROOT/'docs/verification/uvlm-error-certification-result.json';assert not dest.exists();dest.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
