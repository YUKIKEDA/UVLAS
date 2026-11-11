"""Predeclared quadrature-only refinement of the retained planar reduction."""
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import check_uvlm_wp1_three_dimensional as trial

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()


def main():
    points=np.array([[.2,0.,0.],[.75,0.,0.]])
    left,right=trial.planar.weights(1.,.4,points[:,0],.6);rows=[]
    for values in [[1.,0.],[0.,1.]]:
        target=values[0]*left+values[1]*right
        reference=trial.area_reference(trial.SHAPES['flat'],values,points,0.,96)
        results=[]
        for order in [24,48,96]:
            if time.monotonic()-START>60 or trial.core.memory.memory_bytes()>256*1024**2:
                raise RuntimeError('experiment_budget_exceeded')
            velocity=trial.primary(trial.SHAPES['flat'],values,points,0.,surface_order=order)
            results.append(dict(surface_order=order,velocity=velocity.tolist(),
                analytic_difference=trial.error(velocity[:,2],target),
                area_reference_difference=trial.error(velocity,reference),
                reference_eligible=trial.point_eligible(velocity,reference)))
        a=np.array(results[1]['velocity']);b=np.array(results[2]['velocity'])
        rows.append(dict(coefficients=values,cases=results,change_48_96=trial.error(a,b)))
    checks=dict(planar_analytic=all(r['cases'][-1]['analytic_difference']<=1e-9 for r in rows),
        quadrature_change=all(r['change_48_96']<=1e-11 for r in rows),
        independent_area=all(r['cases'][-1]['reference_eligible'] for r in rows))
    source=Path(__file__);parent=ROOT/'docs/verification/uvlm-wp1-three-dimensional-result.json'
    old=json.loads(parent.read_text());names=['check_uvlm_wp1_three_dimensional.py',*old['dependency_sha256']]
    result=dict(protocol='docs/verification/uvlm-wp1-three-dimensional-experiment.md#8',
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={name:hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in names},
        parent_result_sha256=hashlib.sha256(parent.read_bytes()).hexdigest(),cases=rows,
        checks=checks,all_criteria_pass=bool(all(checks.values())),
        numpy=np.__version__,scipy=trial.scipy.__version__,seconds=time.monotonic()-START,
        peak_process_working_bytes=trial.core.memory.memory_bytes(),
        limits='quadrature-only comparison at fixed planar observations; not a general near/self evaluator')
    (ROOT/'docs/verification/uvlm-wp1-planar-reduction-result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result['checks']),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':main()
