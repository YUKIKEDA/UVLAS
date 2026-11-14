"""Predeclared fixed-order near-field diagnostic, not an adaptive evaluator."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib
import json
import sys
import time
import warnings
import numpy as np
import scipy
from scipy.integrate import quad
import check_uvlm_wp1_three_dimensional as sheet

ROOT = rd_paths.ROOT
START = time.monotonic()
ROWS = []

def main():
    patch = sheet.SHAPES['flat']
    length = max(np.linalg.norm(patch[(i+1)%4]-patch[i]) for i in range(4))
    # The physical triangle diameter also includes the fixed diagonal.
    diameter = max(np.linalg.norm(patch[i]-patch[j]) for i in range(4) for j in range(4))
    for distance in [1., .5, .1, .01, .001, .0001, .00001]:
        point = np.array([[1-distance, 0., 0.]])
        actual_distance = 1-point[0, 0]
        analytic = sheet.planar.weights(1., .4, point[0, 0], .6)
        basis = []
        for j in range(2):
            mu = np.eye(2)[j]
            def integrand(t):
                r = actual_distance*np.exp(t)
                fraction = (r-actual_distance)/.4
                weight = 1-fraction if j == 0 else fraction
                return .6*weight/(2*np.pi*r*np.sqrt(r*r+.6**2))
            upper = np.log1p(.4/actual_distance)
            with warnings.catch_warnings():
                warnings.simplefilter('error')
                ref, ref_error = quad(integrand, 0., upper, epsabs=1e-10, epsrel=1e-12, limit=200)
                tight, tight_error = quad(integrand, 0., upper, epsabs=5e-11, epsrel=5e-13, limit=200)
            reference = np.array([0., 0., analytic[j]])
            target = 1e-10/length+1e-8*np.linalg.norm(reference)
            values = []
            for order in [24, 48, 96]:
                v = sheet.primary(patch, mu, point, 0., line_order=order, surface_order=order)[0]
                error = np.linalg.norm(v-reference)
                values.append(dict(order=order, velocity=v.tolist(), absolute_error=float(error),
                    relative_error=float(error/np.linalg.norm(reference)), eligible=bool(error<=target)))
            estimate = 16*np.linalg.norm(np.array(values[-1]['velocity'])-values[-2]['velocity'])
            # Conservative common diameter screen; it cannot accept an unsafe subelement.
            screen = bool(estimate<=target and diameter/actual_distance<=2)
            basis.append(dict(basis=j, analytic=float(analytic[j]), transformed_reference=ref,
                transformed_reported_error=ref_error, transformed_tight_reported_error=tight_error,
                reference_difference=float(abs(ref-analytic[j])/abs(analytic[j])),
                reference_self_difference=float(abs(ref-tight)/abs(tight)),
                target=float(target), orders=values, estimated_48_96_error=float(estimate), screen_eligible=screen))
        c = sheet.primary(patch, [1.,1.], point, 0., line_order=96, surface_order=96)[0]
        summed = np.array(basis[0]['orders'][-1]['velocity'])+basis[1]['orders'][-1]['velocity']
        ROWS.append(dict(distance=distance, finite_surface_distance=actual_distance,
            diameter_distance_ratio=float(diameter/actual_distance), basis=basis,
            constant_basis_difference=float(np.linalg.norm(c-summed)/max(np.linalg.norm(c),1/length))))
    columns = [b for r in ROWS for b in r['basis']]
    accepted = [b for b in columns if b['screen_eligible']]
    checks = dict(
        independent_references=all(b['reference_difference']<=1e-11 and b['reference_self_difference']<=1e-11 for b in columns),
        far_field_accuracy=all(b['orders'][-1]['eligible'] for b in ROWS[0]['basis']),
        near_field_fixed_order_failure_detected=any(v['relative_error']>1e-3 for r in ROWS if r['distance']<=.01 for b in r['basis'] for v in b['orders']),
        screen_no_false_accept=all(b['orders'][-1]['eligible'] for b in accepted),
        screen_accepts_far_field=all(b['screen_eligible'] for b in ROWS[0]['basis']),
        constant_basis_identity=max(r['constant_basis_difference'] for r in ROWS)<=1e-10)
    source = Path(__file__).resolve()
    loaded = sorted({Path(m.__file__).resolve() for m in list(sys.modules.values()) if getattr(m,'__file__',None)
        and rd_paths.is_experiment(m.__file__) and Path(m.__file__).resolve()!=source})
    result = dict(protocol='rd/aero/kernel/records/uvlm-near-field-screening-experiment.md',
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        dependency_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in loaded},
        patch=patch.tolist(), core=0., parent_length=float(length), screen_diameter=float(diameter), cases=ROWS,
        checks=checks, all_criteria_pass=all(checks.values()), numpy=np.__version__, scipy=scipy.__version__,
        seconds=time.monotonic()-START, peak_process_working_bytes=sheet.core.memory.memory_bytes(),
        limits='fixed-order diagnosis and conservative screen only; no adaptive or self evaluator')
    (ROOT/'rd/aero/kernel/records/uvlm-near-field-screening-result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks, accepted_columns=len(accepted),
        max_fixed_relative_error=max(v['relative_error'] for b in columns for v in b['orders']),
        seconds=result['seconds'],peak_bytes=result['peak_process_working_bytes'])),flush=True)
    if not all(checks.values()): raise SystemExit(1)

if __name__=='__main__':
    try: main()
    except Exception as exc:
        (ROOT/'build/analysis/uvlm-near-field-screening-stopped.json').write_text(
            json.dumps(dict(error=repr(exc),partial_cases=ROWS),indent=2,allow_nan=False)+'\n',encoding='utf-8')
        raise
