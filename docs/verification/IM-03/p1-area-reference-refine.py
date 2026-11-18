"""Precision supplement to the frozen independent area reference, not product code."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import mpmath as mp


def evaluate(module, case, order, dps, step_text):
    with mp.workdps(dps):
        point,vertices,radius,beta = module.baseline(case)
        values = point+[x for row in vertices for x in row]+[radius]
        basis = module.area(values,beta,order,dps)
        jacobian = [[[mp.mpf(0) for _ in range(13)] for _ in range(3)] for _ in range(3)]
        step = mp.mpf(step_text)
        for column in range(13):
            perturbed = values.copy()
            perturbed[column] += mp.j*step
            result = module.area(perturbed,beta,order,dps)
            for k in range(3):
                for j in range(3):
                    jacobian[k][j][column] = mp.im(result[k][j])/step
        return basis,jacobian


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--generator',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('output must be new')
    spec = importlib.util.spec_from_file_location('frozen_area_generator',args.generator)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rows = []
    for case in module.cases():
        if case['id'] not in ('flat-ideal-near-1e-05','skew-core-02-face'):
            continue
        previous = evaluate(module,case,96,125,'1e-40')
        for order in (160,256,384,512):
            current = evaluate(module,case,order,125,'1e-40')
            delta = module.difference(current,previous)
            print(case['id'],order,mp.nstr(delta,8),flush=True)
            if delta < mp.mpf('1e-37'):
                break
            previous = current
        else:
            raise AssertionError('high-precision angular refinement unresolved')
        precision = evaluate(module,case,order,150,'1e-40')
        precision_delta = module.difference(precision,current)
        step = evaluate(module,case,order,150,'1e-32')
        step_delta = module.difference(precision,step)
        assert precision_delta < mp.mpf('1e-50')
        assert step_delta < mp.mpf('1e-37')
        basis,jacobian = precision
        row = dict(case)
        row.update(basis=[[mp.nstr(x,80) for x in line] for line in basis],
            jacobian=[[[mp.nstr(x,80) for x in line] for line in plane] for plane in jacobian],
            reference_order=order,reference_dps=150,complex_step='1e-40',
            angular_difference_budget_ratio=mp.nstr(delta,20),
            precision_difference_budget_ratio=mp.nstr(precision_delta,20),
            complex_step_difference_budget_ratio=mp.nstr(step_delta,20))
        rows.append(row)
    assert len(rows)==2
    result = {'id':'p1-area-reference-precision-supplement-v1',
        'scope':'empirical refinement, not a rigorous continuous-integral error bound',
        'sources':{'generator_sha256':hashlib.sha256(args.generator.read_bytes()).hexdigest(),
                   'refinement_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        'cases':rows}
    args.output.write_bytes((json.dumps(result,indent=2)+'\n').encode('utf-8'))


if __name__=='__main__':
    main()
