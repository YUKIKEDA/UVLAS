"""Independent P1 area integral: fixed barycentric fan, exact radial moments.

The boundary-vortex/sheet formula and product code are not imported. Ambient
derivatives use complex-step with the real reference fan held fixed. CLI inputs
are binary64; mpmath receives their exact binary values. Python 3.12.8,
mpmath 1.3.0 and NumPy 2.3.5 are the fixed reference environment.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path
import time

import mpmath as mp
import numpy as np


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


@lru_cache(maxsize=None)
def angular_rule(n, dps):
    with mp.workdps(dps):
        rows = []
        for index in range(1, n + 1):
            z = mp.cos(mp.pi * (mp.mpf(index) - mp.mpf('.25')) / (n + mp.mpf('.5')))
            for _ in range(30):
                previous, current = mp.mpf(1), z
                for degree in range(2, n + 1):
                    previous, current = current, ((2*degree-1)*z*current-(degree-1)*previous)/degree
                derivative = n * (z * current - previous) / (z*z - 1)
                correction = current / derivative
                z -= correction
                if abs(correction) < mp.eps * 16:
                    break
            else:
                raise AssertionError('angular root did not converge')
            # Refresh the derivative at the final root, not its Newton predecessor.
            previous, current = mp.mpf(1), z
            for degree in range(2, n + 1):
                previous, current = current, ((2*degree-1)*z*current-(degree-1)*previous)/degree
            derivative = n * (z * current - previous) / (z*z - 1)
            rows.append(((z+1)/2, 1/((1-z*z)*derivative*derivative)))
        rows.sort()
        assert abs(sum(w for _, w in rows)-1) < mp.mpf('1e-75')
        return rows


def primitives(u, t2):
    radius = mp.sqrt(u*u+t2)
    cube = radius**3
    return [u/(t2*t2*radius)-u**3/(3*t2*t2*cube),
            -1/(3*cube), u**3/(3*t2*cube),
            -1/radius+t2/(3*cube),
            mp.asinh(u/mp.sqrt(t2))-u/radius-u**3/(3*cube)]


def radial_moments(a, b, c):
    shift = b/(2*a)
    t2 = c/a-shift*shift
    lower, upper = primitives(shift, t2), primitives(1+shift, t2)
    return [sum(comb(m,j)*(-shift)**(m-j)*(upper[j]-lower[j])
                for j in range(m+1))/a**mp.mpf('2.5') for m in range(5)]


def baseline(case):
    vertices = [[mp.mpf(value) for value in row] for row in case['vertices']]
    query = case['query']
    if query['kind'] == 'free':
        point = [mp.mpf(value) for value in query['point']]
    elif query['kind'] == 'face':
        u, v = mp.mpf(query['u']), mp.mpf(query['v'])
        point = [(1-u-v)*vertices[0][j]+u*vertices[1][j]+v*vertices[2][j] for j in range(3)]
    elif query['kind'] == 'edge':
        first, second = query['edge_index'], (query['edge_index']+1)%3
        t = mp.mpf(query['t'])
        point = [(1-t)*vertices[first][j]+t*vertices[second][j] for j in range(3)]
    else:
        point = vertices[query['index']].copy()
    e, f = sub(vertices[1], vertices[0]), sub(vertices[2], vertices[0])
    r = sub(point, vertices[0])
    ee, ef, ff = dot(e,e), dot(e,f), dot(f,f)
    determinant = ee*ff-ef*ef
    u, v = (ff*dot(r,e)-ef*dot(r,f))/determinant, (ee*dot(r,f)-ef*dot(r,e))/determinant
    beta = [1-u-v,u,v]
    # An interior anchor avoids signed, individually singular fans for ideal
    # coplanar EXTERIOR points. The radial moments handle general quadratic D.
    if any(value < 0 or value > 1 for value in beta):
        beta = [mp.mpf(1)/3]*3
    return point, vertices, mp.mpf(case['radius']), beta


def area(values, beta, order, dps):
    point = values[:3]
    vertices = [values[3+3*i:6+3*i] for i in range(3)]
    epsilon = values[12]
    e, f = sub(vertices[1],vertices[0]), sub(vertices[2],vertices[0])
    cr = cross(e,f)
    jacobian = mp.sqrt(dot(cr,cr))
    normal = [x/jacobian for x in cr]
    centre = [sum(beta[k]*vertices[k][j] for k in range(3)) for j in range(3)]
    a = sub(point,centre)
    an = dot(a,normal)
    aa = dot(a,a)
    result = [[mp.mpf(0) for _ in range(3)] for _ in range(3)]
    for first, second, opposite in ((0,1,2),(1,2,0),(2,0,1)):
        gamma = beta[opposite]
        if gamma == 0:
            continue
        for q, weight in angular_rule(order,dps):
            angular = [mp.mpf(0)]*3
            angular[first], angular[second] = 1-q,q
            direction = [sum(angular[k]*vertices[k][j] for k in range(3))-centre[j] for j in range(3)]
            ab, bb = dot(a,direction), dot(direction,direction)
            moment = radial_moments(bb,-2*ab,aa+epsilon*epsilon)
            # direction lies in the source plane: direction.normal is identically zero,
            # including its geometry derivative. These are AREA-integrand coefficients.
            c0 = [(aa-2*epsilon*epsilon)*normal[j]-3*a[j]*an for j in range(3)]
            c1 = [-2*ab*normal[j]+3*direction[j]*an for j in range(3)]
            c2 = [bb*normal[j] for j in range(3)]
            factor = jacobian*gamma*weight/(4*mp.pi)
            for k in range(3):
                slope = angular[k]-beta[k]
                for j in range(3):
                    result[k][j] += factor*(beta[k]*c0[j]*moment[1]
                        +(beta[k]*c1[j]+slope*c0[j])*moment[2]
                        +(beta[k]*c2[j]+slope*c1[j])*moment[3]+slope*c2[j]*moment[4])
    return result


def evaluate(case, order, dps):
    with mp.workdps(dps):
        point, vertices, epsilon, beta = baseline(case)
        values = point+[x for row in vertices for x in row]+[epsilon]
        basis = area(values,beta,order,dps)
        derivatives = [[[mp.mpf(0) for _ in range(13)] for _ in range(3)] for _ in range(3)]
        step = mp.mpf('1e-32')
        for column in range(13):
            perturbed = values.copy()
            perturbed[column] += mp.j*step
            actual = area(perturbed,beta,order,dps)
            for k in range(3):
                for j in range(3):
                    derivatives[k][j][column] = mp.im(actual[k][j])/step
        return basis, derivatives


def flatten(result):
    b,j = result
    return [x for row in b for x in row]+[x for plane in j for row in plane for x in row]


def difference(a,b):
    with mp.workdps(110):
        return max(abs(x-y)/max(mp.mpf('1e-10'),mp.mpf('1e-8')*abs(x))
                   for x,y in zip(flatten(a),flatten(b)))


def cases():
    flat = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]]
    skew = [[.25,-.125,.5],[1.25,.25,.625],[.375,1.125,.875]]
    result = []
    for distance in (.01,.001,.00001):
        result.append({'id':f'flat-ideal-near-{distance:g}','vertices':flat,'radius':0.,
                       'query':{'kind':'free','point':[.25,.25,distance]}})
    for radius in (.02,.08):
        result.append({'id':f'flat-core-{radius:g}-face','vertices':flat,'radius':radius,
                       'query':{'kind':'face','u':.25,'v':.25}})
    for kind in ('edge','vertex'):
        query = {'kind':'edge','edge_index':0,'t':.5} if kind=='edge' else {'kind':'vertex','index':0}
        result.append({'id':f'flat-core-08-{kind}','vertices':flat,'radius':.08,'query':query})
    result += [
        {'id':'skew-core-08-free','vertices':skew,'radius':.08,'query':{'kind':'free','point':[.5,.3,1.]}},
        {'id':'skew-core-02-face','vertices':skew,'radius':.02,'query':{'kind':'face','u':.25,'v':.5}},
        {'id':'flat-far-core-08','vertices':flat,'radius':.08,'query':{'kind':'free','point':[5.,4.,6.]}},
        {'id':'flat-ideal-coplanar-exterior','vertices':flat,'radius':0.,'query':{'kind':'free','point':[3.,3.,0.]}},
        {'id':'flat-ideal-line-extension','vertices':flat,'radius':0.,'query':{'kind':'free','point':[3.,0.,0.]}},
        {'id':'skew-ideal-free','vertices':skew,'radius':0.,'query':{'kind':'free','point':[.5,.3,1.]}},
    ]
    quad = [[0.,0.,0.],[1.,.125,.125],[1.,1.,.375],[0.,1.,0.]]
    for index, ids in enumerate(((0,2,1),(0,3,2))):
        for radius in (0.,.02,.08):
            result.append({'id':f'nonplanar-quad-child-{index}-core-{radius:g}','vertices':[quad[i] for i in ids],
                           'radius':radius,'query':{'kind':'free','point':[.375,.375,.875]}})
    return result


def verify_primitives():
    with mp.workdps(90):
        for t2 in (mp.mpf('1e-10'),mp.mpf('.125')):
            for u in (mp.mpf('-.75'),mp.mpf(0),mp.mpf('.25'),mp.mpf('1.5')):
                for degree in range(5):
                    actual = mp.diff(lambda x:primitives(x,t2)[degree],u)
                    want = u**degree/(u*u+t2)**mp.mpf('2.5')
                    assert abs(actual-want) <= mp.mpf('1e-60')*max(1,abs(want))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--only',default='')
    parser.add_argument('--legacy-source',type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('output must be new')
    verify_primitives()
    legacy = None
    if args.legacy_source:
        spec = importlib.util.spec_from_file_location('isolated_area_reference',args.legacy_source)
        assert spec is not None and spec.loader is not None
        legacy = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(legacy)
    selected = [case for case in cases() if not args.only or case['id']==args.only]
    assert selected
    rows = []
    for case in selected:
        start = time.perf_counter()
        previous = evaluate(case,24,90)
        for order in (40,64,96,128,192):
            current = evaluate(case,order,90)
            delta = difference(current,previous)
            print(case['id'],order,mp.nstr(delta,8),flush=True)
            if delta <= mp.mpf('.005'):
                break
            previous = current
        else:
            raise AssertionError(f"unresolved angular convergence: {case['id']}")
        precise = evaluate(case,order,110)
        precision_delta = difference(precise,current)
        assert precision_delta <= mp.mpf('1e-20')
        if legacy is not None and case['id'] in ('flat-core-0.08-face','skew-core-08-free'):
            with mp.workdps(90):
                point,_,_,_ = baseline(case)
                point = np.array([float(value) for value in point])
                original_b = legacy.integrate(point,np.array(case['vertices']),case['radius'],192,False)
                original_j = legacy.jacobian(point,np.array(case['vertices']),case['radius'],192)
                original = ([[mp.mpf(float(x)) for x in row] for row in original_b],
                            [[[mp.mpf(float(x)) for x in row] for row in plane] for plane in original_j])
                assert difference(precise,original) <= mp.mpf('.05')
        basis,jacobian = precise
        row = dict(case)
        row.update(basis=[[mp.nstr(x,65) for x in r] for r in basis],
            jacobian=[[[mp.nstr(x,65) for x in r] for r in plane] for plane in jacobian],
            reference_order=order, reference_dps=110,
            angular_difference_budget_ratio=mp.nstr(delta,16),
            precision_difference_budget_ratio=mp.nstr(precision_delta,16))
        rows.append(row)
        print('completed',case['id'],round(time.perf_counter()-start,2),'seconds',flush=True)
    output = {'id':'p1-area-radial-reference-v1','environment':{'python':'3.12.8','mpmath':mp.__version__,'numpy':np.__version__},
        'method':'independent area integral; fixed reference barycentric fan; exact radial moments; ambient complex-step',
        'complex_step':'1e-32','input_values':'exact binary64; material point evaluated from its original formula',
        'source':{'id':'p1-area-radial-reference-producer-v1','sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        'primitive_derivative_checks':40,'cases':rows}
    if args.legacy_source:
        output['legacy_control'] = {'id':'uvlm-induction-reference-source',
            'sha256':hashlib.sha256(args.legacy_source.read_bytes()).hexdigest(),
            'scope':'two nonsingular control cases; original source imported only from isolated copy'}
    args.output.write_bytes((json.dumps(output,indent=2)+'\n').encode('utf-8'))


if __name__ == '__main__':
    main()
