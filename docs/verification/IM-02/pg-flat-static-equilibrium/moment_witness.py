"""Exact start witness for flat, constant-thickness, isotropic StVK."""

from fractions import Fraction as Q
from itertools import product
from math import factorial
import json

H, L, W = Q(127, 50_000), Q(2), Q(11, 5)
YOUNG, NU = Q(13_000_000_000), Q(2, 7)
MU = YOUNG / (2 * (1 + NU))
LAM = YOUNG * NU / ((1 + NU) * (1 - 2 * NU))
Z = (0, 0, 0)
ONE = {Z: Q(1)}
ZERO = {}
AXES = range(3)
BITS = {"coefficients": 0}


def require(ok, msg):
    if not ok:
        raise AssertionError(msg)


def add(*terms):
    out = {}
    for polynomial in terms:
        for powers, value in polynomial.items():
            out[powers] = out.get(powers, Q()) + value
    return {p: v for p, v in out.items() if v}


def scale(a, poly):
    return {p: a * v for p, v in poly.items()} if a else {}


def multiply(a, b):
    out = {}
    for p, x in a.items():
        for q, y in b.items():
            powers = tuple(i + j for i, j in zip(p, q, strict=True))
            out[powers] = out.get(powers, Q()) + x * y
    return {p: v for p, v in out.items() if v}


def partial(poly, axis):
    result = {}
    for powers, coeff in poly.items():
        if powers[axis]:
            index = list(powers)
            index[axis] -= 1
            result[tuple(index)] = coeff * powers[axis]
    return result


def substitute(poly, args):
    out = ZERO
    for powers, coefficient in poly.items():
        term = {Z: coefficient}
        for axis, exponent in enumerate(powers):
            for _ in range(exponent):
                term = multiply(term, args[axis])
        out = add(out, term)
    return out


BASIS = [{(1, 0, 0): Q(1)}, {(1, 1, 0): Q(1)},
         {(1, 0, 1): Q(1)}, {(1, 0, 2): Q(1)}]
# X=(1+H ζ/2,2+L s,3+W t) in world Cartesian XYZ.
GRADS = tuple(tuple(scale(alpha, partial(poly, axis))
                    for alpha, axis in ((2/H, 2), (1/L, 0), (1/W, 1))) for poly in BASIS)
RANGES = (Q(0), Q(1, 3), Q(1))
LAYER_CUTS = (Q(-1), Q(-1, 3), Q(1))
TRIANGLES = tuple(
    tri for lo, hi in zip(RANGES[:-1], RANGES[1:], strict=True)
    for tri in (((Q(), lo), (Q(1), lo), (Q(1), hi)),
                ((Q(), lo), (Q(1), hi), (Q(), hi)))
)


def affine(a, b, c):
    return {k: v for k, v in {Z: a, (1, 0, 0): b-a, (0, 1, 0): c-a}.items() if v}


def moment(powers, physical_scale):
    i, j, k = powers
    return physical_scale * Q(factorial(i)*factorial(j), factorial(i+j+2)*(k+1))


def exact_integral(poly, physical_scale):
    for value in poly.values():
        BITS["coefficients"] = max(BITS["coefficients"],
                                   value.numerator.bit_length(), value.denominator.bit_length())
    return sum((c * moment(powers, physical_scale) for powers, c in poly.items()), Q())


def material(strain):
    trace = add(*(strain[i][i] for i in AXES))
    return [[add(scale(2*MU, strain[i][j]), scale(LAM if i == j else 0, trace))
             for j in AXES] for i in AXES]


def matrix(a, b):
    return [[add(*(multiply(a[i][k], b[k][j]) for k in AXES))
             for j in AXES] for i in AXES]


def transpose(a):
    return list(zip(*a, strict=True))


def gradient(values, grads, state):
    return [[add(ONE if state and i == j else ZERO,
                 *(scale(q[i], g[j]) for q, g in zip(values, grads, strict=True)))
             for j in AXES] for i in AXES]


def geometry(tri, layer):
    a, b, c = tri
    lo, hi = layer
    da, db = (b[0]-a[0], b[1]-a[1]), (c[0]-a[0], c[1]-a[1])
    ds = da[0]*db[1]-db[0]*da[1]
    require(ds > 0 and hi > lo, 'reference geometry')
    args = (affine(a[0], b[0], c[0]), affine(a[1], b[1], c[1]),
            {Z: lo, (0, 0, 1): hi-lo})
    return tuple(tuple(substitute(g, args) for g in row) for row in GRADS), ds*L*W*H*(hi-lo)/2


def flatten(rows):
    return tuple(value for row in rows for value in row)


def actions(values, direction, cuts=LAYER_CUTS, geometric_tangent=True, energy_only=False):
    energy = Q()
    force = [[Q() for _ in AXES] for _ in BASIS]
    tangent = [[Q() for _ in AXES] for _ in BASIS]
    max_depth_degree = 0
    max_infty = Q()
    for tri, layer in product(TRIANGLES, zip(cuts[:-1], cuts[1:], strict=True)):
        grads, dvol = geometry(tri, layer)
        f = gradient(values, grads, True)
        perturb = max(sum(sum(map(abs, add(f[i][j], scale(-1, ONE) if i == j else ZERO).values()))
                          for j in AXES) for i in AXES)
        require(perturb < 1, 'uniform positive determinant guard for named witness')
        max_infty = max(perturb, max_infty)
        c = matrix(transpose(f), f)
        e = [[scale(Q(1, 2), add(c[i][j], scale(-1, ONE) if i == j else ZERO))
              for j in AXES] for i in AXES]
        s = material(e)
        local_energy = scale(Q(1, 2), add(*(multiply(s[i][j], e[i][j]) for i, j in product(AXES, repeat=2))))
        max_depth_degree = max(max_depth_degree, max((powers[2] for powers in local_energy), default=0))
        energy += exact_integral(local_energy, dvol)
        if energy_only:
            continue
        p = matrix(f, s)
        h = gradient(direction, grads, False)
        q = matrix(transpose(f), h)
        de = [[scale(Q(1, 2), add(q[i][j], q[j][i])) for j in AXES] for i in AXES]
        dp = matrix(f, material(de))
        if geometric_tangent:
            gk = matrix(h, s)
            dp = [[add(dp[i][j], gk[i][j]) for j in AXES] for i in AXES]
        # Integrate stress against each distinct gradient monomial only once.
        support = [set().union(*(gi[a].keys() for gi in grads)) for a in AXES]
        for targets, tensor in ((force, p), (tangent, dp)):
            for i, j in product(AXES, repeat=2):
                moments = {beta: sum((v*moment(tuple(a+b for a,b in zip(gamma,beta,strict=True)), dvol)
                                      for gamma, v in tensor[i][j].items()), Q())
                           for beta in support[j]}
                for k, grad in enumerate(grads):
                    targets[k][i] += sum(v * moments[powers] for powers, v in grad[j].items())
    return energy, flatten(force), flatten(tangent), max_depth_degree, max_infty


def dot(a, b):
    return sum(x*y for x,y in zip(a,b,strict=True))


def numeric_gradient_cartesian(a, g, alpha):
    return tuple(tuple(Q(i == j and alpha) + sum(q[i]*v[j] for q,v in zip(a,g,strict=True))
                       for j in AXES) for i in AXES)


def independent(values, direction):
    # Analytic Cartesian derivative at fixed positions. This does not consume any
    # triangle, affine, strain, stress, or moment object from the candidate.
    # Boole's rule integrates through degree five on either Cartesian planar axis.
    plane = tuple((Q(k, 4), Q(w, 90)) for k,w in enumerate((7,32,12,32,7)))
    # Independent ten-point symmetric positive rational depth rule from the
    # existing start witness, again checked against the full Cartesian interval.
    from trial_quadrature import compact_rule
    thickness = compact_rule()
    require(all(w > 0 for _,w in thickness), 'positive Cartesian depth weights')
    for degree in range(9):
        require(sum(w*x**degree for x,w in thickness) ==
                (Q() if degree % 2 else Q(1,degree+1)), 'Cartesian depth moment')
    u = Q()
    f, kv = [Q() for _ in range(12)], [Q() for _ in range(12)]
    for (s, ws), (t, wt), (zeta, wz) in product(plane,plane,thickness):
        grads = ((Q(), 1/L, Q()),
                 (Q(), t/L, s/W),
                 (2*s/H, zeta/L, Q()),
                 (4*s*zeta/H, zeta*zeta/L, Q()))
        F = numeric_gradient_cartesian(values, grads, True)
        D = numeric_gradient_cartesian(direction, grads, False)
        ft, dt = tuple(zip(*F, strict=True)), tuple(zip(*D, strict=True))
        # B=F^T F-I, C=d_eta B. The scalar energy uses normalized invariants;
        # first/mixed coefficients in eps, eta give force and K·direction.
        B = [[dot(ft[i],ft[j])-Q(i==j) for j in AXES] for i in AXES]
        C = [[dot(ft[i],dt[j])+dot(dt[i],ft[j]) for j in AXES] for i in AXES]
        trB, trC = sum(B[i][i] for i in AXES), sum(C[i][i] for i in AXES)
        vol = H*L*W*ws*wt*wz
        u += vol*(MU/4*sum(x*x for x in flatten(B))+LAM/8*trB*trB)
        for k, g in enumerate(grads):
            for a in AXES:
                A = [[F[a][i]*g[j]+g[i]*F[a][j] for j in AXES] for i in AXES]
                M = [[D[a][i]*g[j]+g[i]*D[a][j] for j in AXES] for i in AXES]
                trA, trM = sum(A[i][i] for i in AXES), sum(M[i][i] for i in AXES)
                value = MU/2*dot(flatten(A),flatten(B))+LAM/4*trB*trA
                mixed = (MU/2*(dot(flatten(A),flatten(C))+dot(flatten(B),flatten(M)))
                         + LAM/4*(trA*trC+trB*trM))
                f[3*k+a] += vol*value
                kv[3*k+a] += vol*mixed
    return u, tuple(f), tuple(kv)


def vectors(seed, divisor):
    return tuple(tuple(Q(((seed+2*i+3*j)*(seed+i+4*j)) % 13 - 6)*H/divisor
                       for j in AXES) for i in range(len(BASIS)))


def move(q, d, alpha):
    return tuple(tuple(x+alpha*y for x,y in zip(a,b,strict=True)) for a,b in zip(q,d,strict=True))


def main():
    entries = []
    for index in (1, 2, 4):
        q, d, v = vectors(index,256), vectors(index+1,90), vectors(index+3,80)
        u, f, kv, degree, inf = actions(q, d)
        check = independent(q,d)
        require((u,f,kv) == check, f'Cartesian energy, all 12 force and tangent output components: {index}')
        other = actions(q,v)
        require(dot(flatten(v),kv) == dot(flatten(d),other[2]), f'bilinear symmetry: {index}')
        require(actions(q,d,cuts=(Q(-1),Q(1)))[:3] == (u,f,kv), f'layer partition: {index}')
        step = Q(1, 8)
        energies = [actions(move(q,d,step*a),d,energy_only=True)[0] for a in (-2,-1,0,1,2)]
        um2, um1, u0, up1, up2 = energies
        require((um2-8*um1+8*up1-up2)/(12*step) == dot(f,flatten(d)), 'exact quartic force stencil')
        require((-up2+16*up1-30*u0+16*um1-um2)/(12*step*step) == dot(kv,flatten(d)), 'exact quartic tangent stencil')
        mutant = actions(q,d,geometric_tangent=False)
        count = sum(a!=b for a,b in zip(kv,mutant[2],strict=True))
        require(count > 0, 'omit HS mutation survives')
        require(degree == 8, 'nonlinear trial does not exercise depth degree eight')
        entries.append(dict(seed=index, energy_joules=str(u), det_guard_infty=str(inf),
                            depth_degree=degree, mutant_tangent_components_detected=count))
    print(json.dumps(dict(fixture='4 shape functions, 12 Cartesian columns, 4 triangles and 2 layers',
                          observations=entries, max_retained_coefficient_integer_bits=BITS['coefficients']),
                     indent=2))


if __name__ == '__main__':
    main()
