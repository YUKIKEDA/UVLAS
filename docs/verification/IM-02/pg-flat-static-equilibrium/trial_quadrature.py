from fractions import Fraction as Q


def ensure(condition, message):
    if not condition:
        raise ValueError(message)


def poly_add(left, right):
    terms = dict(left)
    for key, value in right.items():
        terms[key] = terms.get(key, Q()) + value
    return {key: value for key, value in terms.items() if value}


def poly_mul(left, right):
    result = {}
    for (i, j), a in left.items():
        for (k, l), b in right.items():
            key = (i+k, j+l)
            result[key] = result.get(key, Q()) + a*b
    return result


def scale(value, expression):
    return {key: value * x for key, x in expression.items()}


def energy_terms(a, length, thickness, young):
    # X=(1+H zeta/2, 2+L s, Z), u_Y=a*L*s*zeta**2 on a clamped root.
    shear = {(1, 1): 4*a*length/thickness}
    extension = {(0, 2): a}
    e_xx = scale(Q(1, 2), poly_mul(shear, shear))
    e_yy = poly_add(extension, scale(Q(1, 2), poly_mul(extension, extension)))
    e_xy = scale(Q(1, 2), poly_add(shear, poly_mul(shear, extension)))
    strain2 = poly_add(poly_mul(e_xx, e_xx), poly_mul(e_yy, e_yy))
    strain2 = poly_add(strain2, scale(2, poly_mul(e_xy, e_xy)))
    # Z width equals L; factor integrates area and full thickness with uniform averages.
    return scale(young*length**2*thickness/2, strain2)


def exact_moment(k):
    return Q() if k % 2 else Q(1, k+1)


def energy(coefficients, moment):
    return sum(c*moment(j)/(i+1) for (i, j), c in coefficients.items())


def compact_rule():
    # Solve the rational weights for symmetric dyadic nodes with exact moments 0,2,...,8.
    nodes = [Q(n, 64) for n in (10, 28, 43, 55, 62)]
    table = [[x**(2*i) for x in nodes] + [Q(1, 2*i+1)] for i in range(len(nodes))]
    for col in range(len(nodes)):
        pivot = table[col][col]
        table[col] = [x/pivot for x in table[col]]
        for row in range(len(nodes)):
            if row != col:
                multiplier = table[row][col]
                table[row] = [x-multiplier*y for x,y in zip(table[row], table[col], strict=True)]
    return [(sign*x, row[-1]/2) for x,row in zip(nodes, table, strict=True) for sign in (-1,1)]


def main():
    length, thickness, young = Q(2), Q(127, 50_000), Q(200_000_000_000)
    a = thickness/(10*length)
    coefficients = energy_terms(a, length, thickness, young)
    ensure(coefficients[(0, 8)] > 0, 'quartic Green energy must contain zeta**8')
    analytic = energy(coefficients, exact_moment)
    simpson = energy(coefficients, lambda k: sum(w*x**k for x,w in
                     [(Q(-1),Q(1,6)), (Q(),Q(2,3)), (Q(1),Q(1,6))]))
    ensure(abs(simpson-analytic)/analytic > Q(1,100), 'three thickness points cannot certify the trial')
    rule = compact_rule()
    ensure(all(w > 0 for _,w in rule), 'positive rule')
    for k in range(10):
        ensure(sum(w*x**k for x,w in rule) == exact_moment(k), 'monomial moment')
    computed = energy(coefficients, lambda k: sum(w*x**k for x,w in rule))
    ensure(computed == analytic, 'exact thickness rule must match the analytic energy')
    print(f'trial strain parameter = {a}; det(F) >= 1')
    print(f'analytic Joules = {float(analytic):.12g}')
    print(f'Simpson depth relative error = {float(abs(simpson-analytic)/analytic):.12g}')
    print('positive degree-nine thickness rule reproduces exact energy')


if __name__ == '__main__':
    main()
