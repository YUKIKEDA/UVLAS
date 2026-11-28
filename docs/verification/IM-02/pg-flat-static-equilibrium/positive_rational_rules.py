from fractions import Fraction as Q
from math import cos, pi


def multiply_linear(coeffs, x):
    product = [Q()]*(len(coeffs)+1)
    for j, c in enumerate(coeffs):
        product[j] -= c*x
        product[j+1] += c
    return product


def integrate(coeffs):
    return sum(c/(j+1) for j,c in enumerate(coeffs))


def cardinal_weight(nodes, chosen):
    numerator = [Q(1)]
    denominator = Q(1)
    for k, x in enumerate(nodes):
        if k != chosen:
            numerator = multiply_linear(numerator, x)
            denominator *= nodes[chosen]-x
    return integrate(numerator)/denominator


def rule(degree, fractional_bits=24):
    count = degree+1
    scale = 1 << fractional_bits
    # Round Chebyshev roots strictly inside (0,1). Certify exactness and positivity
    # of the rational rule independently of the floating-point node generation.
    nodes = sorted(Q(round(scale*(1+cos((2*k+1)*pi/(2*count)))/2), scale)
                   for k in range(count))
    weights = [cardinal_weight(nodes, k) for k in range(count)]
    if not all(0 < x < 1 and w > 0 for x,w in zip(nodes, weights, strict=True)):
        raise ArithmeticError('rule fails positivity')
    if any(sum(w*x**k for x,w in zip(nodes, weights, strict=True)) != Q(1,k+1)
           for k in range(degree+1)):
        raise ArithmeticError('rule fails a polynomial moment')
    return nodes, weights


def main():
    for degree in (8, 9, 19, 36, 37):
        nodes,weights=rule(degree)
        max_bits=max(max(w.numerator.bit_length(), w.denominator.bit_length()) for w in weights)
        print(f'degree {degree:2d}, positive nodes {len(nodes):2d}, '
              f'maximum weight integer width {max_bits:3d} bits, '
              f'minimum weight {float(min(weights)):.10g}')


if __name__ == '__main__':
    main()
