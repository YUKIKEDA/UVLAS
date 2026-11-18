"""Independent rational Sturm certificates for Gauss-Legendre 8/12 on [0,1].

No product or R&D code is imported. All root counts, weight bounds and moment
checks use exact integer/Fraction arithmetic. The output is a test fixture.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import gcd, lcm
from pathlib import Path


def trim(p: list[Q]) -> list[Q]:
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def legendre(n: int) -> list[Q]:
    previous, current = [Q(1)], [Q(0), Q(1)]
    for degree in range(2, n + 1):
        result = [Q(0)] * (degree + 1)
        for index, value in enumerate(current):
            result[index + 1] += Q(2 * degree - 1, degree) * value
        for index, value in enumerate(previous):
            result[index] -= Q(degree - 1, degree) * value
        previous, current = current, result
    return current


def derivative(p: list[Q]) -> list[Q]:
    return [index * value for index, value in enumerate(p)][1:]


def remainder(p: list[Q], divisor: list[Q]) -> list[Q]:
    result = p.copy()
    while len(result) >= len(divisor) and any(result):
        degree = len(result) - len(divisor)
        factor = result[-1] / divisor[-1]
        for index, value in enumerate(divisor):
            result[index + degree] -= factor * value
        trim(result)
    return result


def integer_polynomial(p: list[Q]) -> list[int]:
    denominator = lcm(*(value.denominator for value in p))
    integers = [value.numerator * (denominator // value.denominator) for value in p]
    common = gcd(*integers)
    return [value // common for value in integers]


def sturm(p: list[Q]) -> list[list[int]]:
    sequence = [p, derivative(p)]
    while len(sequence[-1]) > 1:
        rest = remainder(sequence[-2], sequence[-1])
        assert any(rest), "Legendre polynomial must be square-free"
        sequence.append([-value for value in rest])
    return [integer_polynomial(row) for row in sequence]


def sign_at(p: list[int], x: Q) -> int:
    value, power = p[-1], 1
    for coefficient in reversed(p[:-1]):
        power *= x.denominator
        value = value * x.numerator + coefficient * power
    return (value > 0) - (value < 0)


def variations(sequence: list[list[int]], x: Q) -> int:
    signs = [sign_at(row, x) for row in sequence]
    nonzero = [sign for sign in signs if sign]
    return sum(a != b for a, b in zip(nonzero, nonzero[1:]))


def isolate(sequence: list[list[int]], n: int) -> list[tuple[Q, Q]]:
    pending = [(Q(-1), Q(1))]
    result = []
    while pending:
        low, high = pending.pop()
        count = variations(sequence, low) - variations(sequence, high)
        if count == 0:
            continue
        if count > 1:
            middle = (low + high) / 2
            assert sign_at(sequence[0], middle) != 0
            pending.extend([(middle, high), (low, middle)])
            continue
        while high - low > Q(1, 1 << 310):
            middle = (low + high) / 2
            assert sign_at(sequence[0], middle) != 0
            if variations(sequence, low) - variations(sequence, middle) == 1:
                high = middle
            else:
                low = middle
        result.append((low, high))
    result.sort()
    assert len(result) == n
    return result


Interval = tuple[Q, Q]


def add(a: Interval, b: Interval) -> Interval:
    return a[0] + b[0], a[1] + b[1]


def mul(a: Interval, b: Interval) -> Interval:
    products = [x * y for x in a for y in b]
    return min(products), max(products)


def square(a: Interval) -> Interval:
    low = Q(0) if a[0] <= 0 <= a[1] else min(a[0] ** 2, a[1] ** 2)
    return low, max(a[0] ** 2, a[1] ** 2)


def polynomial_interval(p: list[Q], x: Interval) -> Interval:
    result = (p[-1], p[-1])
    for coefficient in reversed(p[:-1]):
        result = add(mul(result, x), (coefficient, coefficient))
    return result


def quantize(interval: Interval) -> tuple[int, int]:
    low, high = (value * (1 << 256) for value in interval)
    return low.numerator // low.denominator, -((-high.numerator) // high.denominator)


def decode(interval: tuple[int, int]) -> Interval:
    return Q(interval[0], 1 << 256), Q(interval[1], 1 << 256)


def moments(nodes: list[Interval], weights: list[Interval], n: int) -> list[dict[str, object]]:
    powers = [(Q(1), Q(1)) for _ in nodes]
    result = []
    for degree in range(2 * n):
        total = (Q(0), Q(0))
        for index, weight in enumerate(weights):
            total = add(total, mul(weight, powers[index]))
            powers[index] = mul(powers[index], nodes[index])
        exact = Q(1, degree + 1)
        assert total[0] <= exact <= total[1], (n, degree)
        assert total[1] - total[0] <= Q(1, 1 << 240)
        result.append({"degree": degree, "exact": str(exact), "enclosed": True})
    return result


def table(n: int) -> dict[str, object]:
    p = legendre(n)
    sequence = sturm(p)
    roots = isolate(sequence, n)
    rows, node_intervals, weight_intervals = [], [], []
    for root in roots:
        derivative_range = polynomial_interval(derivative(p), root)
        one_minus_square = (1 - square(root)[1], 1 - square(root)[0])
        denominator = mul(one_minus_square, square(derivative_range))
        assert denominator[0] > 0
        # Mapped weights on [0,1]: w=1/((1-z*z)*P'_n(z)^2).
        weight = (1 / denominator[1], 1 / denominator[0])
        node = ((root[0] + 1) / 2, (root[1] + 1) / 2)
        node_pair, weight_pair = quantize(node), quantize(weight)
        mapped_back = (2 * decode(node_pair)[0] - 1, 2 * decode(node_pair)[1] - 1)
        counts = (variations(sequence, mapped_back[0]), variations(sequence, mapped_back[1]))
        assert counts[0] - counts[1] == 1
        node_intervals.append(decode(node_pair))
        weight_intervals.append(decode(weight_pair))
        rows.append({"node": {"lo": str(node_pair[0]), "hi": str(node_pair[1])},
                     "weight": {"lo": str(weight_pair[0]), "hi": str(weight_pair[1])},
                     "sturm_variations": list(counts)})
    assert all(a[1] < b[0] for a, b in zip(node_intervals, node_intervals[1:]))
    checks = moments(node_intervals, weight_intervals, n)
    damaged_weights = weight_intervals.copy()
    damaged_weights[0] = (Q(0), Q(0))
    assert sum(w[1] for w in damaged_weights) < 1
    return {"order": n, "exponent": -256, "entries": rows,
            "polynomial_integer_coefficients": [str(x) for x in integer_polynomial(p)],
            "moment_certificates": checks, "missing_weight_mutation_rejected": True}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = {"id": "p1-gauss-legendre-reference-v1", "interval": "[0,1]",
              "method": "exact rational Legendre recurrence and Sturm root counts",
              "rounding": "outward to integers times 2^-256",
              "producer": {"id": "p1-gauss-reference-producer-v1",
                           "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
              "tables": [table(8), table(12)]}
    args.output.write_bytes((json.dumps(output, indent=2) + "\n").encode("utf-8"))
    print("2 tables / 20 roots / 40 exact moment enclosures / 2 mutation witnesses passed")


if __name__ == "__main__":
    main()
