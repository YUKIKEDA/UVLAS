"""Independent Fraction oracle: global antiderivatives, then binomial pullback.

No product or R&D code is imported. Each cell denotes its open interior;
the separate one-sided cell convention supplies its boundary trace.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from functools import cache
from math import comb
from pathlib import Path

Poly = dict[tuple[int, int], F]


def clean(p: Poly) -> Poly:
    return {ij: c for ij, c in p.items() if c}


def add(*polys: Poly) -> Poly:
    result: Poly = {}
    for p in polys:
        for ij, c in p.items():
            result[ij] = result.get(ij, F(0)) + c
    return clean(result)


def scale(p: Poly, c: F | int) -> Poly:
    return clean({ij: v * c for ij, v in p.items()})


def mul(a: Poly, b: Poly) -> Poly:
    result: Poly = {}
    for (i, j), x in a.items():
        for (k, ell), y in b.items():
            ij = i + k, j + ell
            result[ij] = result.get(ij, F(0)) + x * y
    return clean(result)


def monomial(i: int, j: int, c: F | int = 1) -> Poly:
    return clean({(i, j): F(c)})


ONE = monomial(0, 0)
S = monomial(1, 0)
T = monomial(0, 1)


def derivative(p: Poly, axis: int) -> Poly:
    result: Poly = {}
    for ij, c in p.items():
        if ij[axis]:
            target = list(ij)
            target[axis] -= 1
            result[tuple(target)] = c * ij[axis]
    return result


def trace(p: Poly, axis: int, x: F) -> Poly:
    result: Poly = {}
    for ij, c in p.items():
        target = list(ij)
        target[axis] = 0
        key = tuple(target)
        result[key] = result.get(key, F(0)) + c * x ** ij[axis]
    return clean(result)


def integrate_from(p: Poly, axis: int, lower: F) -> Poly:
    primitive: Poly = {}
    for ij, c in p.items():
        target = list(ij)
        target[axis] += 1
        primitive[tuple(target)] = c / target[axis]
    return add(primitive, scale(trace(primitive, axis, lower), -1))


def pullback(p: Poly, n: int, cs: int, ct: int) -> Poly:
    h = F(1, 4 * n)
    x, y = F(2 * cs + 1, 4 * n), F(2 * ct + 1, 4 * n)
    result: Poly = {}
    for (i, j), c in p.items():
        for a in range(i + 1):
            for b in range(j + 1):
                value = c * comb(i, a) * comb(j, b) * x ** (i - a) * y ** (j - b) * h ** (a + b)
                result[a, b] = result.get((a, b), F(0)) + value
    return clean(result)


@cache
def legendre(k: int) -> Poly:
    if k == 0:
        return ONE
    if k == 1:
        return S
    return scale(
        add(scale(mul(S, legendre(k - 1)), 2 * k - 1), scale(legendre(k - 2), 1 - k)), F(1, k)
    )


def affine_legendre(k: int, axis: int, a: int, b: int) -> Poly:
    result: Poly = {}
    for (i, _), c in legendre(k).items():
        for j in range(i + 1):
            result = add(
                result,
                monomial(
                    j if axis == 0 else 0,
                    j if axis == 1 else 0,
                    c * comb(i, j) * a**j * b ** (i - j),
                ),
            )
    return result


@cache
def basis(n: int, degree: int, repeat: int, index: int, *, axis: int, cell: int) -> Poly:
    knots = (
        (F(0),) * (degree + 1)
        + tuple(F(j, n) for j in range(1, n) for _ in range(repeat))
        + (F(1),) * (degree + 1)
    )
    midpoint = F(2 * cell + 1, 4 * n)
    coordinate = S if axis == 0 else T

    @cache
    def recurse(i: int, d: int) -> Poly:
        if d == 0:
            return ONE if knots[i] <= midpoint < knots[i + 1] else {}
        result: Poly = {}
        left, right = knots[i + d] - knots[i], knots[i + d + 1] - knots[i + 1]
        if left:
            result = add(
                result,
                scale(mul(add(coordinate, scale(ONE, -knots[i])), recurse(i, d - 1)), 1 / left),
            )
        if right:
            result = add(
                result,
                scale(
                    mul(
                        add(scale(ONE, knots[i + d + 1]), scale(coordinate, -1)),
                        recurse(i + 1, d - 1),
                    ),
                    1 / right,
                ),
            )
        return result

    return recurse(index, degree)


class Oracle:
    def __init__(self, n: int, x0: F, width: F, skew: F, height: Poly):
        self.n, self.x0, self.width, self.skew, self.height = n, x0, width, skew, height
        self.r = (n + 2) ** 2 + 1
        self.d = (2 * n + 1) ** 2 + 3
        self.p = 96 * n
        self.c = (n + 1) ** 2
        self.total = self.r + self.d + self.p + self.c
        self.fx = scale(derivative(height, 0), 1 / width)
        self.fy = add(derivative(height, 1), scale(self.fx, -skew))
        self.g = add(ONE, mul(self.fx, self.fx), mul(self.fy, self.fy))
        self._moment_cache: dict[tuple[int, int, int, int], tuple[Poly, ...]] = {}
        self._lifting_cache: dict[
            tuple[int, int, int, int], tuple[tuple[Poly, ...], tuple[Poly, ...]]
        ] = {}

    def moments(self, axis: int, origin: int, k: int, cell: int) -> tuple[Poly, ...]:
        key = axis, origin, k, cell
        if key in self._moment_cache:
            return self._moment_cache[key]
        if cell < origin:
            return ({},) * 5
        lower, upper = F(origin, 2 * self.n), F(origin + 1, 2 * self.n)
        curvature = affine_legendre(k, axis, 4 * self.n, -(2 * origin + 1))
        result: list[Poly] = []
        for m in range(5):
            integrand = mul(monomial(m if axis == 0 else 0, m if axis == 1 else 0), curvature)
            integral = integrate_from(integrand, axis, lower)
            result.append(integral if cell == origin else trace(integral, axis, upper))
        self._moment_cache[key] = tuple(result)
        return self._moment_cache[key]

    def lifting(
        self, axis: int, origin: int, k: int, cell: int
    ) -> tuple[tuple[Poly, ...], tuple[Poly, ...]]:
        key = axis, origin, k, cell
        if key in self._lifting_cache:
            return self._lifting_cache[key]
        moments = self.moments(axis, origin, k, cell)
        if cell < origin:
            return ({},) * 3, ({},) * 3
        along = S if axis == 0 else T
        j0: Poly = {}
        j1: Poly = {}
        j2: Poly = {}
        for (i, j), c in self.height.items():
            power, cross = (i, j) if axis == 0 else (j, i)
            factor = monomial(0, cross, c) if axis == 0 else monomial(cross, 0, c)
            j0 = add(j0, mul(factor, moments[power]))
            if cross:
                factor = (
                    monomial(0, cross - 1, c * cross)
                    if axis == 0
                    else monomial(cross - 1, 0, c * cross)
                )
                j1 = add(j1, mul(factor, moments[power]))
                j2 = add(j2, mul(factor, moments[power + 1]))
        along_value = add(j0, scale(mul(self.height, moments[0]), -1))
        transverse = add(j2, scale(mul(along, j1), -1))
        first, second = (along_value, transverse) if axis == 0 else (transverse, along_value)
        vx = scale(first, 1 / self.width)
        vy = add(second, scale(vx, -self.skew))
        vz = add(mul(along, moments[0]), scale(moments[1], -1))
        v = vx, vy, vz
        dx = tuple(scale(derivative(x, 0), 1 / self.width) for x in v)
        dy = tuple(add(derivative(x, 1), scale(y, -self.skew)) for x, y in zip(v, dx, strict=True))
        tx = add(scale(mul(self.fx, dx[0]), -1), scale(mul(self.fy, dx[1]), -1), dx[2])
        ty = add(scale(mul(self.fx, dy[0]), -1), scale(mul(self.fy, dy[1]), -1), dy[2])
        p = (
            add(
                scale(mul(add(ONE, mul(self.fy, self.fy)), tx), -1), mul(mul(self.fx, self.fy), ty)
            ),
            add(
                mul(mul(self.fx, self.fy), tx), scale(mul(add(ONE, mul(self.fx, self.fx)), ty), -1)
            ),
            add(scale(mul(self.fx, tx), -1), scale(mul(self.fy, ty), -1)),
        )
        self._lifting_cache[key] = v, p
        return self._lifting_cache[key]

    def column(self, col: int, cs: int, ct: int) -> tuple[Poly, Poly, Poly]:
        if col < self.r - 1:
            i, j = divmod(col, self.n + 2)
            return (
                mul(
                    basis(self.n, 2, 1, i, axis=0, cell=cs), basis(self.n, 2, 1, j, axis=1, cell=ct)
                ),
                {},
                {},
            )
        if col == self.r - 1:
            return self.height, {}, {}
        if col < self.r + self.d:
            d = col - self.r
            if d < self.d - 3:
                i, j = divmod(d, 2 * self.n + 1)
                b = mul(
                    basis(self.n, 2, 2, i, axis=0, cell=cs), basis(self.n, 2, 2, j, axis=1, cell=ct)
                )
                return {}, mul(self.g, b), {}
            p = (scale(mul(self.g, self.fx), -1), scale(mul(self.g, self.fy), -1), self.g)
            return {}, p[d - (self.d - 3)], {}
        if col < self.r + self.d + self.p:
            j, rest = divmod(col - self.r - self.d, 32 * self.n)
            axis, rest = divmod(rest, 16 * self.n)
            origin, k = divmod(rest, 8)
            v, p = self.lifting(axis, origin, k, cs if axis == 0 else ct)
            return v[j], p[j], {}
        i, j = divmod(col - self.r - self.d - self.p, self.n + 1)
        return (
            {},
            {},
            mul(basis(self.n, 1, 1, i, axis=0, cell=cs), basis(self.n, 1, 1, j, axis=1, cell=ct)),
        )

    def anchor_sources(self, field: int) -> list[int]:
        result = []
        for col in range(self.total - self.c):
            if any(trace(self.column(col, 0, ct)[field], 0, F(0)) for ct in range(2 * self.n)):
                result.append(col)
        return result


def terms(p: Poly) -> list[list[int | str]]:
    return [[i, j, str(c.numerator), str(c.denominator)] for (i, j), c in sorted(p.items()) if c]


def make_case(n: int, curved: bool) -> dict[str, object]:
    height = (
        {(i, j): F(i + 2 * j + 1, i + j + 1) for i in range(5) for j in range(5 - i)}
        if curved
        else {}
    )
    x0, width, skew = (F(2, 7), F(4, 5), F(3, 20)) if curved else (F(0), F(1), F(0))
    oracle = Oracle(n, x0, width, skew, height)
    selections: set[tuple[int, int, int, str]] = set()
    cells = [(0, 0), (2 * n - 1, 2 * n - 1)]
    for cs, ct in cells:
        selections.add((-1, cs, ct, "g"))
        for col in [
            0,
            oracle.r - 2,
            oracle.r - 1,
            oracle.r,
            oracle.r + oracle.d - 4,
            *range(oracle.r + oracle.d - 3, oracle.r + oracle.d),
            oracle.total - oracle.c,
            oracle.total - 1,
        ]:
            for field in ("R", "P", "C"):
                selections.add((col, cs, ct, field))
    for j in range(3):
        for axis in range(2):
            for origin in [0, 2 * n - 1]:
                for k in [0, 7]:
                    col = oracle.r + oracle.d + j * 32 * n + axis * 16 * n + 8 * origin + k
                    for cs, ct in cells:
                        for field in ("R", "P"):
                            selections.add((col, cs, ct, field))
    selected = []
    for col, cs, ct, field in sorted(selections):
        p = oracle.g if field == "g" else oracle.column(col, cs, ct)[("R", "P", "C").index(field)]
        selected.append(
            {
                "column": col,
                "cell": [cs, ct],
                "field": field,
                "terms": terms(pullback(p, n, cs, ct)),
            }
        )
    ar, ap = oracle.anchor_sources(0), oracle.anchor_sources(1)
    # Preserve actual appended anchor ordinals and test extension away from s=0.
    for offset, (field, sources) in enumerate([("R", ar), ("P", ap)]):
        if not sources:
            continue
        base = oracle.total + (len(ar) if offset else 0)
        for source in [sources[0], sources[-1]]:
            col = base + sources.index(source)
            for ct in [0, 2 * n - 1]:
                p = mul(
                    add(ONE, scale(S, -1)), trace(oracle.column(source, 0, ct)[offset], 0, F(0))
                )
                selected.append(
                    {
                        "column": col,
                        "cell": [2 * n - 1, ct],
                        "field": field,
                        "terms": terms(pullback(p, n, 2 * n - 1, ct)),
                    }
                )

    def rational(x: F) -> list[str]:
        return [str(x.numerator), str(x.denominator)]

    return {
        "id": ("mixed-quartic" if curved else "flat") + f"-n{n}",
        "n": n,
        "x0": rational(x0),
        "width": rational(width),
        "skew": rational(skew),
        "height": terms(height),
        "original_columns": oracle.total,
        "selected": sorted(selected, key=lambda x: (x["column"], x["cell"], x["field"])),
        "anchor_r_sources": ar,
        "anchor_p_sources": ap,
    }


def self_check() -> None:
    # Hand-derived flat k=0: v_z=s^2/2, p_x=-s; t family swaps x/y.
    oracle = Oracle(1, F(0), F(1), F(0), {})
    vs, ps = oracle.lifting(0, 0, 0, 0)
    vt, pt = oracle.lifting(1, 0, 0, 0)
    assert vs[2] == monomial(2, 0, F(1, 2)) and ps[0] == scale(S, -1)
    assert vt[2] == monomial(0, 2, F(1, 2)) and pt[1] == scale(T, -1)
    # Completed first half interval: b=s/2-1/8, with nonzero carry.
    v, p = oracle.lifting(0, 0, 0, 1)
    assert v[2] == add(scale(S, F(1, 2)), scale(ONE, F(-1, 8)))
    assert p[0] == scale(ONE, F(-1, 2))
    assert legendre(2) == add(scale(monomial(2, 0), F(3, 2)), scale(ONE, F(-1, 2)))
    assert add(*(basis(4, 2, 1, i, axis=0, cell=3) for i in range(6))) == ONE


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--n", type=int, nargs="+", choices=[1, 2, 4], default=[1, 2, 4])
    args = parser.parse_args()
    self_check()
    cases = []
    for n in args.n:
        for curved in [False, True]:
            cases.append(make_case(n, curved))
            print(f"generated {cases[-1]['id']}", flush=True)
    output = {"schema": "uvlas.pg-coefficient-reference/1", "cases": cases}
    args.output.write_bytes(
        (json.dumps(output, ensure_ascii=False, separators=(",", ":")) + "\n").encode()
    )


if __name__ == "__main__":
    main()
