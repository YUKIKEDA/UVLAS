"""Bound the actual rational Gauss approximation on the cut-parabola fixture."""

from decimal import ROUND_HALF_EVEN, localcontext
from fractions import Fraction as F
from functools import cache
from math import isqrt

from reference import unit_gauss

NODE_SCALE = 10**48
MOMENT_SCALE = 2**272
ROOT_SCALE = 2**352
PAIR_TERMS = (
    ("RR", 0, 0, F(1)),
    ("RP", 1, 3, F(1)),
    ("RC", 2, 0, F(1, 2)),
    ("PP", 2, 6, F(1)),
    ("PC", 3, 3, F(1, 2)),
    ("CC", 4, 0, F(1, 4)),
)


@cache
def rational_rule(order: int) -> tuple[tuple[F, F], ...]:
    with localcontext() as context:
        context.prec = 85
        half = []
        for x, w in unit_gauss(order)[: order // 2]:
            nx, nw = (
                int((v * NODE_SCALE).to_integral_value(rounding=ROUND_HALF_EVEN)) for v in (x, w)
            )
            half.append((nx, nw))
    centre_weight = NODE_SCALE - 2 * sum(weight for _, weight in half)
    nodes = [(NODE_SCALE - x, w) for x, w in half]
    nodes.append((NODE_SCALE // 2, centre_weight))
    nodes.extend(reversed(half))
    if any(not 0 < x < NODE_SCALE or w <= 0 for x, w in nodes):
        raise ValueError("Rational rule needs interior nodes with positive weights")
    return tuple((F(x, NODE_SCALE), F(w, NODE_SCALE)) for x, w in nodes)


def enclose(lower: F, upper: F) -> tuple[int, int]:
    low = lower * MOMENT_SCALE
    high = upper * MOMENT_SCALE
    return low.numerator // low.denominator, -((-high.numerator) // high.denominator)


def product(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    lo_a, hi_a = left
    lo_b, hi_b = right
    products = lo_a * lo_b, lo_a * hi_b, hi_a * lo_b, hi_a * hi_b
    return min(products) // MOMENT_SCALE, -((-max(products)) // MOMENT_SCALE)


def powers(x: F, maximum: int) -> list[tuple[int, int]]:
    result = []
    term = F(1)
    for _ in range(maximum + 1):
        result.append(enclose(term, term))
        term *= x
    return result


def depth_moments() -> list[F]:
    layers = ((F(-1), F(-1, 3)), (F(-1, 3), F(1)))
    return [
        sum(
            (hi - lo) * weight * (lo + (hi - lo) * node) ** order
            for lo, hi in layers
            for node, weight in rational_rule(3)
        )
        for order in range(6)
    ]


def value_bounds(s: F, weight: F, depth: list[F]) -> dict[str, tuple[int, int]]:
    half_width = F(1, 4) + s / 8
    g = 1 + s * s / 64
    scaled = g * ROOT_SCALE**2
    q = isqrt(scaled.numerator // scaled.denominator)
    roots = F(q, ROOT_SCALE), F(q + 1, ROOT_SCALE)
    depth_scale = [half_width ** (k + 1) * value for k, value in enumerate(depth)]
    result = {}
    for field, k, half_power, factor in PAIR_TERMS:
        first = weight * factor * depth_scale[k]
        second = -weight * factor * depth_scale[k + 1] / (8 * g)
        # Reduce g^(1/2-a/2) and g^(-1-a/2) to rational + rational*sqrt(g).
        if half_power % 2:
            constant = first / g ** ((half_power - 1) // 2)
            radical = second / g ** ((half_power + 1) // 2)
        else:
            constant = second / g ** (half_power // 2)
            radical = first / g ** (half_power // 2)
        bounds = tuple(constant + radical * root for root in roots)
        result[field] = enclose(min(bounds), max(bounds))
    return result


def sampled_bounds(cell, triangles, keys, order: int):
    moments = {key: [0, 0] for key in keys}
    s_cell, t_cell = cell
    max_s = max(a for _, a, _ in keys)
    max_t = max(b for _, _, b in keys)
    depth = depth_moments()
    for s0, t0, s1, t1, s2, t2 in triangles:
        ds_u, dt_u = s1 - s0, t1 - t0
        ds_v, dt_v = s2 - s0, t2 - t0
        determinant = ds_u * dt_v - ds_v * dt_u
        if determinant <= 0:
            raise ValueError("Nonpositive triangle orientation")
        for u, weight_u in rational_rule(order):
            for raw_v, weight_v in rational_rule(order):
                v = (1 - u) * raw_v
                s = s0 + u * ds_u + v * ds_v
                t = t0 + u * dt_u + v * dt_v
                s_values = powers(4 * s - (2 * s_cell + 1), max_s)
                t_values = powers(4 * t - (2 * t_cell + 1), max_t)
                factors = value_bounds(s, determinant * (1 - u) * weight_u * weight_v, depth)
                in_s = {
                    (field, a): product(factor, s_power)
                    for field, factor in factors.items()
                    for a, s_power in enumerate(s_values)
                }
                for key in keys:
                    field, a, b = key
                    lower, upper = product(in_s[field, a], t_values[b])
                    moments[key][0] += lower
                    moments[key][1] += upper
    return moments


if __name__ == "__main__":
    print(f'constexpr std::string_view rule_denominator = "{NODE_SCALE}";')
    for order in (17, 3):
        print(f"constexpr std::array<std::array<std::string_view, 2>, {order}> g{order}{{{{")
        for coordinate, weight in rational_rule(order):
            print(f'    {{"{coordinate * NODE_SCALE}", "{weight * NODE_SCALE}"}},')
        print("}};")
