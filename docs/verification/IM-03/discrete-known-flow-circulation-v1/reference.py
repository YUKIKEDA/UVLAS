# SPDX-License-Identifier: Apache-2.0
"""Independent closed segment integrals and 2-by-2 circulation solve."""

from decimal import Decimal as D, getcontext

getcontext().prec = 80
PI = D("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628")
EPSILON_SQ = D(1) / 64
POINTS = ((0, 0), (0, 1), (1, 1), (1, 0), (0, 2), (0, 3), (1, 3), (1, 2))
SEED = (D(1) / 4, -D(5) / 32)


def ring(point, panel):
    """Vertical velocity of one finite-core contour from the straight-edge integrals."""
    total = D(0)
    for a, b in zip(panel, panel[1:] + panel[:1], strict=True):
        dx, dy = (j - i for i, j in zip(a, b, strict=True))
        ux, uy = (j - i for i, j in zip(a, point, strict=True))
        vx, vy = (j - i for i, j in zip(b, point, strict=True))
        axis_cross = dx * uy - dy * ux
        begin = (dx * ux + dy * uy) / (ux * ux + uy * uy + EPSILON_SQ).sqrt()
        end = (dx * vx + dy * vy) / (vx * vx + vy * vy + EPSILON_SQ).sqrt()
        total += axis_cross * (begin - end) / (
            4 * PI * (axis_cross**2 + EPSILON_SQ * (dx**2 + dy**2))
        )
    return total


def evaluation(free_dx=D(0), point7_dx=D(0)):
    geometry = [[D(x), D(y)] for x, y in POINTS]
    geometry[6][0] += point7_dx
    # The material row uses facet 203 with fixed weight 3/8 on node 7.
    queries = ((D(3) / 4 + free_dx, D(1) / 2), (D(5) / 8 + 3 * point7_dx / 8, D(19) / 8))
    a, b, c, d = (ring(q, segment) for q in queries for segment in (geometry[:4], geometry[4:]))
    det = a * d - b * c

    def solve(u, v):
        return (u * d - b * v) / det, (a * v - u * c) / det

    sqrt2 = D(2).sqrt()
    # Phase values pi/4 and 0, a=2 m/s, e=(0,0,1). At a perturbed point
    # the affine Taylor term gives the same derivative as the sine.
    relative0 = D("0.5") + sqrt2 + 2 * PI * sqrt2 * free_dx
    relative1 = D("0.25") + 4 * PI * (3 * point7_dx / 8)
    circulation = solve(-relative0, -relative1)
    translation_and_time = solve(-PI * sqrt2 / 4, -PI / 2)
    transpose = ((d * SEED[0] - c * SEED[1]) / det, (-b * SEED[0] + a * SEED[1]) / det)
    objective = sum(q * x for q, x in zip(SEED, circulation, strict=True))
    time_dual = transpose[0] * PI * sqrt2 / 2 + transpose[1] * PI
    return circulation, translation_and_time, transpose, objective, time_dual


def direction_pair(h):
    return (
        (evaluation(free_dx=h)[3] - evaluation(free_dx=-h)[3]) / (2 * h),
        (evaluation(point7_dx=h)[3] - evaluation(point7_dx=-h)[3]) / (2 * h),
    )


if __name__ == "__main__":
    print("primal, tangent, transpose, objective, time dual:", evaluation())
    prior = direction_pair(D("1e-20"))
    reference = direction_pair(D("1e-22"))
    assert max(abs(x - y) for x, y in zip(prior, reference, strict=True)) < D("1e-37")
    print("free x and material point 7 x duals:", reference)