# SPDX-License-Identifier: Apache-2.0
"""Independent one-panel finite-core known-flow load and derivative oracle."""

from decimal import Decimal as D, getcontext

getcontext().prec = 80
PI = D("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628")
EPS2 = D(1) / 64
RHO = D(5) / 4
POINTS = tuple(tuple(map(D, p)) for p in ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0)))
EDGES = tuple(zip(POINTS, POINTS[1:] + POINTS[:1], strict=True))
ORIGIN = (D(1) / 4, -D(1) / 8, D(1) / 2)
FORCE_SEED = tuple(map(D, (1, 2, 3)))
MOMENT_SEED = (D(1) / 8, -D(1) / 4, D(1) / 2)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b, strict=True))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b, strict=True))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b, strict=True))


def mul(a, scale):
    return tuple(x * scale for x in a)


def cross(a, b):
    return tuple(
        a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3] for i in range(3)
    )


def velocity(point, skip=-1):
    total = D(0)
    for i, (a, b) in enumerate(EDGES):
        if i == skip:
            continue
        dl = sub(b, a)
        ra, rb = sub(point, a), sub(point, b)
        z = cross(dl, ra)[2]
        integral = (
            dot(dl, ra) / (dot(ra, ra) + EPS2).sqrt() - dot(dl, rb) / (dot(rb, rb) + EPS2).sqrt()
        )
        total += z * integral / (4 * PI * (z * z + EPS2 * dot(dl, dl)))
    return total


def evaluate(translation, time, rate_direction):
    phase_direction = PI * (translation - time / 2)
    half_root = D(2).sqrt() / 2
    influence = velocity((D(3) / 4, D(1) / 2, D(0)))
    gamma = -(D(3) / 4 + half_root) / influence
    dgamma = half_root * phase_direction / influence
    force = moment = df = dm = (D(0), D(0), D(0))
    delta_x = (translation, D(0), D(0))
    edge_forces, edge_directions = [], []

    def accumulate(x, f, direction):
        nonlocal force, moment, df, dm
        force = add(force, f)
        moment = add(moment, cross(sub(x, ORIGIN), f))
        df = add(df, direction)
        dm = add(dm, add(cross(delta_x, f), cross(sub(x, ORIGIN), direction)))

    for i, ((a, b), s, c) in enumerate(zip(EDGES, (0, 1, 0, 1), (1, 0, -1, 0), strict=True)):
        center = mul(add(a, b), D(1) / 2)
        edge = sub(b, a)
        coeff = velocity(center, i)
        effective = 1 + s + gamma * coeff
        deffective = c * phase_direction + dgamma * coeff
        f = mul(cross((0, 0, effective), edge), RHO * gamma)
        direction = mul(cross((0, 0, dgamma * effective + gamma * deffective), edge), RHO)
        edge_forces.append(f)
        edge_directions.append(direction)
        accumulate(center, f, direction)

    for center in ((D(1) / 3, D(2) / 3, D(0)), (D(2) / 3, D(1) / 3, D(0))):
        accumulate(center, (0, 0, RHO * D(3) / 16), (0, 0, RHO * rate_direction / 2))
    objective_direction = dot(FORCE_SEED, df) + dot(MOMENT_SEED, dm) + edge_directions[0][0]
    return gamma, dgamma, edge_forces, edge_directions, force, moment, df, dm, objective_direction


if __name__ == "__main__":
    for v in evaluate(D(1) / 4, D(1) / 4, D(1) / 8):
        print(v)
    print("time dual", evaluate(D(0), D(1), D(0))[-1])
    print("rate dual", evaluate(D(0), D(0), D(1))[-1])
