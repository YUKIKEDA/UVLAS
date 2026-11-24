# SPDX-License-Identifier: Apache-2.0
"""Independent decimal solves for the fixed discrete circulation sensitivity fixtures."""

from contextlib import redirect_stdout
from decimal import Decimal as D, localcontext
from io import StringIO
from pathlib import Path
from runpy import run_path

# The predecessor evaluates each coefficient from direct physical integrals.
with redirect_stdout(StringIO()):
    model = run_path(
        str(Path(__file__).resolve().parents[1] / "discrete-circulation-v1/reference.py")
    )


def solve(a, b, c, d, rhs):
    """Solve a 2-by-2 system directly by the determinant formula."""
    x, y = rhs
    determinant = a * d - b * c
    return ((d * x - b * y) / determinant, (a * y - c * x) / determinant)


with localcontext() as ctx:
    ctx.prec = 60
    a, b, c, d = (model[key] for key in ("a00", "a01", "a10", "a11"))
    residual = (-D(1) / 8 + D(2) ** -57, D(3) / 16 - D(2) ** -58)
    seed = (D(1) / 4 + D(2) ** -57, -D(5) / 32 - D(2) ** -59)
    tangent = solve(a, b, c, d, tuple(-value for value in residual))
    transpose = solve(a, c, b, d, seed)
    print("attached wake tangent", -D(1) / (8 * (a + model["front"])))
    print("tangent", *tangent)
    print("transpose", *transpose)
    print("duality", sum(x * v + y * w for x, y, v, w in zip(tangent, transpose, seed, residual)))
