"""Generate product Gauss node candidates with Decimal Newton, not the test Sturm oracle."""
from decimal import Decimal, localcontext, ROUND_FLOOR
from math import cos, pi


def legendre(n: int, z: Decimal) -> tuple[Decimal, Decimal]:
    a, b = Decimal(1), z
    for k in range(2, n + 1):
        a, b = b, ((2 * k - 1) * z * b - (k - 1) * a) / k
    return b, n * (z * b - a) / (z * z - 1)


with localcontext() as context:
    context.prec = 300
    lines = [
        "#pragma once", "#include <array>", "#include <string_view>",
        "namespace uvlas::aero::detail::p1 {",
        "// Candidate node brackets; exact root ownership is independently verified in tests.",
        "struct GaussRoot { std::string_view lower, upper; };",
    ]
    for n in (8, 12):
        roots = []
        for i in range(1, n + 1):
            z = Decimal(cos(pi * (i - .25) / (n + .5)))
            for _ in range(32):
                value, derivative = legendre(n, z)
                z -= value / derivative
            roots.append((z + 1) / 2)
        roots.sort()
        lines.append(f"inline constexpr std::array<GaussRoot,{n}> gauss_{n}" + "{{")
        for root in roots:
            lower = int((root * (1 << 256)).to_integral_value(rounding=ROUND_FLOOR))
            lines.append('    {"' + str(lower) + '","' + str(lower + 1) + '"},')
        lines.append("}};")
    lines.append("} // namespace uvlas::aero::detail::p1")
    print("\n".join(lines))