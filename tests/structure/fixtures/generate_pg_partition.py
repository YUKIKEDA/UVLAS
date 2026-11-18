"""Independent PG partition fixture: vertical sections, no product code imports."""

import argparse
import json
from fractions import Fraction as F
from itertools import pairwise
from pathlib import Path


def integral(domain, lo, hi, bottom, top):
    edges = list(zip(domain, domain[1:] + domain[:1], strict=True))
    breaks = {lo, hi}
    for a, b in edges:
        if lo < a[0] < hi:
            breaks.add(a[0])
        if a[1] != b[1]:
            for y in (bottom, top):
                x = a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
                if lo < x < hi and min(a[0], b[0]) <= x <= max(a[0], b[0]):
                    breaks.add(x)
    cuts = sorted(breaks)
    answer = [F(0), F(0), F(0)]
    for left, right in pairwise(cuts):
        mid = (left + right) / 2
        active = []
        for a, b in edges:
            if min(a[0], b[0]) < mid < max(a[0], b[0]):
                slope = (b[1] - a[1]) / (b[0] - a[0])
                intercept = a[1] - slope * a[0]
                active.append((slope * mid + intercept, slope, intercept))
        active.sort()
        assert len(active) % 2 == 0
        for low, high in zip(active[::2], active[1::2], strict=True):
            if high[0] <= bottom or low[0] >= top:
                continue
            lower = (F(0), bottom) if low[0] < bottom else low[1:]
            upper = (F(0), top) if high[0] > top else high[1:]

            def values(x, lower=lower, upper=upper):
                y0, y1 = lower[0] * x + lower[1], upper[0] * x + upper[1]
                return y1 - y0, x * (y1 - y0), (y1 * y1 - y0 * y0) / 2

            samples = [values(x) for x in (left, mid, right)]
            for k in range(3):
                answer[k] += (
                    (right - left) * (samples[0][k] + 4 * samples[1][k] + samples[2][k]) / 6
                )
    return answer


def pair(x):
    return [str(x.numerator), str(x.denominator)]


def generate():
    tiny = F(1, 10**400)
    regions = [
        ("square", [(0, 0), (1, 0), (1, 1), (0, 1)], [F(1), F(1, 2), F(1, 2)]),
        ("trapezoid", [(0, 0), (1, 0), (F(3, 4), 1), (F(1, 4), 1)], [F(3, 4), F(3, 8), F(1, 3)]),
        (
            "two-components",
            [
                (0, 0),
                (1, 0),
                (1, F(1, 8)),
                (F(1, 8), F(1, 8)),
                (F(1, 8), F(3, 8)),
                (1, F(3, 8)),
                (1, F(1, 2)),
                (0, F(1, 2)),
            ],
            [F(9, 32), F(65, 512), F(9, 128)],
        ),
        (
            "tiny-positive",
            [(0, 0), (1, 0), (1, tiny), (0, tiny)],
            [tiny, tiny / 2, tiny * tiny / 2],
        ),
    ]
    result = []
    for name, raw, literal in regions:
        domain = [(F(x), F(y)) for x, y in raw]
        for n in (1, 2, 4):
            cells = []
            total = [F(0), F(0), F(0)]
            for i in range(2 * n):
                for j in range(2 * n):
                    value = integral(
                        domain, F(i, 2 * n), F(i + 1, 2 * n), F(j, 2 * n), F(j + 1, 2 * n)
                    )
                    total = [a + b for a, b in zip(total, value, strict=True)]
                    cells.append({"cell": [i, j], "moments": [pair(x) for x in value]})
            assert total == literal, (name, n, total, literal)
            result.append(
                {
                    "name": name,
                    "n": n,
                    "domain": [[pair(x), pair(y)] for x, y in domain],
                    "total_moments": [pair(x) for x in literal],
                    "cells": cells,
                }
            )
    return {
        "schema": "uvlas.pg-partition-fixture/1",
        "method": "independent vertical-section exact Simpson",
        "tolerance": "exact equality; no floating-point conversion",
        "cases": result,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = generate()
    args.output.write_bytes((json.dumps(output, indent=2) + "\n").encode("utf-8"))
    print(f"cases={len(output['cases'])}; cell moments=1008")
