# SPDX-License-Identifier: Apache-2.0
"""Generate the frozen inputs and observations for the V-AERO harmonic wing checks."""

import argparse
import hashlib
import json
import runpy
from decimal import Decimal as D
from decimal import localcontext
from fractions import Fraction as Q
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
RECORDS = "rd/aero/unsteady-reference/records/"
FINAL = RECORDS + "high-frequency-final-result.json"
INITIAL = RECORDS + "ideal-pitch-frequency-result.json"
GEOMETRY = "docs/verification/IM-03/finite-wing-drag-refinement-v1/reference.py"
POLICY = RECORDS + "ideal-pitch-frequency-assessment-policy.md"
EXTENSION = RECORDS + "high-frequency-final-space-core-protocol.md"
STEADY_GRID = "rd/aero/finite-wing/experiments/uvlm_steady_reference.py"
FREQUENCY_SOLVE = "rd/aero/unsteady-reference/experiments/uvlm_frequency_fast.py"
LOADS = "rd/aero/unsteady-reference/experiments/frequency_motion_columns.py"
SOURCES = {
    FINAL: "76fede7b065f0a71c69d3f181962d74bebd811dbf34140b6804c5c8dcd40e24f",
    INITIAL: "00999fd24040e4827353f1b9f8e5e569c90872d7c9f7dac6527128eeaac717af",
    POLICY: "96af7ca71b5b903495ebc10334365bbd93fcfa1759b60fcc51f0c24ecfda3679",
    EXTENSION: "33fdec27d13eabbbdded627449f7b2efff108f203302917733bb5b5e5c90a929",
    STEADY_GRID: "36c6c6d03b6ef00bfff11697b9036c0130c364b8f9ac4a48534646d3d71e8ad1",
    FREQUENCY_SOLVE: "cf882a9dcdb7de4789b2c6d5d53ac4ccb382ddd346de6cab6489892aaa7ba304",
    LOADS: "a9e53109dd1a2bee79a7b56afe7ed524cc0e969ddffd8ca0eecafe7e8ba975fc",
    GEOMETRY: "5a62eb6ff117b17b760d407d42fe1e9e9d011ff7202ee9053a4fdf958971a69e",
}
FREQUENCIES = tuple(map(D, ("0.05", "0.2", "0.8")))
MOTIONS = ("heave", "pitch")
OBSERVABLES = ("Fz", "My")
REFERENCE_BUDGET = Q("0.0025")
COMBINED_BUDGET = Q("0.01")


def checked_sources():
    result = {}
    for name, digest in SOURCES.items():
        raw = (ROOT / name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError(f"Frozen source has changed: {name}")
        result[name] = raw
    return result


def pair(values):
    if len(values) != 2 or any(type(value) not in (int, D) for value in values):
        raise ValueError("Expected two finite complex components")
    return tuple(map(Q, values))


def accepts(candidate, reference, uncertainty):
    """Test the strict budgets by comparing squared complex magnitudes."""
    reference = pair(reference)
    difference = [a - b for a, b in zip(pair(candidate), reference, strict=True)]
    norm_squared = sum(value**2 for value in reference)
    if norm_squared == 0 or type(uncertainty) is not D or uncertainty < 0:
        raise ValueError("Invalid reference norm or uncertainty")
    remaining = COMBINED_BUDGET - Q(uncertainty)
    return Q(uncertainty) < REFERENCE_BUDGET and (
        sum(value**2 for value in difference) < remaining**2 * norm_squared
    )


def observations(sources):
    final = json.loads(sources[FINAL], parse_float=D)
    initial = json.loads(sources[INITIAL], parse_float=D)
    if final["sources"][INITIAL] != SOURCES[INITIAL]:
        raise ValueError("Final record does not cite the initial record")
    fields = ("k", "motion", "observable")
    result = []
    failed = set()
    keys = product(FREQUENCIES, MOTIONS, OBSERVABLES)
    for key, before, after in zip(keys, initial["records"], final["records"], strict=True):
        if any(tuple(row[field] for field in fields) != key for row in (before, after)):
            raise ValueError(f"Unexpected observation order at {key}")
        reference = after["reference_limit"]
        uncertainty = after["reference_relative_uncertainty"]
        if (reference, uncertainty) != (
            before["reference_limit"],
            before["reference_relative_uncertainty"],
        ):
            raise ValueError(f"Reference changed between records at {key}")
        if not accepts(after["candidate"], reference, uncertainty):
            raise ValueError(f"Final finite candidate misses the frozen budget at {key}")
        initial_passed = accepts(before["candidate"], reference, uncertainty)
        if not initial_passed:
            failed.add(key)
        chord, core = (256, D("1e-6")) if key[0] == FREQUENCIES[-1] else (128, D("1e-5"))
        if (
            after["candidate_chordwise_panels"],
            after["candidate_spanwise_panels"],
            after["candidate_core"],
        ) != (chord, 192, core):
            raise ValueError(f"Wrong finite discretization at {key}")
        result.append(
            {
                **dict(zip(fields, key, strict=True)),
                "chord_panels": chord,
                "span_panels": 192,
                "core_m": core,
                "reference_limit": reference,
                "reference_relative_uncertainty": uncertainty,
                "frozen_finite_value": after["candidate"],
                "initial_finite_value": before["candidate"],
                "initial_passed": initial_passed,
            }
        )
    if failed != {(FREQUENCIES[-1], motion, "My") for motion in MOTIONS}:
        raise ValueError(f"Unexpected initial failures: {failed}")
    return result


def reflect(points, width):
    for start in range(0, len(points), width):
        for left in range(width // 2):
            x, y, z = points[start + width - 1 - left]
            reflected = (x, -y, z)
            original = points[start + left]
            if any(abs(a - b) > D("1e-14") for a, b in zip(original, reflected, strict=True)):
                raise ValueError("UVLM grid points do not agree with their reflection")
            points[start + left] = reflected


def geometries():
    helpers = runpy.run_path(str(ROOT / GEOMETRY))
    frozen = helpers["import_path"]("frozen_geometry", helpers["checked"](helpers["CASE"]))
    result = []
    with localcontext() as context:
        context.prec = 75
        sine, cosine = frozen.trig(35)
        for chord in (128, 256):
            _, _, vertices, collocations = helpers["grids"](frozen, chord, 192, 0, sine / cosine)
            # Keep the physical TE. The steady study's closure is not part of this grid.
            del vertices[(chord + 1) * 193 :]
            reflect(vertices, 193)
            reflect(collocations, 192)
            result.append(
                {
                    "chord_panels": chord,
                    "span_panels": 192,
                    "vertices_m": vertices,
                    "collocations_m": collocations,
                }
            )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = {
        "id": "finite-wing-unsteady-v1",
        "sources": SOURCES,
        "reference_uncertainty_budget": REFERENCE_BUDGET,
        "combined_observation_budget": COMBINED_BUDGET,
        "relative_l2_residual_budget": D("1e-11"),
        "wake_length_m": 128,
        "quadrature_order": 16,
        "density_kg_m3": D("1.2"),
        "inflow_m_s": [20, 0, 0],
        "moment_origin_m": [D("0.375"), 0, 0],
        "reduced_frequency_chord_m": D("1.114285714285714"),
        "observations": observations(checked_sources()),
        "meshes": geometries(),
    }
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(payload, stream, default=float, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
