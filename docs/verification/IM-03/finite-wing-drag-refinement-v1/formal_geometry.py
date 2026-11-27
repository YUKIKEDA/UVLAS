# SPDX-License-Identifier: Apache-2.0
"""Generate the independent geometry input for the formal-grid preparation check."""

import argparse
import hashlib
import json
import runpy
from decimal import localcontext
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    reference = Path(__file__).with_name("reference.py")
    digest = hashlib.sha256(reference.read_bytes()).hexdigest()
    if digest != "5a62eb6ff117b17b760d407d42fe1e9e9d011ff7202ee9053a4fdf958971a69e":
        raise ValueError(f"Frozen refinement generator has changed: {reference}")
    helpers = runpy.run_path(str(reference))
    frozen = helpers["import_path"]("frozen_geometry", helpers["checked"](helpers["CASE"]))
    helpers["checked"](helpers["POLICY"])

    cases = []
    with localcontext() as context:
        context.prec = 75
        sine, cosine = frozen.trig(35)
        sweep = sine / cosine
        for dihedral in (0, 5):
            sine, cosine = frozen.trig(dihedral)
            chords, spans, vertices, queries = helpers["grids"](
                frozen, 128, 192, sine / cosine, sweep
            )
            cases.append(
                {
                    "dihedral_degrees": dihedral,
                    "chord_panels": 128,
                    "span_panels": 192,
                    "chord_fractions": chords,
                    "half_span_fractions": spans,
                    "vertices_m": vertices,
                    "collocations_m": queries,
                }
            )
    payload = {
        "id": "finite-wing-drag-refinement-v1",
        "scope": "geometry_only",
        "source": helpers["SOURCE_SHA256"],
        "closure_m": float(helpers["CLOSURE"]),
        "core_m": float(helpers["RADIUS"]),
        "cases": cases,
    }
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(payload, stream, default=float, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
