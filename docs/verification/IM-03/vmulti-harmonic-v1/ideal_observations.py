# SPDX-License-Identifier: Apache-2.0
"""Extract the frozen Fz and My observations without recomputing the ideal loads."""

import argparse
import json
from pathlib import Path

import numpy as np

from reference import CASES, DENSITY, FLOW, LENGTH, MOTIONS, OMEGA, ORIGIN, archived_ideal


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    archive = archived_ideal(args.root, args.bundle, CASES)
    records = []
    for case in CASES:
        source, digest, values = archive[case]
        observations = values["dual"][:, [2, 4], :]
        if not (abs(observations) > 0).all():
            raise ValueError(f"Zero frozen Fz or My scale for {case}")
        records.append({
            "case": case,
            "source": source,
            "source_sha256": digest,
            "loads": np.stack((observations.real, observations.imag), axis=-1).tolist(),
        })
    with args.output.open("x", encoding="utf-8", newline="\n") as output:
        json.dump({
            "schema": "uvlas.vmulti-harmonic-ideal-observations/1",
            "wake_length_m": LENGTH,
            "omega_rad_s": OMEGA,
            "velocity_m_s": FLOW.tolist(),
            "density_kg_m3": DENSITY,
            "moment_origin_m": ORIGIN.tolist(),
            "motions": MOTIONS,
            "quantities": ["Fz", "My"],
            "records": records,
        }, output, allow_nan=False, separators=(",", ":"))
        output.write("\n")


if __name__ == "__main__":
    main()
