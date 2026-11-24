# SPDX-License-Identifier: Apache-2.0
"""Extract wake cost observations and keep only complete, reproducible mode groups."""

import csv
import sys
from collections import defaultdict
from pathlib import Path

PREPARATION_FIELDS = (
    "wake_panels",
    "contact_ns",
    "case_and_workspaces_ns",
    "contact_planned_bytes",
    "all_preparation_planned_bytes",
    "contact_resident_bytes",
    "preparation_resident_bytes",
    "process_peak_bytes",
)
FIELDS = (
    "wake_panels",
    "precision",
    "operation",
    "replay",
    "particles",
    "source_facet_visits",
    "far_wake_facets",
    "cancelled_analytic_pairs",
    "velocity_charged_steps",
    "elapsed_ns",
    "planned_bytes",
    "before_resident_bytes",
    "after_resident_bytes",
    "process_peak_bytes",
    "result_word_signature",
)


def main(directory):
    preparations, samples = [], []
    repeat = defaultdict(list)
    for panels in (2, 4, 8):
        data = directory.joinpath(f"wake-{panels}.log").read_text(encoding="utf-8-sig")
        assert "All tests passed" in data
        for row in csv.reader(line for line in data.splitlines() if line.startswith("WAKE_")):
            tag, *values = row
            assert int(values[0]) == panels
            match tag:
                case "WAKE_PREP":
                    assert len(values) == len(PREPARATION_FIELDS)
                    preparations.append(values)
                case "WAKE_COST":
                    mapped = dict(zip(FIELDS, values, strict=True))
                    samples.append(values)
                    group = tuple(values[:3])
                    repeat[group].append((int(mapped["replay"]), mapped["result_word_signature"]))
                case _:
                    raise ValueError(tag)
    assert len(preparations) == 3 and len(samples) == 90 and len(repeat) == 30
    for observations in repeat.values():
        indices, signatures = zip(*observations, strict=True)
        assert indices == (0, 1, 2) and len(set(signatures)) == 1
    for name, fields, data in [
        ("preparation.csv", PREPARATION_FIELDS, preparations),
        ("observations.csv", FIELDS, samples),
    ]:
        with directory.joinpath(name).open("w", encoding="utf-8", newline="") as target:
            writer = csv.writer(target, lineterminator="\n")
            writer.writerow(fields)
            writer.writerows(data)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
