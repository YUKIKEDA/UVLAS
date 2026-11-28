# SPDX-License-Identifier: Apache-2.0
"""Verify the published static Case and observation reference without UVLAS."""

import argparse
from decimal import Decimal
from hashlib import sha256
import json
from pathlib import Path

CONDITION_FIELDS = {
    "schema", "units", "model", "flow", "analysis", "numerics", "initial", "restart",
}


def unique_object(pairs):
    result = dict(pairs)
    if len(result) != len(pairs):
        raise ValueError("Duplicate JSON member")
    return result


def reject_constant(value):
    raise ValueError(f"Nonfinite JSON value: {value}")


def read_json(data):
    return json.loads(
        data, parse_int=Decimal, parse_float=Decimal, parse_constant=reject_constant,
        object_pairs_hook=unique_object,
    )


def header(tag, count):
    return tag + str(count).encode("ascii") + b":"


def encode(value):
    if isinstance(value, dict):
        keys = sorted(value, key=lambda key: key.encode("utf-8"))
        return header(b"o", len(keys)) + b"".join(encode(k) + encode(value[k]) for k in keys)
    if isinstance(value, list):
        return header(b"a", len(value)) + b"".join(map(encode, value))
    if isinstance(value, str):
        data = value.encode("utf-8")
        return header(b"s", len(data)) + data
    if isinstance(value, Decimal):
        if value.is_zero():
            return b"d1:0"
        sign, digits, exponent = value.as_tuple()
        coefficient = "".join(map(str, digits)).rstrip("0")
        exponent += len(digits) - len(coefficient)
        data = f"{'-' if sign else ''}{coefficient}e{exponent}".encode("ascii")
        return header(b"d", len(data)) + data
    return {None: b"n", True: b"t", False: b"f"}[value]


def condition(case):
    projection = {key: value for key, value in case.items() if key in CONDITION_FIELDS}
    return sha256(b"uvlas.observation-condition/1\n" + encode(projection)).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    here = Path(__file__).resolve().parent
    parser.add_argument("--case", type=Path, default=here / "static-observation-case.json")
    parser.add_argument("--reference", type=Path, default=here / "parametric-observations.json")
    args = parser.parse_args()
    case = read_json(args.case.read_bytes())
    raw_reference = args.reference.read_bytes()
    reference = read_json(raw_reference)
    if reference["schema"] != "uvlas.observation-reference/1":
        raise ValueError("Unexpected observation reference schema")
    if encode({"x": [Decimal("1.250"), "a"]}) != b"o1:s1:xa2:d6:125e-2s1:a":
        raise ValueError("Declared encoding example differs")
    digest = condition(case)
    if digest != reference["condition_sha256"]:
        raise ValueError(f"Case condition SHA-256 differs: {digest}")
    if digest != condition(dict(reversed(tuple(case.items())))):
        raise ValueError("Case projection depends on input key order")
    reference_hash = sha256(raw_reference).hexdigest()
    fields = {row["id"]: [part["name"] for part in row["components"]]
              for row in reference["observations"]}
    for observation in case["observations"]:
        declared = observation["comparison"]["reference"]
        if declared["source"]["sha256"] != reference_hash:
            raise ValueError("Reference bytes differ from the Case declaration")
        if declared["source"]["path"] != args.reference.name:
            raise ValueError("Reference path differs from the Case declaration")
        if declared["columns"] != fields.pop(observation["id"]):
            raise ValueError("Reference columns differ from the Case declaration")
    if fields:
        raise ValueError("The Case did not declare all reference observations")
    print(f"Case condition SHA-256  {digest}")
    print(f"Reference file SHA-256  {reference_hash}")


if __name__ == "__main__":
    main()
