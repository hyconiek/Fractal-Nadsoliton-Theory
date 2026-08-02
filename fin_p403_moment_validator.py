#!/usr/bin/env python3
"""Validate the P403 synthetic Jordan moment-reconstruction record.

This companion validator checks a frozen atom law, raw event schema,
timestamps, empirical sampling frequencies, and twelve unbiased moment
estimates.  It is deliberately separate from the custody/manifest validator
in ``fin_p403_jsr_validator.py``.  Neither validator can promote a synthetic
record to physical evidence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(protocol_path: Path, events_path: Path) -> dict[str, Any]:
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    atoms = protocol["atoms"]
    required = {"event_id", "run_id", "timestamp_tick", "atom", "node", "sign"}
    samples: list[tuple[int, float, int]] = []
    seen: set[str] = set()
    schema_ok = True
    atom_data_ok = True
    monotone_timestamps = True
    previous_tick = -1
    with events_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        schema_ok = required.issubset(set(reader.fieldnames or []))
        for row in reader:
            event_id = row["event_id"]
            if event_id in seen:
                schema_ok = False
            seen.add(event_id)
            tick = int(row["timestamp_tick"])
            monotone_timestamps &= tick > previous_tick
            previous_tick = tick
            index = int(row["atom"])
            node = float(row["node"])
            sign = int(row["sign"])
            if not 0 <= index < len(atoms):
                atom_data_ok = False
                continue
            reference = atoms[index]
            atom_data_ok &= math.isclose(
                node, float(reference["node"]), rel_tol=0.0, abs_tol=1e-15
            )
            atom_data_ok &= sign == int(reference["sign"])
            samples.append((index, node, sign))
    count = len(samples)
    if count == 0:
        raise ValueError("empty event record")
    total_variation = float(protocol["total_variation"])
    targets = np.asarray(protocol["target_moments"], dtype=float)
    values = np.empty((count, len(targets)), dtype=float)
    for row_index, (_, node, sign) in enumerate(samples):
        values[row_index] = [
            total_variation * sign * node**order
            for order in range(len(targets))
        ]
    estimates = np.mean(values, axis=0)
    standard_errors = np.std(values, axis=0, ddof=1) / math.sqrt(count)
    z_scores = np.abs(estimates - targets) / np.maximum(standard_errors, 1e-30)
    probabilities = np.asarray(
        [float(atom["probability"]) for atom in atoms], dtype=float
    )
    counts = np.bincount(
        np.asarray([sample[0] for sample in samples], dtype=int),
        minlength=len(atoms),
    )
    empirical = counts / count
    probability_z = np.abs(empirical - probabilities) / np.sqrt(
        np.maximum(probabilities * (1.0 - probabilities) / count, 1e-30)
    )
    synthetic = bool(protocol.get("synthetic_reference", False))
    return {
        "event_count": count,
        "schema_ok": bool(schema_ok),
        "atom_data_ok": bool(atom_data_ok),
        "unique_event_ids": len(seen) == count,
        "monotone_timestamps": bool(monotone_timestamps),
        "protocol_probability_sum": float(np.sum(probabilities)),
        "maximum_probability_z_score": float(np.max(probability_z)),
        "maximum_moment_z_score": float(np.max(z_scores)),
        "maximum_absolute_moment_error": float(np.max(np.abs(estimates - targets))),
        "moment_estimates": estimates.tolist(),
        "moment_standard_errors": standard_errors.tolist(),
        "raw_record_sha256": sha256_file(events_path),
        "synthetic_reference": synthetic,
        "physical_evidence_admitted": False,
        "validation_pass": bool(
            schema_ok
            and atom_data_ok
            and len(seen) == count
            and monotone_timestamps
            and math.isclose(float(np.sum(probabilities)), 1.0, abs_tol=2e-15)
            and float(np.max(z_scores)) < 6.0
            and float(np.max(probability_z)) < 6.0
        ),
        "boundary": (
            "A passing synthetic or self-generated record validates software "
            "and statistical reconstruction only. It is not laboratory evidence."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("protocol", type=Path)
    parser.add_argument("events", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = validate(args.protocol, args.events)
    output = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    else:
        print(output, end="")
    raise SystemExit(0 if result["validation_pass"] else 1)


if __name__ == "__main__":
    main()
