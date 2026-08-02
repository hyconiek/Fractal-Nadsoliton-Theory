#!/usr/bin/env python3
"""Structural validator for the P403 Jordan Sampling Realization package.

The validator checks syntax, file integrity, atom typing, and custody fields.
It cannot certify that an apparatus ran or that metadata are truthful.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import tempfile


def digest(path: Path) -> str:
    value = hashlib.sha256()
    value.update(path.read_bytes())
    return value.hexdigest()


def validate_legacy_reference(protocol_path: Path, events_path: Path) -> dict:
    """Validate the archived Release-10.34 synthetic reference format.

    This compatibility path never admits physical evidence.
    """
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    atoms = {int(atom["atom"]): atom for atom in protocol["atoms"]}
    with events_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        schema_ok = reader.fieldnames == protocol["required_event_fields"]
        rows = list(reader)
    atom_data_ok = True
    identifiers = []
    timestamps = []
    counts = {index: 0 for index in atoms}
    moment_sum = [0.0] * 12
    moment_square_sum = [0.0] * 12
    total_variation = float(protocol["total_variation"])
    for row in rows:
        try:
            atom_id = int(row["atom"])
            atom = atoms[atom_id]
            node = float(row["node"])
            sign = int(row["sign"])
            atom_data_ok &= math.isclose(node, float(atom["node"]), abs_tol=1e-15)
            atom_data_ok &= sign == int(atom["sign"])
            identifiers.append(row["event_id"])
            timestamps.append(int(row["timestamp_tick"]))
            counts[atom_id] += 1
            for order in range(12):
                score = total_variation * sign * node**order
                moment_sum[order] += score
                moment_square_sum[order] += score * score
        except (KeyError, ValueError):
            atom_data_ok = False
    count = len(rows)
    probability_z = []
    for index, atom in atoms.items():
        probability = float(atom["probability"])
        standard = math.sqrt(max(count * probability * (1 - probability), 1e-30))
        probability_z.append(abs(counts[index] - count * probability) / standard)
    estimates = [value / count for value in moment_sum]
    standard_errors = []
    for total, square in zip(moment_sum, moment_square_sum):
        variance = max((square - total * total / count) / (count - 1), 0.0)
        standard_errors.append(math.sqrt(variance / count))
    moment_z = [
        abs(estimate - float(target)) / max(error, 1e-30)
        for estimate, target, error in zip(
            estimates, protocol["target_moments"], standard_errors
        )
    ]
    result = {
        "event_count": count,
        "schema_ok": schema_ok,
        "atom_data_ok": atom_data_ok,
        "unique_event_ids": len(identifiers) == len(set(identifiers)),
        "monotone_timestamps": all(left < right for left, right in zip(timestamps, timestamps[1:])),
        "protocol_probability_sum": sum(float(atom["probability"]) for atom in atoms.values()),
        "maximum_probability_z_score": max(probability_z),
        "maximum_moment_z_score": max(moment_z),
        "maximum_absolute_moment_error": max(abs(estimate - float(target)) for estimate, target in zip(estimates, protocol["target_moments"])),
        "moment_estimates": estimates,
        "moment_standard_errors": standard_errors,
        "raw_record_sha256": digest(events_path),
        "synthetic_reference": bool(protocol.get("synthetic_reference")),
        "physical_evidence_admitted": False,
    }
    result["validation_pass"] = (
        result["schema_ok"]
        and result["atom_data_ok"]
        and result["unique_event_ids"]
        and result["monotone_timestamps"]
        and result["maximum_probability_z_score"] < 6
        and result["maximum_moment_z_score"] < 6
    )
    result["boundary"] = "Passing synthetic data validate software only, never physical evidence."
    return result


def validate(spec_path: Path, manifest_path: Path, events_path: Path | None = None) -> dict:
    if events_path is None:
        return validate_legacy_reference(spec_path, manifest_path)
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    atoms = {int(atom["atom_id"]): atom for atom in spec["atoms"]}
    failures: list[str] = []
    missing = [key for key in spec["manifest_required"] if key not in manifest]
    if missing:
        failures.append(f"missing manifest fields: {missing}")
    if manifest.get("raw_events_sha256") != digest(events_path):
        failures.append("raw event SHA-256 mismatch")
    roles = [manifest.get(name) for name in ("provider", "registrar", "analyst")]
    if len(set(roles)) != 3 or any(not role for role in roles):
        failures.append("provider, registrar, and analyst must be nonempty and distinct")
    if manifest.get("holdout_frozen_utc", "") >= manifest.get("unblinding_utc", ""):
        failures.append("hold-out must be frozen before unblinding")
    probability_sum = sum(float(atom["sampling_probability"]) for atom in atoms.values())
    if not math.isclose(probability_sum, 1.0, abs_tol=spec["acceptance"]["probability_sum_tolerance"]):
        failures.append("atom probabilities do not sum to one")
    with events_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != spec["event_schema"]:
            failures.append("event columns differ from the frozen schema")
        rows = list(reader)
    accepted = 0
    for line, row in enumerate(rows, start=2):
        try:
            atom = atoms[int(row["atom_id"])]
            if int(row["sign"]) != int(atom["sign"]):
                failures.append(f"line {line}: sign mismatch")
            if not math.isclose(float(row["node"]), float(atom["node"]), abs_tol=spec["acceptance"]["node_match_tolerance"]):
                failures.append(f"line {line}: node mismatch")
            accepted += row["accepted"].lower() in ("1", "true", "yes")
        except (KeyError, ValueError):
            failures.append(f"line {line}: invalid typed event")
    if accepted < int(spec["acceptance"]["minimum_accepted_events"]):
        failures.append("too few accepted events")
    return {
        "structurally_valid": not failures,
        "failures": failures,
        "rows": len(rows),
        "accepted": accepted,
        "physical_evidence_certified": False,
    }


def self_test(spec_path: Path) -> None:
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        events = root / "synthetic_structure_only.csv"
        with events.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=spec["event_schema"])
            writer.writeheader()
            for index in range(spec["acceptance"]["minimum_accepted_events"]):
                atom = spec["atoms"][index % len(spec["atoms"])]
                writer.writerow({
                    "event_id": index,
                    "timestamp_utc": "2000-01-01T00:00:00Z",
                    "run_id": "SYNTHETIC-STRUCTURE-TEST",
                    "blinded_config": "B",
                    "atom_id": atom["atom_id"],
                    "node": atom["node"],
                    "sign": atom["sign"],
                    "detector_channel": index % 12,
                    "accepted": True,
                })
        manifest = root / "manifest.json"
        manifest.write_text(json.dumps({
            "experiment_id": "SYNTHETIC-STRUCTURE-TEST-NOT-EVIDENCE",
            "apparatus_id": "SELF-TEST",
            "calibration_ids": ["SELF-TEST"],
            "provider": "role-a",
            "registrar": "role-b",
            "analyst": "role-c",
            "raw_events_sha256": digest(events),
            "holdout_frozen_utc": "2000-01-01T00:00:00Z",
            "unblinding_utc": "2000-01-02T00:00:00Z",
            "physical_units": "none: schema test",
            "pre_registered_tests": ["schema only"],
        }), encoding="utf-8")
        result = validate(spec_path, manifest, events)
        if not result["structurally_valid"]:
            raise SystemExit(json.dumps(result, indent=2))
        print(json.dumps(result | {"self_test_only": True}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--events", type=Path)
    parser.add_argument("--self-test", action="store_true")
    arguments = parser.parse_args()
    if arguments.self_test:
        self_test(arguments.spec)
        return
    if arguments.manifest is None or arguments.events is None:
        parser.error("--manifest and --events are required outside --self-test")
    result = validate(arguments.spec, arguments.manifest, arguments.events)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["structurally_valid"] else 1)


if __name__ == "__main__":
    main()
