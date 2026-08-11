#!/usr/bin/env python3
"""ST189 external-record interface.

The validator can check an externally supplied JSONL event record and a prior
registration packet.  It deliberately cannot create custody, laboratory
execution, calibration, or an external signature.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


REQUIRED_EVENT_FIELDS = ("timestamp", "outcome", "config", "run_id")
REQUIRED_ROLES = ("provider", "registrar", "analyst")


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def digest_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def digest_json(value: Any) -> str:
    return digest_bytes(canonical_bytes(value))


def parse_jsonl_bytes(payload: bytes) -> list[dict]:
    rows = []
    for lineno, line in enumerate(payload.decode("utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        value = json.loads(line)
        if not isinstance(value, dict):
            raise ValueError(f"line {lineno} is not a JSON object")
        rows.append(value)
    return rows


def validate_payload(event_bytes: bytes, registration: dict, *, external_paths_supplied: bool) -> dict:
    try:
        events = parse_jsonl_bytes(event_bytes)
        parse_ok = True
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        events = []
        parse_ok = False
        parse_error = str(exc)
    else:
        parse_error = None

    event_fields_ok = bool(events) and all(all(field in row for field in REQUIRED_EVENT_FIELDS) for row in events)
    event_hash_ok = registration.get("event_file_sha256") == digest_bytes(event_bytes)
    holdout = registration.get("holdout", {})
    holdout_hash_ok = registration.get("holdout_sha256") == digest_json(holdout)
    roles = registration.get("roles", {})
    role_values = [roles.get(role) for role in REQUIRED_ROLES]
    roles_distinct = None not in role_values and len(set(role_values)) == len(REQUIRED_ROLES)
    registration_frozen = registration.get("frozen_before_unblinding") is True
    calibration_hash_present = isinstance(registration.get("calibration_sha256"), str) and len(registration["calibration_sha256"]) == 64
    custody_attested_externally = registration.get("independent_custody_attested") is True
    laboratory_execution_attested = registration.get("laboratory_execution_attested") is True

    structural = all((parse_ok, event_fields_ok, event_hash_ok, holdout_hash_ok, roles_distinct, registration_frozen))
    physical = structural and external_paths_supplied and calibration_hash_present and custody_attested_externally and laboratory_execution_attested
    checks = {
        "jsonl_parse": parse_ok,
        "nonempty_required_event_fields": event_fields_ok,
        "event_file_hash": event_hash_ok,
        "holdout_hash": holdout_hash_ok,
        "roles_declared_distinct": roles_distinct,
        "registration_frozen_before_unblinding": registration_frozen,
        "external_paths_supplied": external_paths_supplied,
        "calibration_hash_present": calibration_hash_present,
        "independent_custody_attested_externally": custody_attested_externally,
        "laboratory_execution_attested": laboratory_execution_attested,
    }
    return {
        "checks": checks,
        "event_count": len(events),
        "parse_error": parse_error,
        "structural_record_valid": structural,
        "physical_execution_valid": physical,
        "verdict": "PHYSICAL_EXECUTION_RECORD_VALID" if physical else (
            "EXTERNAL_RECORD_STRUCTURALLY_VALID_PHYSICAL_ATTESTATION_BLOCKED" if structural
            else "EXTERNAL_RECORD_INVALID"
        ),
        "boundary": (
            "This program validates syntax, hashes, a frozen holdout declaration, and distinct role labels. "
            "It cannot verify human identity, create independent custody, certify calibration, or prove that laboratory events occurred."
        ),
    }


def validate_files(event_path: Path, registration_path: Path) -> dict:
    event_bytes = event_path.read_bytes()
    registration = json.loads(registration_path.read_text(encoding="utf-8"))
    result = validate_payload(event_bytes, registration, external_paths_supplied=True)
    result["event_path"] = str(event_path)
    result["registration_path"] = str(registration_path)
    return result


def synthetic_self_test() -> dict:
    events = [
        {"timestamp": 0.0, "outcome": 1, "config": "synthetic_A", "run_id": "SELF_TEST"},
        {"timestamp": 0.1, "outcome": 0, "config": "synthetic_B", "run_id": "SELF_TEST"},
    ]
    payload = b"".join(canonical_bytes(row) + b"\n" for row in events)
    holdout = {"run_ids": ["SELF_TEST"], "rule": "synthetic validator self-test only"}
    registration = {
        "event_file_sha256": digest_bytes(payload),
        "holdout": holdout,
        "holdout_sha256": digest_json(holdout),
        "roles": {"provider": "SYNTHETIC_PROVIDER", "registrar": "SYNTHETIC_REGISTRAR", "analyst": "SYNTHETIC_ANALYST"},
        "frozen_before_unblinding": True,
        "calibration_sha256": None,
        "independent_custody_attested": False,
        "laboratory_execution_attested": False,
    }
    clean = validate_payload(payload, registration, external_paths_supplied=False)
    tampered = validate_payload(payload + b" ", registration, external_paths_supplied=False)
    return {
        "clean": clean,
        "tampered": tampered,
        "synthetic_only": True,
        "no_external_record_supplied": True,
        "test_passed": clean["structural_record_valid"] and not clean["physical_execution_valid"] and not tampered["structural_record_valid"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("events", nargs="?", type=Path)
    parser.add_argument("registration", nargs="?", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        result = synthetic_self_test()
    elif args.events is not None and args.registration is not None:
        result = validate_files(args.events, args.registration)
    else:
        parser.error("provide EVENTS REGISTRATION or --self-test")
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False))


if __name__ == "__main__":
    main()
