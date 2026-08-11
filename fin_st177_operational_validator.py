#!/usr/bin/env python3
"""Executable ST177 validator for the conditional FIN operational tuple."""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "FIN_ST177_Operational_Validator_Record.json"


ATOMS = ["ALG", "GAUGE", "PREP", "DYN", "CLOCK", "INST", "ENV", "APP", "SELECT", "UNITS", "RECORD"]


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def template_packet() -> dict:
    raw_events: list[dict] = []
    holdout = {"rule": "single frozen analysis after registration", "indices": [], "status": "empty_no_laboratory_data"}
    return {
        "A_strict": {"source": "strict decimal kernel", "dimension": 12, "dimensionless": True},
        "ALG": {"algebra": "typed finite observable algebra", "composition_embeddings": "declared"},
        "GAUGE": {"equivalence": "declared operational carrier equivalence"},
        "PREP": {"preparations": ["vertex_basis", "informationally_complete_placeholder"], "realized": False},
        "DYN": {"channels": ["unitary", "heat"], "selected_channel": None},
        "CLOCK": {"parameter": "tau", "calibration_to_seconds": None},
        "INST": {"instrument": "vertex POVM placeholder", "tomography": None},
        "ENV": {"bath_model": None, "dilation": None},
        "APP": {"apparatus": None, "calibration_record": None},
        "SELECT": {"branch_rule": None, "strict_source": False},
        "UNITS": {"mode": "dimensionless_only", "energy_scale": None, "time_scale": None, "length_scale": None},
        "RECORD": {"raw_events": raw_events, "raw_events_sha256": digest(raw_events),
                   "holdout": holdout, "holdout_sha256": digest(holdout)},
        "roles": {"provider": "UNASSIGNED_PROVIDER", "registrar": "UNASSIGNED_REGISTRAR", "analyst": "UNASSIGNED_ANALYST"},
        "independent_custody_verified": False,
        "laboratory_execution_verified": False,
    }


def validate(packet: dict) -> dict:
    checks: dict[str, bool] = {}
    for atom in ATOMS:
        checks[f"atom_{atom}"] = atom in packet and packet[atom] is not None
    record = packet.get("RECORD") or {}
    checks["record_hash"] = record.get("raw_events_sha256") == digest(record.get("raw_events", []))
    checks["holdout_hash"] = record.get("holdout_sha256") == digest(record.get("holdout", {}))
    roles = packet.get("roles", {})
    role_values = [roles.get(k) for k in ["provider", "registrar", "analyst"]]
    checks["roles_declared_distinct"] = None not in role_values and len(set(role_values)) == 3
    checks["physical_units_supplied"] = packet.get("UNITS", {}).get("mode") == "physical_calibrated"
    checks["clock_calibrated"] = packet.get("CLOCK", {}).get("calibration_to_seconds") is not None
    checks["apparatus_calibrated"] = packet.get("APP", {}).get("calibration_record") is not None
    checks["independent_custody"] = packet.get("independent_custody_verified") is True
    checks["laboratory_execution"] = packet.get("laboratory_execution_verified") is True
    structural = all(checks[f"atom_{a}"] for a in ATOMS) and checks["record_hash"] and checks["holdout_hash"] and checks["roles_declared_distinct"]
    physical = structural and all(checks[k] for k in ["physical_units_supplied", "clock_calibrated", "apparatus_calibrated", "independent_custody", "laboratory_execution"])
    return {"checks": checks, "structural_packet_valid": structural, "physical_execution_valid": physical,
            "verdict": "PHYSICAL_EXECUTION_VALID" if physical else ("MATHEMATICAL_PACKET_READY_PHYSICAL_EXECUTION_BLOCKED" if structural else "STRUCTURAL_PACKET_INVALID")}


def deletion_audit(packet: dict) -> list[dict]:
    rows = []
    for atom in ATOMS:
        trial = deepcopy(packet); trial.pop(atom, None)
        result = validate(trial)
        rows.append({"deleted_atom": atom, "corresponding_check": f"atom_{atom}",
                     "corresponding_check_passed": result["checks"][f"atom_{atom}"],
                     "structural_packet_valid": result["structural_packet_valid"]})
    return rows


def run(write: bool = True) -> dict:
    packet = template_packet()
    validation = validate(packet)
    result = {
        "packet": packet,
        "packet_sha256": digest(packet),
        "validation": validation,
        "deletion_audit": deletion_audit(packet),
        "all_eleven_deletions_detected": all(not r["corresponding_check_passed"] for r in deletion_audit(packet)),
        "boundary": "Code verifies schema, hashes, declared role separation, and deletion behavior. It cannot create physical units, apparatus calibration, laboratory execution, or independent custody.",
    }
    if write:
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(run(True), indent=2, sort_keys=True, ensure_ascii=False))
