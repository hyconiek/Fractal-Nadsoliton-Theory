#!/usr/bin/env python3
"""P2057 S1007: policy-regret audit.

Quantify regret relative to oracle-best policy at same compute budget, with
bootstrap CI for regret. Audit-only; C3 remains OPEN.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2057_s1007_strict_same_scheme_policy_regret_audit.json"
MD = GEN / "p2057_s1007_strict_same_scheme_policy_regret_audit.md"

SCHEMA_VERSION = "p2057_s1007_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2056 = load("p2056_s1006_strict_same_scheme_allocation_policy_calibration_audit.json")

    ready = p2056.get("result_kind") == "PASS_ALLOCATION_POLICY_CALIBRATION_AUDIT_WITH_BOOTSTRAP_GAIN_VARIANCE__C3_STILL_OPEN"
    rows = p2056.get("policy_results") or []

    # Oracle by max bootstrap mean gain.
    oracle = max(rows, key=lambda r: float(r.get("mean", 0.0))) if rows else None
    selected = p2056.get("selected_policy") or {}

    selected_mean = float(selected.get("mean", 0.0))
    oracle_mean = float(oracle.get("mean", 0.0)) if oracle else 0.0
    regret_mean = oracle_mean - selected_mean

    selected_ci = (selected.get("ci95") or {"low": 0.0, "high": 0.0})
    oracle_ci = (oracle.get("ci95") if oracle else {"low": 0.0, "high": 0.0})

    # Conservative interval arithmetic for regret = oracle - selected
    regret_ci = {
        "low": float(oracle_ci.get("low", 0.0)) - float(selected_ci.get("high", 0.0)),
        "high": float(oracle_ci.get("high", 0.0)) - float(selected_ci.get("low", 0.0)),
    }

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2057",
        "stage_id": "S1007",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_POLICY_REGRET_AUDIT_WITH_BOOTSTRAP_CI__C3_STILL_OPEN"
            if ready and rows and oracle
            else "OPEN_POLICY_REGRET_AUDIT_BLOCKED"
        ),
        "depends_on": {"p2056_present": p2056.get("_missing") is None, "preconditions_ready": ready},
        "input_hashes": {
            "p2056_json_sha256": file_sha256(GEN / "p2056_s1006_strict_same_scheme_allocation_policy_calibration_audit.json"),
        },
        "regret_definition": "oracle_mean_gain - selected_policy_mean_gain",
        "oracle_policy": oracle,
        "selected_policy": selected,
        "regret_summary": {
            "oracle_mean_gain": oracle_mean,
            "selected_mean_gain": selected_mean,
            "regret_mean": regret_mean,
            "regret_ci95": regret_ci,
        },
        "c3_gate_update": {
            "C3_policy_regret_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "theorem-grade operator identity across full background family",
                "cross-background finite-part transport identity theorem",
                "global finite-part lock/cocycle theorem",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "policies_nonempty": len(rows) > 0,
            "oracle_present": oracle is not None,
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2057 S1007: policy-regret audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Regret mean: `{regret_mean}`",
                f"- Regret CI95: `[{regret_ci['low']}, {regret_ci['high']}]`",
                "",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
