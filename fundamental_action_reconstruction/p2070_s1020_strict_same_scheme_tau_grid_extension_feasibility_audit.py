#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2070_s1020_strict_same_scheme_tau_grid_extension_feasibility_audit.json"
MD = GEN / "p2070_s1020_strict_same_scheme_tau_grid_extension_feasibility_audit.md"

SCHEMA_VERSION = "p2070_s1020_v1"
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

    p2069 = load("p2069_s1019_strict_same_scheme_tau_remediation_scenario_audit.json")
    ready = p2069.get("result_kind") == "PASS_TAU_REMEDIATION_SCENARIO_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    tau_block = p2069.get("tau_remediation") or {}
    radius_star = float(tau_block.get("radius_star", 0.0))
    current_tau_amp = float(tau_block.get("current_tau_amplitude_proxy", 0.0))

    # Honest extension below previous lower bound 0.05, focused around radius_star.
    extension_grid = sorted(set([0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04]))
    extension_rows = []
    for tau in extension_grid:
        feasible = tau <= radius_star
        extension_rows.append(
            {
                "tau_candidate": tau,
                "transportable_under_radius_star": feasible,
                "margin_to_radius_star": radius_star - tau,
            }
        )

    feasible_rows = [r for r in extension_rows if r["transportable_under_radius_star"]]
    maximal_feasible_tau = max(feasible_rows, key=lambda r: r["tau_candidate"]) if feasible_rows else None

    if maximal_feasible_tau is not None and current_tau_amp > 0.0:
        tau_target = float(maximal_feasible_tau["tau_candidate"])
        reduction_factor = tau_target / current_tau_amp
    else:
        tau_target = 0.0 if current_tau_amp > 0.0 else current_tau_amp
        reduction_factor = 0.0 if current_tau_amp > 0.0 else 1.0

    reduction_percent = max(0.0, min(100.0, (1.0 - reduction_factor) * 100.0))

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2070",
        "stage_id": "S1020",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_TAU_GRID_EXTENSION_FEASIBILITY_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_TAU_GRID_EXTENSION_FEASIBILITY_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2069_present": p2069.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2069_json_sha256": file_sha256(GEN / "p2069_s1019_strict_same_scheme_tau_remediation_scenario_audit.json"),
        },
        "tau_grid_extension_audit": {
            "radius_star": radius_star,
            "current_tau_amplitude_proxy": current_tau_amp,
            "extension_grid": extension_grid,
            "extension_rows": extension_rows,
            "maximal_feasible_tau_in_extension": maximal_feasible_tau,
            "required_tau_reduction_factor": reduction_factor,
            "required_tau_reduction_percent": reduction_percent,
        },
        "c3_gate_update": {
            "C3_tau_grid_extension_feasibility_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "extension_rows_nonempty": len(extension_rows) > 0,
            "maximal_feasible_exists": maximal_feasible_tau is not None,
            "required_tau_reduction_factor_in_0_1": 0.0 <= reduction_factor <= 1.0,
            "required_tau_reduction_percent_nonnegative": reduction_percent >= 0.0,
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
                "# P2070 S1020: tau-grid extension feasibility audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- radius_star: `{radius_star}`",
                f"- maximal_feasible_tau_in_extension: `{maximal_feasible_tau}`",
                f"- required_tau_reduction_factor: `{reduction_factor}`",
                f"- required_tau_reduction_percent: `{reduction_percent}`",
                "",
                "This stage extends tau below prior grid minimum and locates first feasible strict scenario.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
