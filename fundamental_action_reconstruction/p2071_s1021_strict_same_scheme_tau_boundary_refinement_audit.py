#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2071_s1021_strict_same_scheme_tau_boundary_refinement_audit.json"
MD = GEN / "p2071_s1021_strict_same_scheme_tau_boundary_refinement_audit.md"

SCHEMA_VERSION = "p2071_s1021_v1"
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
    p2070 = load("p2070_s1020_strict_same_scheme_tau_grid_extension_feasibility_audit.json")
    ready = p2070.get("result_kind") == "PASS_TAU_GRID_EXTENSION_FEASIBILITY_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    block = p2070.get("tau_grid_extension_audit") or {}
    radius_star = float(block.get("radius_star", 0.0))
    current_tau = float(block.get("current_tau_amplitude_proxy", 0.0))

    # local refinement around boundary tau ~ radius_star with conservative ±0.005 window
    refine_grid = sorted(set([round(radius_star + d, 6) for d in (-0.005, -0.0025, -0.001, 0.0, 0.001, 0.0025, 0.005) if radius_star + d > 0]))

    rows = []
    for tau in refine_grid:
        rows.append({
            "tau_candidate": tau,
            "transportable_under_radius_star": tau <= radius_star,
            "signed_margin": radius_star - tau,
        })

    feasible = [r for r in rows if r["transportable_under_radius_star"]]
    tight_feasible = max(feasible, key=lambda r: r["tau_candidate"]) if feasible else None

    if tight_feasible is not None and current_tau > 0:
        req_factor = tight_feasible["tau_candidate"] / current_tau
    else:
        req_factor = 0.0 if current_tau > 0 else 1.0
    req_percent = max(0.0, min(100.0, (1.0 - req_factor) * 100.0))

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2071",
        "stage_id": "S1021",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_TAU_BOUNDARY_REFINEMENT_AUDIT_WITH_TRACE__C3_STILL_OPEN" if ready else "OPEN_TAU_BOUNDARY_REFINEMENT_AUDIT_BLOCKED",
        "depends_on": {"p2070_present": p2070.get("_missing") is None, "preconditions_ready": ready},
        "input_hashes": {
            "p2070_json_sha256": file_sha256(GEN / "p2070_s1020_strict_same_scheme_tau_grid_extension_feasibility_audit.json"),
        },
        "tau_boundary_refinement": {
            "radius_star": radius_star,
            "current_tau_amplitude_proxy": current_tau,
            "refine_grid": refine_grid,
            "rows": rows,
            "tight_feasible_tau": tight_feasible,
            "required_tau_reduction_factor": req_factor,
            "required_tau_reduction_percent": req_percent,
        },
        "c3_gate_update": {
            "C3_tau_boundary_refinement_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "rows_nonempty": len(rows) > 0,
            "tight_feasible_exists": tight_feasible is not None,
            "required_tau_reduction_factor_in_0_1": 0.0 <= req_factor <= 1.0,
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2071 S1021: tau boundary refinement audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- radius_star: `{radius_star}`",
            f"- tight_feasible_tau: `{tight_feasible}`",
            f"- required_tau_reduction_factor: `{req_factor}`",
            f"- required_tau_reduction_percent: `{req_percent}`",
            "",
            "This stage refines tau candidates around the transportability boundary.",
            "C3 remains OPEN (not discharged).",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
