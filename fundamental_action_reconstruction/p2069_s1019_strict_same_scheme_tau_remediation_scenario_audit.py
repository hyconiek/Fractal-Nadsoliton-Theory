#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2069_s1019_strict_same_scheme_tau_remediation_scenario_audit.json"
MD = GEN / "p2069_s1019_strict_same_scheme_tau_remediation_scenario_audit.md"

SCHEMA_VERSION = "p2069_s1019_v1"
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

    p2068 = load("p2068_s1018_strict_same_scheme_transportability_improvement_lever_audit.json")
    p2058 = load("p2058_s1008_strict_same_scheme_policy_regime_switch_audit.json")

    ready = p2068.get("result_kind") == "PASS_TRANSPORTABILITY_IMPROVEMENT_LEVER_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    comps = p2068.get("transportability_components") or {}
    radius_star = float(comps.get("radius_star", 0.0))

    tau_grid_raw = ((p2058.get("regime_grid") or {}).get("softmax_tau") or []) if isinstance(p2058, dict) else []
    tau_grid: list[float] = []
    for t in tau_grid_raw:
        try:
            tau_grid.append(abs(float(t)))
        except Exception:
            pass

    tau_grid = sorted(set(tau_grid))
    current_tau_amp = max(tau_grid, default=0.0)

    scenario_rows = []
    for tau in tau_grid:
        transportable = tau <= radius_star
        scenario_rows.append(
            {
                "tau_candidate": tau,
                "transportable_under_radius_star": transportable,
                "margin_to_radius_star": radius_star - tau,
            }
        )

    feasible = [r for r in scenario_rows if r["transportable_under_radius_star"]]
    minimal_feasible = max(feasible, key=lambda r: r["tau_candidate"]) if feasible else None

    if minimal_feasible is not None and current_tau_amp > 0:
        tau_needed = float(minimal_feasible["tau_candidate"])
        reduction_factor = tau_needed / current_tau_amp
    else:
        tau_needed = 0.0 if current_tau_amp > 0 else current_tau_amp
        reduction_factor = 0.0 if current_tau_amp > 0 else 1.0

    reduction_percent = max(0.0, min(100.0, (1.0 - reduction_factor) * 100.0))

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2069",
        "stage_id": "S1019",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_TAU_REMEDIATION_SCENARIO_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_TAU_REMEDIATION_SCENARIO_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2068_present": p2068.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2068_json_sha256": file_sha256(GEN / "p2068_s1018_strict_same_scheme_transportability_improvement_lever_audit.json"),
            "p2058_json_sha256": file_sha256(GEN / "p2058_s1008_strict_same_scheme_policy_regime_switch_audit.json"),
        },
        "tau_remediation": {
            "radius_star": radius_star,
            "current_tau_amplitude_proxy": current_tau_amp,
            "scenario_count": len(scenario_rows),
            "scenario_rows": scenario_rows,
            "minimal_feasible_tau_scenario": minimal_feasible,
            "required_tau_reduction_factor": reduction_factor,
            "required_tau_reduction_percent": reduction_percent,
        },
        "c3_gate_update": {
            "C3_tau_remediation_scenario_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "scenario_rows_nonempty": len(scenario_rows) > 0,
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
                "# P2069 S1019: tau remediation scenario audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- radius_star: `{radius_star}`",
                f"- current_tau_amplitude_proxy: `{current_tau_amp}`",
                f"- required_tau_reduction_factor: `{reduction_factor}`",
                f"- required_tau_reduction_percent: `{reduction_percent}`",
                "",
                "Discrete tau scenarios are audited against radius_star for strict transportability.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
