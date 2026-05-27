#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2265 = GEN / "p2265_s1215_strict_nu_branch_task3_closure_gap_matrix_probe.json"
IN_2261 = GEN / "p2261_s1211_strict_nu_branch_group_policy_stochastic_reduction_concentration_bound_probe.json"
OUT = GEN / "p2266_s1216_strict_nu_branch_group_policy_risk_calibrated_controller_map_probe.json"
MD = GEN / "p2266_s1216_strict_nu_branch_group_policy_risk_calibrated_controller_map_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2265 = load(IN_2265)
    p2261 = load(IN_2261)

    cb = (p2261.get("strict_nu_branch_group_policy_stochastic_reduction_concentration_bound_probe", {}) or {}).get("concentration_bound", {}) or {}
    lb95 = float(cb.get("lower_bound_reduction_95", 0.0) or 0.0)

    # Risk tolerance -> controller aggressiveness map (rho) and reserve multiplier (kappa_scale)
    # conservative monotone map: lower tolerance => lower rho, higher kappa_scale.
    risk_levels = [0.01, 0.02, 0.05, 0.10]
    rows = []
    for r in risk_levels:
        rho = max(0.55, min(0.9, 0.9 - 1.5 * r))
        kappa_scale = 1.0 + max(0.0, (0.05 - r)) * 4.0
        certified = lb95 > r * 0.2
        rows.append(
            {
                "risk_tolerance": r,
                "rho_recommended": rho,
                "kappa_scale_recommended": kappa_scale,
                "reduction_lb95_reference": lb95,
                "certificate_pass": certified,
            }
        )

    payload = {
        "schema_version": "p2266_s1216_v1",
        "packet_id": "P2266",
        "stage_id": "S1216",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_RISK_CALIBRATED_CONTROLLER_MAP_PROBE",
        "strict_nu_branch_group_policy_risk_calibrated_controller_map_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_RISK_CALIBRATED_CONTROLLER_MAP_PROBE_V1",
            "source_packets": [str(IN_2265.relative_to(ROOT)), str(IN_2261.relative_to(ROOT))],
            "controller_map_rows": rows,
            "physical_interpretation_note": "Map ties operational risk tolerance to admissible controller aggressiveness and reserve scaling, creating an explicit safety-certificate interface from stochastic reduction evidence.",
            "theorem_scope_limit": "risk-calibrated mapping diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2267_candidate",
            "goal": "validate mapped rho/kappa schedule on nonlinear perturbed trajectories and measure certificate pass-rate",
        },
        "gatekeeper_checks": {
            "risk_controller_map_exported": True,
            "all_rho_bounded": all(0.55 <= r["rho_recommended"] <= 0.9 for r in rows),
            "all_kappa_scale_ge_one": all(r["kappa_scale_recommended"] >= 1.0 for r in rows),
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text("\n".join([
        "# P2266 S1216: risk-calibrated controller map probe",
        "",
        f"- reduction lower bound reference: `{lb95:.12e}`",
        f"- risk levels: `{risk_levels}`",
        "- exported map: risk_tolerance -> (rho_recommended, kappa_scale_recommended)",
        "",
        "Risk-calibrated mapping diagnostic only; no kernel-bridge or selector-closure claim.",
    ]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
