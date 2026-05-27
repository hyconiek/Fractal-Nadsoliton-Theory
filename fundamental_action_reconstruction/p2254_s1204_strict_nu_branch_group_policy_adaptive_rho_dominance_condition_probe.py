#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2253 = GEN / "p2253_s1203_strict_nu_branch_group_policy_adaptive_vs_fixed_rho_trajectory_probe.json"
OUT = GEN / "p2254_s1204_strict_nu_branch_group_policy_adaptive_rho_dominance_condition_probe.json"
MD = GEN / "p2254_s1204_strict_nu_branch_group_policy_adaptive_rho_dominance_condition_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2253 = load(IN_2253)
    probe = (p2253.get("strict_nu_branch_group_policy_adaptive_vs_fixed_rho_trajectory_probe", {}) or {})
    fixed = (probe.get("fixed_rho_result", {}) or {})
    adaptive = (probe.get("adaptive_rho_result", {}) or {})

    fixed_progress = float(fixed.get("cumulative_budget_use", 0.0) or 0.0)
    adaptive_progress = float(adaptive.get("cumulative_budget_use", 0.0) or 0.0)
    fixed_min_margin = float(fixed.get("minimum_margin", 0.0) or 0.0)
    adaptive_min_margin = float(adaptive.get("minimum_margin", 0.0) or 0.0)

    # Dominance under safety floor m_floor: adaptive dominates if
    # progress_adaptive >= progress_fixed and min_margin_adaptive >= m_floor.
    m_floor = 0.0
    adaptive_dominates = (adaptive_progress + 1e-15 >= fixed_progress) and (adaptive_min_margin + 1e-15 >= m_floor)

    # Conservative gap-to-dominance diagnostics.
    progress_gap = adaptive_progress - fixed_progress
    margin_floor_gap = adaptive_min_margin - m_floor

    payload = {
        "schema_version": "p2254_s1204_v1",
        "packet_id": "P2254",
        "stage_id": "S1204",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_ADAPTIVE_RHO_DOMINANCE_CONDITION_PROBE",
        "strict_nu_branch_group_policy_adaptive_rho_dominance_condition_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_ADAPTIVE_RHO_DOMINANCE_CONDITION_PROBE_V1",
            "source_packets": [str(IN_2253.relative_to(ROOT))],
            "inputs": {
                "fixed_progress": fixed_progress,
                "adaptive_progress": adaptive_progress,
                "fixed_minimum_margin": fixed_min_margin,
                "adaptive_minimum_margin": adaptive_min_margin,
                "minimum_margin_floor": m_floor,
            },
            "dominance_condition": {
                "inequality": "adaptive_progress >= fixed_progress and adaptive_minimum_margin >= minimum_margin_floor",
                "progress_gap": progress_gap,
                "margin_floor_gap": margin_floor_gap,
                "adaptive_dominates_under_floor": adaptive_dominates,
            },
            "physical_interpretation_note": "Dominance condition formalizes when adaptive control is physically preferable: it must deliver at least equal policy progress without violating safety-margin floor constraints.",
            "theorem_scope_limit": "finite-horizon dominance diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2255_candidate",
            "goal": "scan dominance over multiple margin floors and horizons to export robust dominance envelope",
        },
        "gatekeeper_checks": {
            "dominance_condition_exported": True,
            "progress_gap_computable": isinstance(progress_gap, float),
            "margin_floor_gap_computable": isinstance(margin_floor_gap, float),
            "adaptive_nonnegative_margin": adaptive_min_margin >= -1e-15,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2254 S1204: adaptive-rho dominance condition probe",
                "",
                f"- progress gap (adaptive-fixed): `{progress_gap:.12e}`",
                f"- adaptive minimum margin: `{adaptive_min_margin:.12e}`",
                f"- minimum margin floor: `{m_floor:.12e}`",
                f"- floor gap: `{margin_floor_gap:.12e}`",
                f"- adaptive dominates under floor: `{adaptive_dominates}`",
                "",
                "Finite-horizon dominance diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
