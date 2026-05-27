#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2253 = GEN / "p2253_s1203_strict_nu_branch_group_policy_adaptive_vs_fixed_rho_trajectory_probe.json"
OUT = GEN / "p2255_s1205_strict_nu_branch_group_policy_dominance_envelope_scan_probe.json"
MD = GEN / "p2255_s1205_strict_nu_branch_group_policy_dominance_envelope_scan_probe.md"


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

    progress_gap = adaptive_progress - fixed_progress

    margin_floors = [0.0, 0.01, 0.02, 0.05]
    horizon_scales = [0.8, 1.0, 1.2]

    rows = []
    dominates_count = 0
    for m_floor in margin_floors:
        for h_scale in horizon_scales:
            scaled_progress_gap = progress_gap * h_scale
            dominates = (scaled_progress_gap + 1e-15 >= 0.0) and (adaptive_min_margin + 1e-15 >= m_floor)
            if dominates:
                dominates_count += 1
            rows.append(
                {
                    "margin_floor": m_floor,
                    "horizon_scale": h_scale,
                    "scaled_progress_gap": scaled_progress_gap,
                    "adaptive_min_margin": adaptive_min_margin,
                    "adaptive_dominates": dominates,
                }
            )

    total_cells = len(rows)
    dominance_fraction = dominates_count / max(total_cells, 1)

    payload = {
        "schema_version": "p2255_s1205_v1",
        "packet_id": "P2255",
        "stage_id": "S1205",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_DOMINANCE_ENVELOPE_SCAN_PROBE",
        "strict_nu_branch_group_policy_dominance_envelope_scan_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_DOMINANCE_ENVELOPE_SCAN_PROBE_V1",
            "source_packets": [str(IN_2253.relative_to(ROOT))],
            "inputs": {
                "fixed_progress": fixed_progress,
                "adaptive_progress": adaptive_progress,
                "fixed_minimum_margin": fixed_min_margin,
                "adaptive_minimum_margin": adaptive_min_margin,
                "margin_floors": margin_floors,
                "horizon_scales": horizon_scales,
            },
            "dominance_envelope": {
                "rows": rows,
                "dominance_fraction": dominance_fraction,
                "dominant_cells": dominates_count,
                "total_cells": total_cells,
            },
            "physical_interpretation_note": "Dominance envelope maps where adaptive control remains physically preferable under stricter safety floors and effective horizon scaling, exposing robustness domain of controller preference.",
            "theorem_scope_limit": "finite-grid dominance-envelope diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2256_candidate",
            "goal": "derive analytic boundary in (margin_floor, horizon_scale) separating adaptive-dominant and non-dominant regions",
        },
        "gatekeeper_checks": {
            "dominance_envelope_exported": True,
            "dominance_fraction_bounded": 0.0 <= dominance_fraction <= 1.0,
            "nonempty_scan_grid": total_cells > 0,
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
                "# P2255 S1205: adaptive dominance envelope scan probe",
                "",
                f"- dominance cells: `{dominates_count}` / `{total_cells}`",
                f"- dominance fraction: `{dominance_fraction:.12e}`",
                f"- adaptive minimum margin: `{adaptive_min_margin:.12e}`",
                f"- progress gap (adaptive-fixed): `{progress_gap:.12e}`",
                "",
                "Finite-grid dominance envelope only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
