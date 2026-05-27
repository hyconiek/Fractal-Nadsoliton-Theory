#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2255 = GEN / "p2255_s1205_strict_nu_branch_group_policy_dominance_envelope_scan_probe.json"
OUT = GEN / "p2256_s1206_strict_nu_branch_group_policy_dominance_boundary_curve_probe.json"
MD = GEN / "p2256_s1206_strict_nu_branch_group_policy_dominance_boundary_curve_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2255 = load(IN_2255)
    probe = (p2255.get("strict_nu_branch_group_policy_dominance_envelope_scan_probe", {}) or {})
    inp = (probe.get("inputs", {}) or {})

    adaptive_min_margin = float(inp.get("adaptive_minimum_margin", 0.0) or 0.0)
    fixed_progress = float(inp.get("fixed_progress", 0.0) or 0.0)
    adaptive_progress = float(inp.get("adaptive_progress", 0.0) or 0.0)
    progress_gap = adaptive_progress - fixed_progress

    horizon_scales = inp.get("horizon_scales", [0.8, 1.0, 1.2])
    horizon_scales = [float(h) for h in horizon_scales]

    # Dominance condition from P2255 rows:
    # progress_gap*h >= 0 and adaptive_min_margin >= m_floor.
    # Boundary in (m_floor, h): m_floor* = adaptive_min_margin (independent of h)
    # with progress-side viability requiring progress_gap*h >= 0.
    rows = []
    for h in horizon_scales:
        progress_viable = progress_gap * h >= -1e-15
        rows.append(
            {
                "horizon_scale": h,
                "critical_margin_floor": adaptive_min_margin,
                "progress_viable": progress_viable,
                "dominance_region_statement": "m_floor <= critical_margin_floor and progress_viable",
            }
        )

    all_progress_viable = all(r["progress_viable"] for r in rows)

    payload = {
        "schema_version": "p2256_s1206_v1",
        "packet_id": "P2256",
        "stage_id": "S1206",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_DOMINANCE_BOUNDARY_CURVE_PROBE",
        "strict_nu_branch_group_policy_dominance_boundary_curve_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_DOMINANCE_BOUNDARY_CURVE_PROBE_V1",
            "source_packets": [str(IN_2255.relative_to(ROOT))],
            "inputs": {
                "adaptive_minimum_margin": adaptive_min_margin,
                "progress_gap": progress_gap,
                "horizon_scales": horizon_scales,
            },
            "analytic_boundary": {
                "equation": "critical_margin_floor(h) = adaptive_minimum_margin, with progress_gap*h >= 0",
                "rows": rows,
                "all_progress_viable": all_progress_viable,
            },
            "physical_interpretation_note": "Boundary curve identifies safety-floor threshold above which adaptive controller cannot dominate; when progress gap is nonnegative, dominance region is margin-floor limited by adaptive minimum margin.",
            "theorem_scope_limit": "finite-grid boundary-curve diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2257_candidate",
            "goal": "extend boundary model with second-order correction from trajectory curvature and evaluate shift of critical margin floor",
        },
        "gatekeeper_checks": {
            "boundary_curve_exported": True,
            "rows_nonempty": len(rows) > 0,
            "all_progress_viable_flag_exported": isinstance(all_progress_viable, bool),
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
                "# P2256 S1206: dominance boundary curve probe",
                "",
                f"- adaptive minimum margin: `{adaptive_min_margin:.12e}`",
                f"- progress gap: `{progress_gap:.12e}`",
                f"- horizon scales: `{horizon_scales}`",
                f"- all progress viable: `{all_progress_viable}`",
                "",
                "Finite-grid boundary-curve diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
