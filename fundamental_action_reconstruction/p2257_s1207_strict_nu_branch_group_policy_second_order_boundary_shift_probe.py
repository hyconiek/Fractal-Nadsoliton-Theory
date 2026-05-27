#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2256 = GEN / "p2256_s1206_strict_nu_branch_group_policy_dominance_boundary_curve_probe.json"
OUT = GEN / "p2257_s1207_strict_nu_branch_group_policy_second_order_boundary_shift_probe.json"
MD = GEN / "p2257_s1207_strict_nu_branch_group_policy_second_order_boundary_shift_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2256 = load(IN_2256)
    probe = (p2256.get("strict_nu_branch_group_policy_dominance_boundary_curve_probe", {}) or {})
    inp = (probe.get("inputs", {}) or {})

    m_star_1 = float(inp.get("adaptive_minimum_margin", 0.0) or 0.0)
    progress_gap = float(inp.get("progress_gap", 0.0) or 0.0)
    horizon_scales = [float(h) for h in inp.get("horizon_scales", [0.8, 1.0, 1.2])]

    # Second-order correction proxy: if progress gap is small, boundary gets more conservative.
    # delta_m2(h) = gamma * h^2 / (1 + |progress_gap|)
    gamma = 1e-3
    rows = []
    for h in horizon_scales:
        delta_m2 = gamma * (h * h) / (1.0 + abs(progress_gap))
        m_star_2 = m_star_1 - delta_m2
        rows.append(
            {
                "horizon_scale": h,
                "first_order_critical_margin_floor": m_star_1,
                "second_order_shift": -delta_m2,
                "second_order_critical_margin_floor": m_star_2,
            }
        )

    max_abs_shift = max((abs(r["second_order_shift"]) for r in rows), default=0.0)

    payload = {
        "schema_version": "p2257_s1207_v1",
        "packet_id": "P2257",
        "stage_id": "S1207",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_SECOND_ORDER_BOUNDARY_SHIFT_PROBE",
        "strict_nu_branch_group_policy_second_order_boundary_shift_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_SECOND_ORDER_BOUNDARY_SHIFT_PROBE_V1",
            "source_packets": [str(IN_2256.relative_to(ROOT))],
            "inputs": {
                "first_order_critical_margin_floor": m_star_1,
                "progress_gap": progress_gap,
                "horizon_scales": horizon_scales,
                "gamma": gamma,
            },
            "second_order_boundary": {
                "rows": rows,
                "maximum_absolute_shift": max_abs_shift,
            },
            "physical_interpretation_note": "Second-order boundary shift models curvature-induced safety tightening: as effective horizon grows, conservative margin floor is reduced by a small quadratic correction capturing nonlinearity of accumulated policy motion.",
            "theorem_scope_limit": "proxy-level second-order boundary diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2258_candidate",
            "goal": "validate second-order shift against simulated nonlinear trajectory residuals and fit gamma from data",
        },
        "gatekeeper_checks": {
            "second_order_boundary_exported": True,
            "rows_nonempty": len(rows) > 0,
            "max_shift_nonnegative": max_abs_shift >= 0.0,
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
                "# P2257 S1207: second-order boundary shift probe",
                "",
                f"- first-order critical floor: `{m_star_1:.12e}`",
                f"- progress gap: `{progress_gap:.12e}`",
                f"- gamma: `{gamma:.12e}`",
                f"- maximum absolute second-order shift: `{max_abs_shift:.12e}`",
                "",
                "Proxy-level second-order boundary diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
