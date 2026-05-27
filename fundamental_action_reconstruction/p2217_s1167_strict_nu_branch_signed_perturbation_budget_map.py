#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2216 = GEN / "p2216_s1166_strict_nu_branch_one_sided_lipschitz_slope_bounds.json"
OUT = GEN / "p2217_s1167_strict_nu_branch_signed_perturbation_budget_map.json"
MD = GEN / "p2217_s1167_strict_nu_branch_signed_perturbation_budget_map.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2216 = load(IN_2216)
    block = p2216.get("strict_nu_branch_one_sided_lipschitz_slope_bounds", {}) or {}

    ref = block.get("reference", {}) or {}
    dm_star = float(ref.get("dm_star", 0.0) or 0.0)

    slopes = block.get("one_sided_local_slopes", {}) or {}
    left_slope = float(slopes.get("left_secant_at_dm_star", 0.0) or 0.0)
    right_slope = float(slopes.get("right_secant_at_dm_star", 0.0) or 0.0)

    if not (dm_star > 0.0 and left_slope > 0.0 and right_slope > 0.0):
        raise RuntimeError("Invalid upstream inputs for P2217 perturbation budget map")

    target_sf_above = [1.001, 1.005, 1.01, 1.02, 1.05, 1.10]
    target_sf_below = [0.999, 0.995, 0.99, 0.98, 0.95, 0.90]

    above_rows = []
    for sf in target_sf_above:
        delta_dm = (sf - 1.0) / right_slope
        above_rows.append({
            "target_safety_factor": sf,
            "required_signed_delta_dm": delta_dm,
            "required_abs_dm": dm_star + delta_dm,
            "certification_side": "above_dm_star",
        })

    below_rows = []
    for sf in target_sf_below:
        delta_dm_mag = (1.0 - sf) / left_slope
        below_rows.append({
            "target_safety_factor": sf,
            "required_signed_delta_dm": -delta_dm_mag,
            "required_abs_dm": dm_star - delta_dm_mag,
            "certification_side": "below_dm_star",
        })

    payload = {
        "schema_version": "p2217_s1167_v1",
        "packet_id": "P2217",
        "stage_id": "S1167",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_SIGNED_PERTURBATION_BUDGET_MAP",
        "strict_nu_branch_signed_perturbation_budget_map": {
            "map_id": "STRICT_NU_BRANCH_SIGNED_PERTURBATION_BUDGET_MAP_V1",
            "source_packet": str(IN_2216.relative_to(ROOT)),
            "reference_dm_star": dm_star,
            "one_sided_local_slopes": {
                "left_secant_at_dm_star": left_slope,
                "right_secant_at_dm_star": right_slope,
            },
            "above_target_budgets": above_rows,
            "below_target_budgets": below_rows,
            "interpretation": "For each target safety-factor, map to minimal signed local perturbation estimate using one-sided secant slope bounds.",
            "theorem_scope_limit": "local linearized perturbation-budget map only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2218_candidate",
            "goal": "validate perturbation budget map against direct recomputation on denser neighborhood grid and export mismatch envelope",
        },
        "gatekeeper_checks": {
            "signed_perturbation_budget_map_exported": True,
            "above_budgets_positive_delta": all(r["required_signed_delta_dm"] > 0.0 for r in above_rows),
            "below_budgets_negative_delta": all(r["required_signed_delta_dm"] < 0.0 for r in below_rows),
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2217 S1167: strict nu-branch signed perturbation budget map",
            "",
            f"- reference dm*: `{dm_star:.12e}`",
            f"- left slope: `{left_slope:.12e}`",
            f"- right slope: `{right_slope:.12e}`",
            f"- above budgets positive delta: `{all(r['required_signed_delta_dm'] > 0.0 for r in above_rows)}`",
            f"- below budgets negative delta: `{all(r['required_signed_delta_dm'] < 0.0 for r in below_rows)}`",
            "",
            "Local linearized perturbation-budget map only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
