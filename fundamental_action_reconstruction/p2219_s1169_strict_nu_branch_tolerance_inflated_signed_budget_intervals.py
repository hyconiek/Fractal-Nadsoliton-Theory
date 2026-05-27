#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2217 = GEN / "p2217_s1167_strict_nu_branch_signed_perturbation_budget_map.json"
IN_2218 = GEN / "p2218_s1168_strict_nu_branch_budget_map_direct_recompute_mismatch_envelope.json"
OUT = GEN / "p2219_s1169_strict_nu_branch_tolerance_inflated_signed_budget_intervals.json"
MD = GEN / "p2219_s1169_strict_nu_branch_tolerance_inflated_signed_budget_intervals.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2217 = load(IN_2217)
    p2218 = load(IN_2218)

    bmap = p2217.get("strict_nu_branch_signed_perturbation_budget_map", {}) or {}
    dm_star = float(bmap.get("reference_dm_star", 0.0) or 0.0)
    slopes = bmap.get("one_sided_local_slopes", {}) or {}
    left = float(slopes.get("left_secant_at_dm_star", 0.0) or 0.0)
    right = float(slopes.get("right_secant_at_dm_star", 0.0) or 0.0)

    env = (p2218.get("strict_nu_branch_budget_map_direct_recompute_mismatch_envelope", {}) or {}).get("mismatch_envelopes", {}) or {}
    eps_sf = float(env.get("max_abs_mismatch", 0.0) or 0.0)

    above = bmap.get("above_target_budgets", []) or []
    below = bmap.get("below_target_budgets", []) or []

    if not (dm_star > 0 and left > 0 and right > 0 and eps_sf >= 0 and above and below):
        raise RuntimeError("Invalid upstream inputs for P2219 inflated budget intervals")

    inflated_above = []
    for r in above:
        sf = float(r["target_safety_factor"])
        d_lo = (sf - 1.0 - eps_sf) / right
        d_hi = (sf - 1.0 + eps_sf) / right
        inflated_above.append({
            "target_safety_factor": sf,
            "signed_delta_dm_interval": [d_lo, d_hi],
            "abs_dm_interval": [dm_star + d_lo, dm_star + d_hi],
        })

    inflated_below = []
    for r in below:
        sf = float(r["target_safety_factor"])
        mag_lo = (1.0 - sf - eps_sf) / left
        mag_hi = (1.0 - sf + eps_sf) / left
        d_lo = -mag_hi
        d_hi = -mag_lo
        inflated_below.append({
            "target_safety_factor": sf,
            "signed_delta_dm_interval": [d_lo, d_hi],
            "abs_dm_interval": [dm_star + d_lo, dm_star + d_hi],
        })

    payload = {
        "schema_version": "p2219_s1169_v1",
        "packet_id": "P2219",
        "stage_id": "S1169",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_TOLERANCE_INFLATED_SIGNED_BUDGET_INTERVALS",
        "strict_nu_branch_tolerance_inflated_signed_budget_intervals": {
            "interval_id": "STRICT_NU_BRANCH_TOLERANCE_INFLATED_SIGNED_BUDGET_INTERVALS_V1",
            "source_packets": [str(IN_2217.relative_to(ROOT)), str(IN_2218.relative_to(ROOT))],
            "reference_dm_star": dm_star,
            "inflation_abs_safety_factor_tolerance": eps_sf,
            "inflated_above_target_intervals": inflated_above,
            "inflated_below_target_intervals": inflated_below,
            "theorem_scope_limit": "local tolerance-inflated budget intervals only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2220_candidate",
            "goal": "validate interval coverage against direct neighborhood recompute table and export pass-rate certificate",
        },
        "gatekeeper_checks": {
            "inflated_interval_map_exported": True,
            "interval_endpoints_ordered": all(x["signed_delta_dm_interval"][0] <= x["signed_delta_dm_interval"][1] for x in inflated_above + inflated_below),
            "above_intervals_include_positive_side": all(x["signed_delta_dm_interval"][1] > 0 for x in inflated_above),
            "below_intervals_include_negative_side": all(x["signed_delta_dm_interval"][0] < 0 for x in inflated_below),
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
            "# P2219 S1169: strict nu-branch tolerance-inflated signed budget intervals",
            "",
            f"- reference dm*: `{dm_star:.12e}`",
            f"- inflation abs safety-factor tol: `{eps_sf:.12e}`",
            f"- interval endpoints ordered: `{all(x['signed_delta_dm_interval'][0] <= x['signed_delta_dm_interval'][1] for x in inflated_above + inflated_below)}`",
            f"- above intervals include positive side: `{all(x['signed_delta_dm_interval'][1] > 0 for x in inflated_above)}`",
            f"- below intervals include negative side: `{all(x['signed_delta_dm_interval'][0] < 0 for x in inflated_below)}`",
            "",
            "Local tolerance-inflated map only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
