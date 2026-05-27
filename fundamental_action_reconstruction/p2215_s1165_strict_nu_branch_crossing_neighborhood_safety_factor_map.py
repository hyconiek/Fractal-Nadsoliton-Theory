#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2209 = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
IN_2214 = GEN / "p2214_s1164_strict_nu_branch_closed_form_threshold_consistency_certificate.json"
OUT = GEN / "p2215_s1165_strict_nu_branch_crossing_neighborhood_safety_factor_map.json"
MD = GEN / "p2215_s1165_strict_nu_branch_crossing_neighborhood_safety_factor_map.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def eff(dm: float, c1: float, c2: float) -> float:
    return max(c1 * dm, c2 * dm * dm)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2209 = load(IN_2209)
    p2214 = load(IN_2214)

    inputs = (p2209.get("strict_nu_branch_transport_nonvanishing_threshold_map", {}) or {}).get("inputs", {}) or {}
    c1 = float(inputs.get("linear_lower_bound_coeff_c1", 0.0) or 0.0)
    c2 = float(inputs.get("quadratic_lower_bound_coeff_c2", 0.0) or 0.0)
    min_margin = float(inputs.get("min_majorant_margin", 0.0) or 0.0)

    cref = (p2214.get("strict_nu_branch_closed_form_threshold_consistency_certificate", {}) or {})
    dm_star = float((cref.get("bisection_reference", {}) or {}).get("dm_star_bisection", 0.0) or 0.0)

    if not (c1 > 0 and c2 > 0 and min_margin > 0 and dm_star > 0):
        raise RuntimeError("Invalid upstream inputs for P2215 safety-factor map")

    rel_offsets = [-0.20, -0.10, -0.05, -0.02, -0.01, 0.0, 0.01, 0.02, 0.05, 0.10, 0.20]
    rows = []
    for r in rel_offsets:
        dm = dm_star * (1.0 + r)
        lower = eff(dm, c1, c2)
        sf = lower / min_margin
        rows.append({
            "relative_offset_from_dm_star": r,
            "abs_dm": dm,
            "effective_lower_bound": lower,
            "safety_factor_vs_margin": sf,
            "margin_gap": lower - min_margin,
            "certifies_margin_exceedance": lower > min_margin,
        })

    below = [x for x in rows if x["relative_offset_from_dm_star"] < 0]
    above = [x for x in rows if x["relative_offset_from_dm_star"] > 0]

    payload = {
        "schema_version": "p2215_s1165_v1",
        "packet_id": "P2215",
        "stage_id": "S1165",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_CROSSING_NEIGHBORHOOD_SAFETY_FACTOR_MAP",
        "strict_nu_branch_crossing_neighborhood_safety_factor_map": {
            "map_id": "STRICT_NU_BRANCH_CROSSING_NEIGHBORHOOD_SAFETY_FACTOR_MAP_V1",
            "source_packets": [str(IN_2209.relative_to(ROOT)), str(IN_2214.relative_to(ROOT))],
            "reference_dm_star": dm_star,
            "rows": rows,
            "local_transition_summary": {
                "max_safety_factor_below_dm_star": max(x["safety_factor_vs_margin"] for x in below),
                "min_safety_factor_above_dm_star": min(x["safety_factor_vs_margin"] for x in above),
                "all_below_not_certified": all(not x["certifies_margin_exceedance"] for x in below),
                "all_above_certified": all(x["certifies_margin_exceedance"] for x in above),
            },
            "theorem_scope_limit": "local neighborhood map around modeled dm* only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2216_candidate",
            "goal": "export asymmetric one-sided local Lipschitz slope bounds around dm* for robustness-to-perturbation accounting",
        },
        "gatekeeper_checks": {
            "safety_factor_map_exported": True,
            "below_dm_star_not_certified": all(not x["certifies_margin_exceedance"] for x in below),
            "above_dm_star_certified": all(x["certifies_margin_exceedance"] for x in above),
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
            "# P2215 S1165: strict nu-branch crossing neighborhood safety-factor map",
            "",
            f"- reference dm*: `{dm_star:.12e}`",
            f"- all below dm* not certified: `{payload['gatekeeper_checks']['below_dm_star_not_certified']}`",
            f"- all above dm* certified: `{payload['gatekeeper_checks']['above_dm_star_certified']}`",
            "",
            "Local modeled-neighborhood map only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
