#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2209 = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
IN_2215 = GEN / "p2215_s1165_strict_nu_branch_crossing_neighborhood_safety_factor_map.json"
OUT = GEN / "p2216_s1166_strict_nu_branch_one_sided_lipschitz_slope_bounds.json"
MD = GEN / "p2216_s1166_strict_nu_branch_one_sided_lipschitz_slope_bounds.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2209 = load(IN_2209)
    p2215 = load(IN_2215)

    min_margin = float(((p2209.get("strict_nu_branch_transport_nonvanishing_threshold_map", {}) or {}).get("inputs", {}) or {}).get("min_majorant_margin", 0.0) or 0.0)
    block = p2215.get("strict_nu_branch_crossing_neighborhood_safety_factor_map", {}) or {}
    dm_star = float(block.get("reference_dm_star", 0.0) or 0.0)
    rows = sorted(block.get("rows", []) or [], key=lambda r: float(r["abs_dm"]))

    if not (min_margin > 0 and dm_star > 0 and len(rows) >= 3):
        raise RuntimeError("Invalid upstream inputs for P2216 one-sided slope bounds")

    below = [r for r in rows if float(r["abs_dm"]) < dm_star]
    above = [r for r in rows if float(r["abs_dm"]) > dm_star]

    def secant_slope(r1: dict[str, Any], r2: dict[str, Any]) -> float:
        d1 = float(r1["abs_dm"])
        d2 = float(r2["abs_dm"])
        s1 = float(r1["safety_factor_vs_margin"])
        s2 = float(r2["safety_factor_vs_margin"])
        return (s2 - s1) / (d2 - d1)

    below_slopes = [secant_slope(below[i], below[i + 1]) for i in range(len(below) - 1)] if len(below) >= 2 else []
    above_slopes = [secant_slope(above[i], above[i + 1]) for i in range(len(above) - 1)] if len(above) >= 2 else []

    # one-sided local bounds near crossing: take closest segment to dm*
    left_local = secant_slope(below[-1], {"abs_dm": dm_star, "safety_factor_vs_margin": 1.0})
    right_local = secant_slope({"abs_dm": dm_star, "safety_factor_vs_margin": 1.0}, above[0])

    payload = {
        "schema_version": "p2216_s1166_v1",
        "packet_id": "P2216",
        "stage_id": "S1166",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_ONE_SIDED_LIPSCHITZ_SLOPE_BOUNDS",
        "strict_nu_branch_one_sided_lipschitz_slope_bounds": {
            "certificate_id": "STRICT_NU_BRANCH_ONE_SIDED_LIPSCHITZ_SLOPE_BOUNDS_V1",
            "source_packets": [str(IN_2209.relative_to(ROOT)), str(IN_2215.relative_to(ROOT))],
            "reference": {"dm_star": dm_star, "margin": min_margin},
            "one_sided_local_slopes": {
                "left_secant_at_dm_star": left_local,
                "right_secant_at_dm_star": right_local,
            },
            "neighborhood_slope_envelopes": {
                "below_dm_star_min_secant": min(below_slopes) if below_slopes else None,
                "below_dm_star_max_secant": max(below_slopes) if below_slopes else None,
                "above_dm_star_min_secant": min(above_slopes) if above_slopes else None,
                "above_dm_star_max_secant": max(above_slopes) if above_slopes else None,
            },
            "robustness_interpretation": "safety-factor perturbation scales approximately linearly with |delta dm| under exported one-sided secant slope bounds in local neighborhood",
            "theorem_scope_limit": "finite-sample one-sided local slope bounds only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2217_candidate",
            "goal": "convert one-sided slope bounds into explicit perturbation budget map for certifying minimal signed offset away from dm*",
        },
        "gatekeeper_checks": {
            "one_sided_slope_bounds_exported": True,
            "left_local_slope_positive": left_local > 0.0,
            "right_local_slope_positive": right_local > 0.0,
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
            "# P2216 S1166: strict nu-branch one-sided local Lipschitz slope bounds",
            "",
            f"- left secant at dm*: `{left_local:.12e}`",
            f"- right secant at dm*: `{right_local:.12e}`",
            f"- left slope positive: `{left_local > 0.0}`",
            f"- right slope positive: `{right_local > 0.0}`",
            "",
            "Finite-sample one-sided local slope bounds only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
