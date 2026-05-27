#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2209 = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
IN_2217 = GEN / "p2217_s1167_strict_nu_branch_signed_perturbation_budget_map.json"
IN_2218 = GEN / "p2218_s1168_strict_nu_branch_budget_map_direct_recompute_mismatch_envelope.json"
OUT = GEN / "p2221_s1171_strict_nu_branch_min_inflation_factor_for_full_coverage.json"
MD = GEN / "p2221_s1171_strict_nu_branch_min_inflation_factor_for_full_coverage.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def eff(dm: float, c1: float, c2: float) -> float:
    return max(c1 * dm, c2 * dm * dm)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2209 = load(IN_2209)
    p2217 = load(IN_2217)
    p2218 = load(IN_2218)

    inputs = (p2209.get("strict_nu_branch_transport_nonvanishing_threshold_map", {}) or {}).get("inputs", {}) or {}
    c1 = float(inputs.get("linear_lower_bound_coeff_c1", 0.0) or 0.0)
    c2 = float(inputs.get("quadratic_lower_bound_coeff_c2", 0.0) or 0.0)
    margin = float(inputs.get("min_majorant_margin", 0.0) or 0.0)

    bmap = p2217.get("strict_nu_branch_signed_perturbation_budget_map", {}) or {}
    dm_star = float(bmap.get("reference_dm_star", 0.0) or 0.0)
    slopes = bmap.get("one_sided_local_slopes", {}) or {}
    left = float(slopes.get("left_secant_at_dm_star", 0.0) or 0.0)
    right = float(slopes.get("right_secant_at_dm_star", 0.0) or 0.0)
    above = bmap.get("above_target_budgets", []) or []
    below = bmap.get("below_target_budgets", []) or []

    env = (p2218.get("strict_nu_branch_budget_map_direct_recompute_mismatch_envelope", {}) or {}).get("mismatch_envelopes", {}) or {}
    eps_base = float(env.get("max_abs_mismatch", 0.0) or 0.0)

    if not (c1 > 0 and c2 > 0 and margin > 0 and dm_star > 0 and left > 0 and right > 0 and above and below and eps_base >= 0):
        raise RuntimeError("Invalid upstream inputs for P2221 min-inflation factor search")

    comp_rows = (p2218.get("strict_nu_branch_budget_map_direct_recompute_mismatch_envelope", {}) or {}).get("comparison_rows", []) or []
    required_eps = [float(r["abs_mismatch"]) for r in comp_rows]
    eps_required = max(required_eps) if required_eps else 0.0

    def coverage_ok(factor: float) -> tuple[bool, bool]:
        eps_use = eps_base * factor

        def covered(sf_target: float, slope: float, side: str) -> bool:
            if side == "above":
                d_lo = (sf_target - 1.0 - eps_use) / slope
                d_hi = (sf_target - 1.0 + eps_use) / slope
            else:
                mag_lo = (1.0 - sf_target - eps_use) / slope
                mag_hi = (1.0 - sf_target + eps_use) / slope
                d_lo, d_hi = -mag_hi, -mag_lo
            a, b = dm_star + d_lo, dm_star + d_hi
            lo, hi = (a, b) if a <= b else (b, a)
            sf_lo = eff(lo, c1, c2) / margin
            sf_hi = eff(hi, c1, c2) / margin
            return min(sf_lo, sf_hi) <= sf_target <= max(sf_lo, sf_hi)

        above_ok = all(covered(float(r["target_safety_factor"]), right, "above") for r in above)
        below_ok = all(covered(float(r["target_safety_factor"]), left, "below") for r in below)
        return above_ok, below_ok

    lo_f, hi_f = 1.0, 1.0
    aok, bok = coverage_ok(hi_f)
    while not (aok and bok) and hi_f < 1e6:
        hi_f *= 2.0
        aok, bok = coverage_ok(hi_f)

    if aok and bok:
        for _ in range(80):
            mid = 0.5 * (lo_f + hi_f)
            ma, mb = coverage_ok(mid)
            if ma and mb:
                hi_f = mid
            else:
                lo_f = mid
        min_inflation_factor = hi_f
        above_ok, below_ok = coverage_ok(min_inflation_factor)
    else:
        min_inflation_factor = hi_f
        above_ok, below_ok = aok, bok

    payload = {
        "schema_version": "p2221_s1171_v1",
        "packet_id": "P2221",
        "stage_id": "S1171",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_MIN_INFLATION_FACTOR_FOR_FULL_COVERAGE",
        "strict_nu_branch_min_inflation_factor_for_full_coverage": {
            "certificate_id": "STRICT_NU_BRANCH_MIN_INFLATION_FACTOR_FOR_FULL_COVERAGE_V1",
            "source_packets": [str(IN_2209.relative_to(ROOT)), str(IN_2217.relative_to(ROOT)), str(IN_2218.relative_to(ROOT))],
            "base_abs_mismatch_envelope": eps_base,
            "required_abs_mismatch_envelope": eps_required,
            "min_inflation_factor": min_inflation_factor,
            "verification": {
                "all_above_targets_covered": above_ok,
                "all_below_targets_covered": below_ok,
                "all_targets_covered": above_ok and below_ok,
            },
            "theorem_scope_limit": "local endpoint-band inflation calibration only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2222_candidate",
            "goal": "export calibrated inflated intervals with min_inflation_factor and compare width-cost vs base inflation",
        },
        "gatekeeper_checks": {
            "min_inflation_factor_exported": True,
            "inflation_factor_ge_1": min_inflation_factor >= 1.0,
            "all_targets_covered_with_min_factor": above_ok and below_ok,
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
            "# P2221 S1171: strict nu-branch min inflation factor for full coverage",
            "",
            f"- base eps: `{eps_base:.12e}`",
            f"- required eps: `{eps_required:.12e}`",
            f"- min inflation factor: `{min_inflation_factor:.12e}`",
            f"- all targets covered with min factor: `{above_ok and below_ok}`",
            "",
            "Local inflation calibration only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
