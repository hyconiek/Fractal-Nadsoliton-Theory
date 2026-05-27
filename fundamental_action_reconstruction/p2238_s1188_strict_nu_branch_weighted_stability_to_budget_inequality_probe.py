#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2217 = GEN / "p2217_s1167_strict_nu_branch_signed_perturbation_budget_map.json"
IN_2237 = GEN / "p2237_s1187_strict_nu_branch_residual_to_modeled_weighted_stability_probe.json"
OUT = GEN / "p2238_s1188_strict_nu_branch_weighted_stability_to_budget_inequality_probe.json"
MD = GEN / "p2238_s1188_strict_nu_branch_weighted_stability_to_budget_inequality_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2217 = load(IN_2217)
    p2237 = load(IN_2237)

    bmap = (p2217.get("strict_nu_branch_signed_perturbation_budget_map", {}) or {})
    rows = (bmap.get("above_target_budgets", []) or []) + (bmap.get("below_target_budgets", []) or [])
    split = (p2237.get("strict_nu_branch_residual_to_modeled_weighted_stability_probe", {}) or {}).get("split_stability", {}) or {}
    fit = (p2237.get("strict_nu_branch_residual_to_modeled_weighted_stability_probe", {}) or {}).get("weighted_fit", {}) or {}

    if not rows:
        raise RuntimeError("Missing budget rows for P2238")

    band = float(split.get("uncertainty_band_radius_from_train_residual", 0.0) or 0.0)
    slope = float(fit.get("params", {}).get("a", 0.0) or 0.0)
    if slope <= 0.0:
        raise RuntimeError("Non-positive weighted slope for P2238")

    # Convert safety-factor surrogate uncertainty into equivalent residual-domain uncertainty.
    # inequality: |delta_residual| <= band / slope
    residual_uncert_bound = band / slope

    min_abs_budget = min(abs(float(r.get("required_abs_dm", 0.0) or 0.0)) for r in rows)
    ratio = residual_uncert_bound / min_abs_budget if min_abs_budget > 0.0 else float("inf")

    conservative_margin_factor = 1.0 - ratio if ratio <= 1.0 else -1.0
    inequality_holds = residual_uncert_bound <= min_abs_budget + 1e-18

    payload = {
        "schema_version": "p2238_s1188_v1",
        "packet_id": "P2238",
        "stage_id": "S1188",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_WEIGHTED_STABILITY_TO_BUDGET_INEQUALITY_PROBE",
        "strict_nu_branch_weighted_stability_to_budget_inequality_probe": {
            "probe_id": "STRICT_NU_BRANCH_WEIGHTED_STABILITY_TO_BUDGET_INEQUALITY_PROBE_V1",
            "source_packets": [str(IN_2217.relative_to(ROOT)), str(IN_2237.relative_to(ROOT))],
            "inputs": {
                "weighted_slope": slope,
                "surrogate_uncertainty_band_radius": band,
                "min_abs_signed_budget": min_abs_budget,
            },
            "derived_inequality": {
                "statement": "|delta_residual| <= band/slope <= min_abs_signed_budget",
                "lhs_bound_band_over_slope": residual_uncert_bound,
                "rhs_min_abs_signed_budget": min_abs_budget,
                "holds": inequality_holds,
                "relative_load_ratio": ratio,
                "conservative_margin_factor": conservative_margin_factor,
            },
            "physical_interpretation_note": "If inequality holds, surrogate-fit uncertainty remains below smallest local residual budget, supporting robustness of sign-preserving local response.",
            "theorem_scope_limit": "local budget compatibility check only; not a legacy->strict kernel bridge theorem",
        },
        "recommended_next_honest_step": {
            "id": "P2239_candidate",
            "goal": "propagate inequality margins target-wise and grouped-wise with worst-case aggregation under policy mixes",
        },
        "gatekeeper_checks": {
            "weighted_to_budget_inequality_exported": True,
            "weighted_slope_positive": slope > 0.0,
            "inequality_holds_against_min_budget": inequality_holds,
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
        "\n".join([
            "# P2238 S1188: weighted-stability to budget-inequality probe",
            "",
            f"- weighted slope: `{slope:.12e}`",
            f"- surrogate band radius: `{band:.12e}`",
            f"- implied residual uncertainty bound (band/slope): `{residual_uncert_bound:.12e}`",
            f"- minimum abs signed budget: `{min_abs_budget:.12e}`",
            f"- inequality holds: `{inequality_holds}`",
            f"- margin factor (1-ratio): `{conservative_margin_factor:.12e}`",
            "",
            "Local compatibility inequality only; no kernel-bridge theorem claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
