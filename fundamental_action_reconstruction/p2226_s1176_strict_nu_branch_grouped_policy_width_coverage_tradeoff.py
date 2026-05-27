#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2221 = GEN / "p2221_s1171_strict_nu_branch_min_inflation_factor_for_full_coverage.json"
IN_2223 = GEN / "p2223_s1173_strict_nu_branch_targetwise_min_inflation_optimization.json"
OUT = GEN / "p2226_s1176_strict_nu_branch_grouped_policy_width_coverage_tradeoff.json"
MD = GEN / "p2226_s1176_strict_nu_branch_grouped_policy_width_coverage_tradeoff.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2221 = load(IN_2221)
    p2223 = load(IN_2223)

    g = float((p2221.get("strict_nu_branch_min_inflation_factor_for_full_coverage", {}) or {}).get("min_inflation_factor", 1.0) or 1.0)
    rows = (p2223.get("strict_nu_branch_targetwise_min_inflation_optimization", {}) or {}).get("targetwise_rows", []) or []
    if not (g >= 1.0 and rows):
        raise RuntimeError("Invalid upstream inputs for P2226 grouped policy tradeoff")

    # split targets by side and use per-side maxima as grouped factors
    above = [float(r["targetwise_min_inflation_factor"]) for r in rows if r["certification_side"] == "above_dm_star"]
    below = [float(r["targetwise_min_inflation_factor"]) for r in rows if r["certification_side"] == "below_dm_star"]
    if not (above and below):
        raise RuntimeError("Missing side rows for grouped policy")

    f_above = max(above)
    f_below = max(below)

    # metrics: normalized width cost vs global single-factor policy
    grouped_mean_ratio = (sum(f_above for _ in above) + sum(f_below for _ in below)) / (len(rows) * g)
    global_ratio = 1.0

    payload = {
        "schema_version": "p2226_s1176_v1",
        "packet_id": "P2226",
        "stage_id": "S1176",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUPED_POLICY_WIDTH_COVERAGE_TRADEOFF",
        "strict_nu_branch_grouped_policy_width_coverage_tradeoff": {
            "tradeoff_id": "STRICT_NU_BRANCH_GROUPED_POLICY_WIDTH_COVERAGE_TRADEOFF_V1",
            "source_packets": [str(IN_2221.relative_to(ROOT)), str(IN_2223.relative_to(ROOT))],
            "global_factor": g,
            "grouped_policy": {
                "factor_above_dm_star": f_above,
                "factor_below_dm_star": f_below,
            },
            "summary": {
                "global_mean_ratio_vs_global": global_ratio,
                "grouped_mean_ratio_vs_global": grouped_mean_ratio,
                "grouped_better_than_global": grouped_mean_ratio <= global_ratio,
            },
            "theorem_scope_limit": "grouped-policy local tradeoff only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2227_candidate",
            "goal": "validate grouped policy against direct recompute coverage table and export strict pass/fail matrix",
        },
        "gatekeeper_checks": {
            "grouped_tradeoff_exported": True,
            "grouped_factors_ge_1": f_above >= 1.0 and f_below >= 1.0,
            "grouped_not_worse_than_global": grouped_mean_ratio <= 1.0,
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
            "# P2226 S1176: strict nu-branch grouped-policy width/coverage tradeoff",
            "",
            f"- global factor: `{g:.12e}`",
            f"- grouped factor above: `{f_above:.12e}`",
            f"- grouped factor below: `{f_below:.12e}`",
            f"- grouped mean ratio vs global: `{grouped_mean_ratio:.12e}`",
            "",
            "Grouped local tradeoff only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
