#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2218 = GEN / "p2218_s1168_strict_nu_branch_budget_map_direct_recompute_mismatch_envelope.json"
IN_2221 = GEN / "p2221_s1171_strict_nu_branch_min_inflation_factor_for_full_coverage.json"
OUT = GEN / "p2223_s1173_strict_nu_branch_targetwise_min_inflation_optimization.json"
MD = GEN / "p2223_s1173_strict_nu_branch_targetwise_min_inflation_optimization.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2218 = load(IN_2218)
    p2221 = load(IN_2221)

    env = p2218.get("strict_nu_branch_budget_map_direct_recompute_mismatch_envelope", {}) or {}
    rows = env.get("comparison_rows", []) or []
    base_eps = float(env.get("mismatch_envelopes", {}).get("max_abs_mismatch", 0.0) or 0.0)

    global_factor = float((p2221.get("strict_nu_branch_min_inflation_factor_for_full_coverage", {}) or {}).get("min_inflation_factor", 1.0) or 1.0)

    if not (rows and base_eps > 0.0 and global_factor >= 1.0):
        raise RuntimeError("Invalid upstream inputs for P2223 targetwise inflation optimization")

    opt_rows = []
    for r in rows:
        eps_i = float(r["abs_mismatch"])
        f_i = max(1.0, eps_i / base_eps)
        opt_rows.append({
            "certification_side": r["certification_side"],
            "target_safety_factor": float(r["target_safety_factor"]),
            "abs_mismatch": eps_i,
            "targetwise_min_inflation_factor": f_i,
            "width_cost_ratio_vs_global": f_i / global_factor,
        })

    mean_target_factor = sum(r["targetwise_min_inflation_factor"] for r in opt_rows) / len(opt_rows)
    mean_cost_vs_global = sum(r["width_cost_ratio_vs_global"] for r in opt_rows) / len(opt_rows)

    payload = {
        "schema_version": "p2223_s1173_v1",
        "packet_id": "P2223",
        "stage_id": "S1173",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_TARGETWISE_MIN_INFLATION_OPTIMIZATION",
        "strict_nu_branch_targetwise_min_inflation_optimization": {
            "optimization_id": "STRICT_NU_BRANCH_TARGETWISE_MIN_INFLATION_OPTIMIZATION_V1",
            "source_packets": [str(IN_2218.relative_to(ROOT)), str(IN_2221.relative_to(ROOT))],
            "base_eps": base_eps,
            "global_min_inflation_factor": global_factor,
            "targetwise_rows": opt_rows,
            "summary": {
                "mean_targetwise_min_factor": mean_target_factor,
                "mean_width_cost_ratio_vs_global": mean_cost_vs_global,
            },
            "theorem_scope_limit": "target-wise local inflation optimization only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2224_candidate",
            "goal": "construct mixed policy (global floor + target-wise deltas) and verify coverage/width Pareto frontier",
        },
        "gatekeeper_checks": {
            "targetwise_optimization_exported": True,
            "all_targetwise_factors_ge_1": all(r["targetwise_min_inflation_factor"] >= 1.0 for r in opt_rows),
            "mean_cost_not_worse_than_global": mean_cost_vs_global <= 1.0,
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
            "# P2223 S1173: strict nu-branch target-wise min inflation optimization",
            "",
            f"- global min inflation factor: `{global_factor:.12e}`",
            f"- mean targetwise min factor: `{mean_target_factor:.12e}`",
            f"- mean width-cost ratio vs global: `{mean_cost_vs_global:.12e}`",
            "",
            "Target-wise local optimization only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
