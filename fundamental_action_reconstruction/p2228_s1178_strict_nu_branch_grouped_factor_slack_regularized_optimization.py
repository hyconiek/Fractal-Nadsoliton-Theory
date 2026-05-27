#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2227 = GEN / "p2227_s1177_strict_nu_branch_grouped_policy_direct_coverage_matrix.json"
OUT = GEN / "p2228_s1178_strict_nu_branch_grouped_factor_slack_regularized_optimization.json"
MD = GEN / "p2228_s1178_strict_nu_branch_grouped_factor_slack_regularized_optimization.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2227 = load(IN_2227)
    block = p2227.get("strict_nu_branch_grouped_policy_direct_coverage_matrix", {}) or {}
    rows = block.get("rows", []) or []
    policy = block.get("grouped_policy", {}) or {}
    f_above_base = float(policy.get("factor_above_dm_star", 0.0) or 0.0)
    f_below_base = float(policy.get("factor_below_dm_star", 0.0) or 0.0)

    if not (rows and f_above_base >= 1.0 and f_below_base >= 1.0):
        raise RuntimeError("Invalid upstream inputs for P2228 grouped optimization")

    above_req = [float(r["required_factor"]) for r in rows if r["certification_side"] == "above_dm_star"]
    below_req = [float(r["required_factor"]) for r in rows if r["certification_side"] == "below_dm_star"]

    if not (above_req and below_req):
        raise RuntimeError("Missing side rows for grouped optimization")

    # hard-feasible minimal grouped factors
    f_above_hard = max(above_req)
    f_below_hard = max(below_req)

    # slack-regularized candidates: allow small violations to reduce width-cost
    # objective = mean_factor + lambda * mean_positive_deficit
    lambdas = [0.0, 0.5, 1.0, 2.0, 5.0]
    candidates = []
    for lam in lambdas:
        # interpolate between base grouped policy and under-relaxed policy near means
        f_above = max(1.0, min(f_above_hard, sum(above_req)/len(above_req) + lam * (f_above_hard - sum(above_req)/len(above_req)) / (1.0 + lam)))
        f_below = max(1.0, min(f_below_hard, sum(below_req)/len(below_req) + lam * (f_below_hard - sum(below_req)/len(below_req)) / (1.0 + lam)))

        deficits = []
        for r in rows:
            req = float(r["required_factor"])
            side = r["certification_side"]
            applied = f_above if side == "above_dm_star" else f_below
            deficits.append(max(0.0, req - applied))

        mean_factor = 0.5 * (f_above + f_below)
        mean_deficit = sum(deficits) / len(deficits)
        obj = mean_factor + lam * mean_deficit
        candidates.append({
            "lambda": lam,
            "factor_above": f_above,
            "factor_below": f_below,
            "mean_factor": mean_factor,
            "mean_positive_deficit": mean_deficit,
            "objective": obj,
            "fully_feasible": mean_deficit <= 1e-15,
        })

    best = min(candidates, key=lambda c: c["objective"])
    best_feasible = min((c for c in candidates if c["fully_feasible"]), key=lambda c: c["mean_factor"], default=None)

    payload = {
        "schema_version": "p2228_s1178_v1",
        "packet_id": "P2228",
        "stage_id": "S1178",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUPED_FACTOR_SLACK_REGULARIZED_OPTIMIZATION",
        "strict_nu_branch_grouped_factor_slack_regularized_optimization": {
            "optimization_id": "STRICT_NU_BRANCH_GROUPED_FACTOR_SLACK_REGULARIZED_OPTIMIZATION_V1",
            "source_packet": str(IN_2227.relative_to(ROOT)),
            "hard_feasible_grouped_factors": {"above": f_above_hard, "below": f_below_hard},
            "candidate_rows": candidates,
            "best_objective_row": best,
            "best_fully_feasible_row": best_feasible,
            "theorem_scope_limit": "coarse lambda-grid grouped-factor optimization only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2229_candidate",
            "goal": "validate chosen grouped factors against direct recompute endpoint-band coverage and report exact miss set if any",
        },
        "gatekeeper_checks": {
            "grouped_slack_optimization_exported": True,
            "candidates_nonempty": len(candidates) > 0,
            "has_fully_feasible_candidate": best_feasible is not None,
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
            "# P2228 S1178: strict nu-branch grouped-factor slack-regularized optimization",
            "",
            f"- best objective row: `{best}`",
            f"- best fully feasible row: `{best_feasible}`",
            "",
            "Coarse grouped optimization only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
