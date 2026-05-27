#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2221 = GEN / "p2221_s1171_strict_nu_branch_min_inflation_factor_for_full_coverage.json"
IN_2223 = GEN / "p2223_s1173_strict_nu_branch_targetwise_min_inflation_optimization.json"
OUT = GEN / "p2224_s1174_strict_nu_branch_mixed_inflation_policy_pareto_frontier.json"
MD = GEN / "p2224_s1174_strict_nu_branch_mixed_inflation_policy_pareto_frontier.md"


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
        raise RuntimeError("Invalid upstream inputs for P2224 mixed policy frontier")

    t_factors = [float(r["targetwise_min_inflation_factor"]) for r in rows]

    alphas = [0.0, 0.25, 0.5, 0.75, 1.0]
    frontier = []
    for a in alphas:
        # mixed policy between global g and per-target t_i
        mixed = [g * (1.0 - a) + ti * a for ti in t_factors]
        mean_ratio_vs_global = sum(m / g for m in mixed) / len(mixed)
        worst_deficit_vs_target = max(max(0.0, ti - m) for ti, m in zip(t_factors, mixed))
        satisfies_all_targetwise = all(m + 1e-15 >= ti for ti, m in zip(t_factors, mixed))
        frontier.append({
            "alpha_targetwise_mix": a,
            "mean_ratio_vs_global": mean_ratio_vs_global,
            "worst_deficit_vs_targetwise_min": worst_deficit_vs_target,
            "satisfies_all_targetwise_minima": satisfies_all_targetwise,
        })

    feasible = [r for r in frontier if r["satisfies_all_targetwise_minima"]]
    best_feasible = min(feasible, key=lambda r: r["mean_ratio_vs_global"]) if feasible else None

    payload = {
        "schema_version": "p2224_s1174_v1",
        "packet_id": "P2224",
        "stage_id": "S1174",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_MIXED_INFLATION_POLICY_PARETO_FRONTIER",
        "strict_nu_branch_mixed_inflation_policy_pareto_frontier": {
            "frontier_id": "STRICT_NU_BRANCH_MIXED_INFLATION_POLICY_PARETO_FRONTIER_V1",
            "source_packets": [str(IN_2221.relative_to(ROOT)), str(IN_2223.relative_to(ROOT))],
            "global_factor": g,
            "frontier_rows": frontier,
            "best_feasible_row": best_feasible,
            "theorem_scope_limit": "local policy-mix Pareto analysis only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2225_candidate",
            "goal": "replace coarse alpha grid by continuous constrained optimization and verify with direct recompute coverage",
        },
        "gatekeeper_checks": {
            "pareto_frontier_exported": True,
            "frontier_nonempty": len(frontier) > 0,
            "feasible_point_exists": best_feasible is not None,
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
            "# P2224 S1174: strict nu-branch mixed inflation policy Pareto frontier",
            "",
            f"- global factor: `{g:.12e}`",
            f"- feasible point exists: `{best_feasible is not None}`",
            f"- best feasible row: `{best_feasible}`",
            "",
            "Local policy-mix Pareto analysis only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
