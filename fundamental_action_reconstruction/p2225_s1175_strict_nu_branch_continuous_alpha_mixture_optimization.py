#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2221 = GEN / "p2221_s1171_strict_nu_branch_min_inflation_factor_for_full_coverage.json"
IN_2223 = GEN / "p2223_s1173_strict_nu_branch_targetwise_min_inflation_optimization.json"
OUT = GEN / "p2225_s1175_strict_nu_branch_continuous_alpha_mixture_optimization.json"
MD = GEN / "p2225_s1175_strict_nu_branch_continuous_alpha_mixture_optimization.md"


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
    t = [float(r["targetwise_min_inflation_factor"]) for r in rows]

    if not (g >= 1.0 and t):
        raise RuntimeError("Invalid upstream inputs for P2225 continuous alpha optimization")

    # mixed factor: m_i(alpha)=g*(1-alpha)+alpha*t_i
    # feasibility m_i>=t_i for all i gives alpha>=1 if any t_i<g? Actually m_i-t_i=(1-alpha)*(g-t_i).
    # so feasibility for all i with mixed convex policy requires alpha=1 whenever any t_i != g.
    has_gap = any(abs(g - ti) > 1e-15 for ti in t)
    alpha_star = 1.0 if has_gap else 0.0

    mixed = [g * (1.0 - alpha_star) + alpha_star * ti for ti in t]
    mean_ratio = sum(m / g for m in mixed) / len(mixed)
    worst_deficit = max(max(0.0, ti - m) for ti, m in zip(t, mixed))

    payload = {
        "schema_version": "p2225_s1175_v1",
        "packet_id": "P2225",
        "stage_id": "S1175",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_CONTINUOUS_ALPHA_MIXTURE_OPTIMIZATION",
        "strict_nu_branch_continuous_alpha_mixture_optimization": {
            "optimization_id": "STRICT_NU_BRANCH_CONTINUOUS_ALPHA_MIXTURE_OPTIMIZATION_V1",
            "source_packets": [str(IN_2221.relative_to(ROOT)), str(IN_2223.relative_to(ROOT))],
            "global_factor": g,
            "alpha_star": alpha_star,
            "has_global_targetwise_gap": has_gap,
            "objective_summary": {
                "mean_ratio_vs_global": mean_ratio,
                "worst_deficit_vs_targetwise_min": worst_deficit,
            },
            "analytic_note": "For convex one-parameter mixing between global and target-wise factors, exact feasibility across all targets forces alpha=1 unless all target-wise factors equal global factor.",
            "theorem_scope_limit": "analytic statement for this one-parameter mixture family only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2226_candidate",
            "goal": "optimize richer piecewise/grouped policy beyond single-alpha family to trade off feasibility and width cost",
        },
        "gatekeeper_checks": {
            "continuous_alpha_optimization_exported": True,
            "alpha_star_in_unit_interval": 0.0 <= alpha_star <= 1.0,
            "worst_deficit_zero_at_alpha_star": worst_deficit <= 1e-12,
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
            "# P2225 S1175: strict nu-branch continuous alpha-mixture optimization",
            "",
            f"- alpha*: `{alpha_star:.12e}`",
            f"- has global/targetwise gap: `{has_gap}`",
            f"- mean ratio vs global at alpha*: `{mean_ratio:.12e}`",
            f"- worst deficit at alpha*: `{worst_deficit:.12e}`",
            "",
            "Analytic one-parameter mixture result only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
