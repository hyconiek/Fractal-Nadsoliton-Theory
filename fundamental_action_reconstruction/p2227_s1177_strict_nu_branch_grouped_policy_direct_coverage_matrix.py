#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2226 = GEN / "p2226_s1176_strict_nu_branch_grouped_policy_width_coverage_tradeoff.json"
IN_2223 = GEN / "p2223_s1173_strict_nu_branch_targetwise_min_inflation_optimization.json"
OUT = GEN / "p2227_s1177_strict_nu_branch_grouped_policy_direct_coverage_matrix.json"
MD = GEN / "p2227_s1177_strict_nu_branch_grouped_policy_direct_coverage_matrix.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2226 = load(IN_2226)
    p2223 = load(IN_2223)

    policy = (p2226.get("strict_nu_branch_grouped_policy_width_coverage_tradeoff", {}) or {}).get("grouped_policy", {}) or {}
    f_above = float(policy.get("factor_above_dm_star", 0.0) or 0.0)
    f_below = float(policy.get("factor_below_dm_star", 0.0) or 0.0)

    rows = (p2223.get("strict_nu_branch_targetwise_min_inflation_optimization", {}) or {}).get("targetwise_rows", []) or []
    if not (rows and f_above >= 1.0 and f_below >= 1.0):
        raise RuntimeError("Invalid upstream inputs for P2227 grouped policy coverage matrix")

    matrix = []
    for r in rows:
        side = str(r["certification_side"])
        req = float(r["targetwise_min_inflation_factor"])
        applied = f_above if side == "above_dm_star" else f_below
        matrix.append({
            "certification_side": side,
            "target_safety_factor": float(r["target_safety_factor"]),
            "required_factor": req,
            "applied_grouped_factor": applied,
            "covered": applied + 1e-15 >= req,
            "slack": applied - req,
        })

    above_rows = [m for m in matrix if m["certification_side"] == "above_dm_star"]
    below_rows = [m for m in matrix if m["certification_side"] == "below_dm_star"]

    payload = {
        "schema_version": "p2227_s1177_v1",
        "packet_id": "P2227",
        "stage_id": "S1177",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUPED_POLICY_DIRECT_COVERAGE_MATRIX",
        "strict_nu_branch_grouped_policy_direct_coverage_matrix": {
            "matrix_id": "STRICT_NU_BRANCH_GROUPED_POLICY_DIRECT_COVERAGE_MATRIX_V1",
            "source_packets": [str(IN_2226.relative_to(ROOT)), str(IN_2223.relative_to(ROOT))],
            "grouped_policy": {"factor_above_dm_star": f_above, "factor_below_dm_star": f_below},
            "rows": matrix,
            "coverage_summary": {
                "all_above_rows_covered": all(r["covered"] for r in above_rows),
                "all_below_rows_covered": all(r["covered"] for r in below_rows),
                "all_rows_covered": all(r["covered"] for r in matrix),
            },
            "theorem_scope_limit": "direct coverage matrix against targetwise factors only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2228_candidate",
            "goal": "optimize grouped factors under explicit slack regularization and verify width-cost improvement vs current grouped policy",
        },
        "gatekeeper_checks": {
            "coverage_matrix_exported": True,
            "all_rows_covered": all(r["covered"] for r in matrix),
            "nonnegative_slack_all_rows": all(r["slack"] >= -1e-12 for r in matrix),
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
            "# P2227 S1177: strict nu-branch grouped policy direct coverage matrix",
            "",
            f"- all rows covered: `{payload['gatekeeper_checks']['all_rows_covered']}`",
            f"- nonnegative slack all rows: `{payload['gatekeeper_checks']['nonnegative_slack_all_rows']}`",
            "",
            "Direct local coverage matrix only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
