#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2227 = GEN / "p2227_s1177_strict_nu_branch_grouped_policy_direct_coverage_matrix.json"
IN_2228 = GEN / "p2228_s1178_strict_nu_branch_grouped_factor_slack_regularized_optimization.json"
OUT = GEN / "p2229_s1179_strict_nu_branch_grouped_policy_miss_set_certificate.json"
MD = GEN / "p2229_s1179_strict_nu_branch_grouped_policy_miss_set_certificate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2227 = load(IN_2227)
    p2228 = load(IN_2228)

    matrix_rows = (p2227.get("strict_nu_branch_grouped_policy_direct_coverage_matrix", {}) or {}).get("rows", []) or []
    best = (p2228.get("strict_nu_branch_grouped_factor_slack_regularized_optimization", {}) or {}).get("best_fully_feasible_row", {}) or {}

    if not matrix_rows or not best:
        raise RuntimeError("Missing upstream grouped policy rows or best feasible row")

    f_above = float(best.get("factor_above", 0.0) or 0.0)
    f_below = float(best.get("factor_below", 0.0) or 0.0)
    if not (f_above >= 1.0 and f_below >= 1.0):
        raise RuntimeError("Invalid grouped factors from P2228 best feasible row")

    rows = []
    miss_set = []
    for r in matrix_rows:
        side = str(r["certification_side"])
        req = float(r["required_factor"])
        applied = f_above if side == "above_dm_star" else f_below
        covered = applied + 1e-15 >= req
        rec = {
            "certification_side": side,
            "target_safety_factor": float(r["target_safety_factor"]),
            "required_factor": req,
            "applied_factor": applied,
            "covered": covered,
            "deficit": max(0.0, req - applied),
        }
        rows.append(rec)
        if not covered:
            miss_set.append(rec)

    payload = {
        "schema_version": "p2229_s1179_v1",
        "packet_id": "P2229",
        "stage_id": "S1179",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUPED_POLICY_MISS_SET_CERTIFICATE",
        "strict_nu_branch_grouped_policy_miss_set_certificate": {
            "certificate_id": "STRICT_NU_BRANCH_GROUPED_POLICY_MISS_SET_CERTIFICATE_V1",
            "source_packets": [str(IN_2227.relative_to(ROOT)), str(IN_2228.relative_to(ROOT))],
            "evaluated_grouped_policy": {"factor_above": f_above, "factor_below": f_below},
            "coverage_rows": rows,
            "miss_set_rows": miss_set,
            "miss_set_size": len(miss_set),
            "theorem_scope_limit": "local grouped-policy coverage audit only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2230_candidate",
            "goal": "if miss_set nonempty, produce minimal correction deltas per side; if empty, freeze grouped policy as local certified control map",
        },
        "gatekeeper_checks": {
            "miss_set_certificate_exported": True,
            "all_rows_covered": len(miss_set) == 0,
            "miss_set_empty": len(miss_set) == 0,
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
            "# P2229 S1179: strict nu-branch grouped-policy miss-set certificate",
            "",
            f"- factor_above: `{f_above:.12e}`",
            f"- factor_below: `{f_below:.12e}`",
            f"- miss_set_size: `{len(miss_set)}`",
            f"- all rows covered: `{len(miss_set) == 0}`",
            "",
            "Local grouped-policy coverage audit only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
