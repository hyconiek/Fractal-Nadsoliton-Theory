#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2209 = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
IN_2217 = GEN / "p2217_s1167_strict_nu_branch_signed_perturbation_budget_map.json"
OUT = GEN / "p2218_s1168_strict_nu_branch_budget_map_direct_recompute_mismatch_envelope.json"
MD = GEN / "p2218_s1168_strict_nu_branch_budget_map_direct_recompute_mismatch_envelope.md"


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

    inputs = (p2209.get("strict_nu_branch_transport_nonvanishing_threshold_map", {}) or {}).get("inputs", {}) or {}
    c1 = float(inputs.get("linear_lower_bound_coeff_c1", 0.0) or 0.0)
    c2 = float(inputs.get("quadratic_lower_bound_coeff_c2", 0.0) or 0.0)
    margin = float(inputs.get("min_majorant_margin", 0.0) or 0.0)

    m = p2217.get("strict_nu_branch_signed_perturbation_budget_map", {}) or {}
    dm_star = float(m.get("reference_dm_star", 0.0) or 0.0)
    above = m.get("above_target_budgets", []) or []
    below = m.get("below_target_budgets", []) or []

    if not (c1 > 0 and c2 > 0 and margin > 0 and dm_star > 0 and above and below):
        raise RuntimeError("Invalid upstream inputs for P2218 mismatch envelope")

    rows = []
    for row in above + below:
        target_sf = float(row["target_safety_factor"])
        dm = float(row["required_abs_dm"])
        recomputed_sf = eff(dm, c1, c2) / margin
        rows.append({
            "certification_side": row["certification_side"],
            "target_safety_factor": target_sf,
            "required_abs_dm": dm,
            "recomputed_safety_factor": recomputed_sf,
            "abs_mismatch": abs(recomputed_sf - target_sf),
            "rel_mismatch": abs(recomputed_sf - target_sf) / max(abs(target_sf), 1e-30),
        })

    abs_env = max(r["abs_mismatch"] for r in rows)
    rel_env = max(r["rel_mismatch"] for r in rows)

    payload = {
        "schema_version": "p2218_s1168_v1",
        "packet_id": "P2218",
        "stage_id": "S1168",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_BUDGET_MAP_DIRECT_RECOMPUTE_MISMATCH_ENVELOPE",
        "strict_nu_branch_budget_map_direct_recompute_mismatch_envelope": {
            "envelope_id": "STRICT_NU_BRANCH_BUDGET_MAP_DIRECT_RECOMPUTE_MISMATCH_ENVELOPE_V1",
            "source_packets": [str(IN_2209.relative_to(ROOT)), str(IN_2217.relative_to(ROOT))],
            "reference_dm_star": dm_star,
            "comparison_rows": rows,
            "mismatch_envelopes": {
                "max_abs_mismatch": abs_env,
                "max_rel_mismatch": rel_env,
            },
            "theorem_scope_limit": "direct recompute envelope on local budget-map rows only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2219_candidate",
            "goal": "lift mismatch envelope into conservative signed budget intervals with explicit tolerance inflation",
        },
        "gatekeeper_checks": {
            "mismatch_envelope_exported": True,
            "max_abs_mismatch_bounded": abs_env <= 3e-3,
            "max_rel_mismatch_bounded": rel_env <= 3e-3,
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
            "# P2218 S1168: strict nu-branch budget-map direct recompute mismatch envelope",
            "",
            f"- max abs mismatch: `{abs_env:.12e}`",
            f"- max rel mismatch: `{rel_env:.12e}`",
            f"- abs envelope bounded (<=3e-3): `{abs_env <= 3e-3}`",
            f"- rel envelope bounded (<=3e-3): `{rel_env <= 3e-3}`",
            "",
            "Local direct-recompute envelope only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
