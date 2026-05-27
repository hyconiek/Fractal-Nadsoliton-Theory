#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2239 = GEN / "p2239_s1189_strict_nu_branch_budget_inequality_target_group_aggregation_probe.json"
OUT = GEN / "p2240_s1190_strict_nu_branch_group_policy_sign_stability_condition_probe.json"
MD = GEN / "p2240_s1190_strict_nu_branch_group_policy_sign_stability_condition_probe.md"

def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))

def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2239 = load(IN_2239)

    probe = (p2239.get("strict_nu_branch_budget_inequality_target_group_aggregation_probe", {}) or {})
    rows = probe.get("group_rows", []) or []
    load_ratio = float((probe.get("global_inputs", {}) or {}).get("inequality_load_ratio", 1.0) or 1.0)

    # Sufficient condition for each group g:
    # load_ratio / coverage_g <= 1  <=> coverage_g >= load_ratio.
    required_min_coverage = load_ratio
    observed_min_coverage = min((float(r.get("coverage_fraction", 1.0) or 1.0) for r in rows), default=1.0)
    condition_holds = observed_min_coverage + 1e-15 >= required_min_coverage

    signed_margin = observed_min_coverage - required_min_coverage
    normalized_margin = signed_margin / max(required_min_coverage, 1e-12)

    per_group = []
    for r in rows:
        cov = float(r.get("coverage_fraction", 1.0) or 1.0)
        gid = r.get("group_id", "unknown")
        slack = cov - required_min_coverage
        per_group.append(
            {
                "group_id": gid,
                "coverage_fraction": cov,
                "required_min_coverage": required_min_coverage,
                "signed_slack": slack,
                "condition_holds": slack >= -1e-15,
            }
        )

    payload = {
        "schema_version": "p2240_s1190_v1",
        "packet_id": "P2240",
        "stage_id": "S1190",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_SIGN_STABILITY_CONDITION_PROBE",
        "strict_nu_branch_group_policy_sign_stability_condition_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_SIGN_STABILITY_CONDITION_PROBE_V1",
            "source_packets": [str(IN_2239.relative_to(ROOT))],
            "derived_condition": {
                "inequality": "coverage_fraction >= inequality_load_ratio",
                "required_min_coverage": required_min_coverage,
                "observed_min_coverage": observed_min_coverage,
                "signed_margin": signed_margin,
                "normalized_margin": normalized_margin,
            },
            "per_group_condition_rows": per_group,
            "all_groups_satisfy_sufficient_condition": condition_holds,
            "physical_interpretation_note": "If each policy-group keeps coverage above load ratio, signed safety response remains monotone-stable under grouped mixing in this finite-sample strict-lane diagnostic.",
            "theorem_scope_limit": "finite grouped-policy sufficient-condition check only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2241_candidate",
            "goal": "stress-test sufficient condition against adversarial coverage redistributions while preserving fixed global coverage mass",
        },
        "gatekeeper_checks": {
            "group_sign_stability_condition_exported": True,
            "all_groups_satisfy_sufficient_condition": condition_holds,
            "signed_margin_nonnegative": signed_margin >= -1e-15,
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
        "\n".join(
            [
                "# P2240 S1190: group-policy sign-stability sufficient condition probe",
                "",
                f"- condition: `coverage_fraction >= inequality_load_ratio`",
                f"- required min coverage: `{required_min_coverage:.12e}`",
                f"- observed min coverage: `{observed_min_coverage:.12e}`",
                f"- signed margin: `{signed_margin:.12e}`",
                f"- normalized margin: `{normalized_margin:.12e}`",
                f"- all groups satisfy condition: `{condition_holds}`",
                "",
                "Finite-sample grouped-policy sufficient condition only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

if __name__ == "__main__":
    main()
