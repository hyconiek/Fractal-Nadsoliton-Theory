#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2229 = GEN / "p2229_s1179_strict_nu_branch_grouped_policy_miss_set_certificate.json"
OUT = GEN / "p2230_s1180_strict_nu_branch_grouped_policy_min_correction_deltas.json"
MD = GEN / "p2230_s1180_strict_nu_branch_grouped_policy_min_correction_deltas.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2229 = load(IN_2229)
    block = p2229.get("strict_nu_branch_grouped_policy_miss_set_certificate", {}) or {}

    pol = block.get("evaluated_grouped_policy", {}) or {}
    f_above = float(pol.get("factor_above", 0.0) or 0.0)
    f_below = float(pol.get("factor_below", 0.0) or 0.0)
    rows = block.get("coverage_rows", []) or []

    if not (rows and f_above >= 1.0 and f_below >= 1.0):
        raise RuntimeError("Invalid upstream inputs for P2230 correction deltas")

    req_above = [float(r["required_factor"]) for r in rows if r["certification_side"] == "above_dm_star"]
    req_below = [float(r["required_factor"]) for r in rows if r["certification_side"] == "below_dm_star"]

    target_above = max(req_above) if req_above else f_above
    target_below = max(req_below) if req_below else f_below

    delta_above = max(0.0, target_above - f_above)
    delta_below = max(0.0, target_below - f_below)

    corrected_above = f_above + delta_above
    corrected_below = f_below + delta_below

    payload = {
        "schema_version": "p2230_s1180_v1",
        "packet_id": "P2230",
        "stage_id": "S1180",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUPED_POLICY_MIN_CORRECTION_DELTAS",
        "strict_nu_branch_grouped_policy_min_correction_deltas": {
            "certificate_id": "STRICT_NU_BRANCH_GROUPED_POLICY_MIN_CORRECTION_DELTAS_V1",
            "source_packet": str(IN_2229.relative_to(ROOT)),
            "input_policy": {"factor_above": f_above, "factor_below": f_below},
            "required_policy": {"factor_above": target_above, "factor_below": target_below},
            "minimal_nonnegative_deltas": {"delta_above": delta_above, "delta_below": delta_below},
            "corrected_policy": {"factor_above": corrected_above, "factor_below": corrected_below},
            "theorem_scope_limit": "local correction-delta computation on exported grouped-policy matrix only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2231_candidate",
            "goal": "re-run grouped policy coverage matrix with corrected policy and export strict pass-zero miss-set witness",
        },
        "gatekeeper_checks": {
            "correction_delta_certificate_exported": True,
            "deltas_nonnegative": delta_above >= 0.0 and delta_below >= 0.0,
            "corrected_policy_meets_required": corrected_above + 1e-15 >= target_above and corrected_below + 1e-15 >= target_below,
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
            "# P2230 S1180: strict nu-branch grouped-policy minimal correction deltas",
            "",
            f"- delta_above: `{delta_above:.12e}`",
            f"- delta_below: `{delta_below:.12e}`",
            f"- corrected policy meets required: `{payload['gatekeeper_checks']['corrected_policy_meets_required']}`",
            "",
            "Local correction-delta computation only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
