#!/usr/bin/env python3
"""P2054 S1004: adaptive budget allocator audit.

Compare fixed-vs-adaptive budget allocation using P2053 recommendation signals,
while keeping theorem-grade C3 claims explicitly open.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2054_s1004_strict_same_scheme_adaptive_budget_allocator_audit.json"
MD = GEN / "p2054_s1004_strict_same_scheme_adaptive_budget_allocator_audit.md"

SCHEMA_VERSION = "p2054_s1004_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
TOTAL_EXTRA_UNITS = 9000.0


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2053 = load("p2053_s1003_strict_same_scheme_censoring_aware_budget_profile_recommendation_audit.json")

    ready = p2053.get("result_kind") == "PASS_CENSORING_AWARE_BUDGET_PROFILE_RECOMMENDATION_AUDIT_WITH_TRACE__C3_STILL_OPEN"
    rows = p2053.get("ranked_profiles") or []

    n = len(rows)
    fixed_extra = TOTAL_EXTRA_UNITS / n if n > 0 else 0.0

    # Fixed allocation: equal extra compute to all profiles.
    fixed_rows = []
    for r in rows:
        fixed_rows.append(
            {
                "profile_id": r.get("profile_id"),
                "base_compute_units": float(r.get("compute_units_proxy", math.inf)),
                "extra_units_allocated": fixed_extra,
                "post_allocation_units": float(r.get("compute_units_proxy", 0.0)) + fixed_extra,
            }
        )

    # Adaptive allocation: allocate all extra units to the top ranked profile.
    adaptive_rows = []
    top_profile_id = rows[0].get("profile_id") if rows else None
    for r in rows:
        extra = TOTAL_EXTRA_UNITS if r.get("profile_id") == top_profile_id else 0.0
        adaptive_rows.append(
            {
                "profile_id": r.get("profile_id"),
                "base_compute_units": float(r.get("compute_units_proxy", math.inf)),
                "extra_units_allocated": extra,
                "post_allocation_units": float(r.get("compute_units_proxy", 0.0)) + extra,
                "priority_score": float(r.get("censoring_aware_priority_score", 0.0)),
            }
        )

    # "information/compute" proxy from P2053 score.
    fixed_info_proxy = sum(float(r.get("censoring_aware_priority_score", 0.0)) for r in rows)
    adaptive_info_proxy = fixed_info_proxy
    # Efficiency proxy: emphasize concentrated follow-up on highest-score profile.
    fixed_followup_signal = fixed_info_proxy / max(TOTAL_EXTRA_UNITS, 1.0)
    adaptive_followup_signal = (
        float(rows[0].get("censoring_aware_priority_score", 0.0)) / max(TOTAL_EXTRA_UNITS, 1.0)
        if rows
        else 0.0
    )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2054",
        "stage_id": "S1004",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_ADAPTIVE_BUDGET_ALLOCATOR_AUDIT_FIXED_VS_ADAPTIVE_WITH_TRACE__C3_STILL_OPEN"
            if ready and rows
            else "OPEN_ADAPTIVE_BUDGET_ALLOCATOR_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2053_present": p2053.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2053_json_sha256": file_sha256(GEN / "p2053_s1003_strict_same_scheme_censoring_aware_budget_profile_recommendation_audit.json"),
        },
        "allocation_parameters": {
            "total_extra_units": TOTAL_EXTRA_UNITS,
            "adaptive_rule": "allocate all extra units to highest censoring_aware_priority_score profile",
            "fixed_rule": "equal split across profiles",
            "information_compute_proxy": "censoring_aware_priority_score per added compute unit",
        },
        "fixed_allocation": fixed_rows,
        "adaptive_allocation": adaptive_rows,
        "comparison_summary": {
            "profiles_total": n,
            "top_profile_id": top_profile_id,
            "fixed_info_proxy": fixed_info_proxy,
            "adaptive_info_proxy": adaptive_info_proxy,
            "fixed_followup_signal": fixed_followup_signal,
            "adaptive_followup_signal": adaptive_followup_signal,
        },
        "c3_gate_update": {
            "C3_adaptive_budget_allocator_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "theorem-grade operator identity across full background family",
                "cross-background finite-part transport identity theorem",
                "global finite-part lock/cocycle theorem",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "rows_nonempty": len(rows) > 0,
            "budget_conserved_fixed": math.isclose(sum(r["extra_units_allocated"] for r in fixed_rows), TOTAL_EXTRA_UNITS, rel_tol=0, abs_tol=1e-9),
            "budget_conserved_adaptive": math.isclose(sum(r["extra_units_allocated"] for r in adaptive_rows), TOTAL_EXTRA_UNITS, rel_tol=0, abs_tol=1e-9),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2054 S1004: adaptive budget allocator audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Top profile: `{top_profile_id}`",
                "",
                "Compared fixed vs adaptive extra-budget allocation under censoring-aware ranking.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
