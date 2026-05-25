#!/usr/bin/env python3
"""P2053 S1003: censoring-aware budget profile recommendation audit.

Consumes P2052 and exports a conservative recommendation map that separates
"no break detected" from "insufficient compute" regimes.
Audit-level only; C3 remains OPEN.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2053_s1003_strict_same_scheme_censoring_aware_budget_profile_recommendation_audit.json"
MD = GEN / "p2053_s1003_strict_same_scheme_censoring_aware_budget_profile_recommendation_audit.md"

SCHEMA_VERSION = "p2053_s1003_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


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
    p2052 = load("p2052_s1002_strict_same_scheme_budget_normalized_frontier_efficiency_audit.json")

    ready = p2052.get("result_kind") == "PASS_BUDGET_NORMALIZED_FRONTIER_EFFICIENCY_AUDIT_WITH_TRACE__C3_STILL_OPEN"
    rows = p2052.get("profile_efficiency_rows") or []

    costs = [float(r.get("compute_units_proxy", math.inf)) for r in rows]
    finite_costs = [c for c in costs if math.isfinite(c) and c > 0]
    cmin = min(finite_costs) if finite_costs else 1.0
    cmax = max(finite_costs) if finite_costs else 1.0

    ranked = []
    for r in rows:
        cost = float(r.get("compute_units_proxy", math.inf))
        detected = int(r.get("detected_frontier_flag", 0))
        # normalized cost in [0,1] for finite costs
        if math.isfinite(cost) and cmax > cmin:
            norm_cost = (cost - cmin) / (cmax - cmin)
        elif math.isfinite(cost):
            norm_cost = 0.0
        else:
            norm_cost = 1.0

        # conservative score: detection dominates; otherwise prefer lower cost.
        score = (1.0 if detected else 0.0) - 0.25 * norm_cost
        ranked.append(
            {
                "profile_id": r.get("profile_id"),
                "resolution_kind": r.get("resolution_kind"),
                "compute_units_proxy": cost,
                "detected_frontier_flag": detected,
                "normalized_compute_cost": norm_cost,
                "censoring_aware_priority_score": score,
            }
        )

    ranked.sort(key=lambda x: (x["detected_frontier_flag"], x["censoring_aware_priority_score"]), reverse=True)
    recommended = ranked[0] if ranked else None

    detection_count = sum(int(r.get("detected_frontier_flag", 0)) for r in rows)
    recommendation_mode = "frontier_detected" if detection_count > 0 else "cost_minimization_under_censoring"

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2053",
        "stage_id": "S1003",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_CENSORING_AWARE_BUDGET_PROFILE_RECOMMENDATION_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready and ranked
            else "OPEN_CENSORING_AWARE_BUDGET_PROFILE_RECOMMENDATION_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2052_present": p2052.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2052_json_sha256": file_sha256(GEN / "p2052_s1002_strict_same_scheme_budget_normalized_frontier_efficiency_audit.json"),
        },
        "recommendation_policy": {
            "mode": recommendation_mode,
            "score_formula": "detected_frontier_flag - 0.25 * normalized_compute_cost",
            "interpretation": "if no frontier is detected, minimize compute under censoring while preserving explicit OPEN status",
        },
        "ranked_profiles": ranked,
        "recommended_profile": recommended,
        "c3_gate_update": {
            "C3_censoring_aware_budget_recommendation": "COMPUTED",
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
            "rows_nonempty": len(ranked) > 0,
            "has_recommendation": recommended is not None,
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
                "# P2053 S1003: censoring-aware budget profile recommendation audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Recommendation mode: `{recommendation_mode}`",
                "",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
