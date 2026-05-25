#!/usr/bin/env python3
"""P2052 S1002: budget-normalized frontier efficiency audit.

Consumes P2051 budget profiles and exports cost/effect metrics per profile,
including frontier_detection_per_compute_unit.

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
OUT = GEN / "p2052_s1002_strict_same_scheme_budget_normalized_frontier_efficiency_audit.json"
MD = GEN / "p2052_s1002_strict_same_scheme_budget_normalized_frontier_efficiency_audit.md"

SCHEMA_VERSION = "p2052_s1002_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
BASIS = ["R2_bar", "Ric2_bar", "Riem2_bar"]


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_bool(v: Any) -> bool:
    return bool(v is True)


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2051 = load("p2051_s1001_strict_same_scheme_budget_scaling_audit.json")
    checks_2051 = p2051.get("gatekeeper_checks") or {}

    ready = (
        p2051.get("result_kind")
        == "PASS_BUDGET_SCALING_AUDIT_WITH_SENSITIVITY_MAP__C3_STILL_OPEN"
        and as_bool(checks_2051.get("profiles_nonempty"))
        and as_bool(checks_2051.get("counts_consistent"))
    )

    profiles = p2051.get("budget_profiles") or []

    rows = []
    best_eff = -math.inf
    best_profile = None

    for row in profiles:
        scan_steps = int(row.get("scan_steps_used", 0))
        sign_budget = int(row.get("final_sign_pattern_budget", 0))
        vectors_per_bucket = int(row.get("vectors_per_bucket", 0))

        # proxy compute units: scan_steps * sign_budget * vectors_per_bucket
        compute_units = max(scan_steps * sign_budget * max(vectors_per_bucket, 1), 1)
        detected = 1 if row.get("resolution_kind") == "detected_frontier" else 0
        eff = detected / compute_units

        out = {
            "profile_id": row.get("profile_id"),
            "resolution_kind": row.get("resolution_kind"),
            "scan_steps_used": scan_steps,
            "final_sign_pattern_budget": sign_budget,
            "vectors_per_bucket": vectors_per_bucket,
            "compute_units_proxy": compute_units,
            "detected_frontier_flag": detected,
            "frontier_detection_per_compute_unit": eff,
            "horizon_reached": bool(row.get("horizon_reached", False)),
            "pattern_budget_reached": bool(row.get("pattern_budget_reached", False)),
        }
        rows.append(out)

        if eff > best_eff:
            best_eff = eff
            best_profile = out

    detected_count = sum(r["detected_frontier_flag"] for r in rows)
    mean_eff = sum(r["frontier_detection_per_compute_unit"] for r in rows) / len(rows) if rows else math.inf

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2052",
        "stage_id": "S1002",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_BUDGET_NORMALIZED_FRONTIER_EFFICIENCY_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready and rows
            else "OPEN_BUDGET_NORMALIZED_FRONTIER_EFFICIENCY_AUDIT_BLOCKED"
        ),
        "route": "strict_only_budget_normalized_efficiency_audit",
        "depends_on": {
            "p2051_present": p2051.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2051_json_sha256": file_sha256(GEN / "p2051_s1001_strict_same_scheme_budget_scaling_audit.json"),
        },
        "audit_scope": {
            "controlled_background_pair": (p2051.get("audit_scope") or {}).get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "profiles_compared": len(rows),
            "compute_unit_formula": "scan_steps_used * final_sign_pattern_budget * vectors_per_bucket",
        },
        "profile_efficiency_rows": rows,
        "efficiency_summary": {
            "detected_frontier_count": detected_count,
            "profiles_total": len(rows),
            "mean_frontier_detection_per_compute_unit": mean_eff,
            "best_profile": best_profile,
        },
        "c3_gate_update": {
            "C3_budget_normalized_efficiency_audit": "COMPUTED",
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
            "detected_count_consistent": detected_count <= len(rows),
            "mean_eff_finite": math.isfinite(mean_eff),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2052 S1002: budget-normalized frontier efficiency audit",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Exported cost/effect rows per budget profile and metric",
        "`frontier_detection_per_compute_unit`.",
        "",
        "## Gate update",
        "",
        "- `C3`: efficiency audit computed.",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
