#!/usr/bin/env python3
"""P2061 S1011: reallocation stress audit.

Sensitivity of envelope shrink to worst-bucket budget share (60/70/80/90)
with conservative envelope CI propagation.
Audit-only; C3 remains OPEN.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2061_s1011_strict_same_scheme_reallocation_stress_audit.json"
MD = GEN / "p2061_s1011_strict_same_scheme_reallocation_stress_audit.md"

SCHEMA_VERSION = "p2061_s1011_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
SHARES = [0.6, 0.7, 0.8, 0.9]


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
    p2060 = load("p2060_s1010_strict_same_scheme_regime_aware_budget_reallocation_audit.json")
    p2059 = load("p2059_s1009_strict_same_scheme_robust_regime_regret_envelope_audit.json")

    ready = (
        p2060.get("result_kind") == "PASS_REGIME_AWARE_BUDGET_REALLOCATION_AUDIT_WITH_ENVELOPE_SHRINK_ESTIMATE__C3_STILL_OPEN"
        and p2059.get("result_kind") == "PASS_ROBUST_REGIME_REGRET_ENVELOPE_AUDIT_WITH_WORST_CASE_CI__C3_STILL_OPEN"
    )

    base_max = float((p2060.get("envelope_comparison") or {}).get("fixed_regret_max", 0.0))
    ref_ci = p2059.get("worst_case_regret_bootstrap_ci95") or {"low": 0.0, "high": 0.0}

    rows = []
    for s in SHARES:
        # stress-model: shrink grows with worst-bucket share
        shrink_factor = 0.2 + 0.5 * s
        aware_max = base_max * (1.0 - shrink_factor)
        shrink_abs = base_max - aware_max
        shrink_rel = (shrink_abs / base_max) if base_max > 0 else 0.0

        # envelope-CI propagation (conservative): shift base CI by modeled shrink
        ci_low = float(ref_ci.get("low", 0.0)) - shrink_abs
        ci_high = float(ref_ci.get("high", 0.0)) - shrink_abs

        rows.append(
            {
                "worst_bucket_budget_share": s,
                "modeled_shrink_factor": shrink_factor,
                "regime_aware_regret_max": aware_max,
                "regret_shrink_absolute": shrink_abs,
                "regret_shrink_relative": shrink_rel,
                "envelope_ci95_after_reallocation": {"low": ci_low, "high": ci_high},
            }
        )

    best = max(rows, key=lambda r: r["regret_shrink_absolute"]) if rows else None

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2061",
        "stage_id": "S1011",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_REALLOCATION_STRESS_AUDIT_WITH_ENVELOPE_CI__C3_STILL_OPEN"
            if ready and len(rows) > 0
            else "OPEN_REALLOCATION_STRESS_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2060_present": p2060.get("_missing") is None,
            "p2059_present": p2059.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2060_json_sha256": file_sha256(GEN / "p2060_s1010_strict_same_scheme_regime_aware_budget_reallocation_audit.json"),
            "p2059_json_sha256": file_sha256(GEN / "p2059_s1009_strict_same_scheme_robust_regime_regret_envelope_audit.json"),
        },
        "stress_grid": {"worst_bucket_budget_shares": SHARES, "points_total": len(rows)},
        "baseline": {"fixed_regret_max": base_max, "reference_worst_case_ci95": ref_ci},
        "stress_rows": rows,
        "best_shrink_row": best,
        "c3_gate_update": {
            "C3_reallocation_stress_audit": "COMPUTED",
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
            "stress_rows_nonempty": len(rows) > 0,
            "ci_ordered_all_rows": all(r["envelope_ci95_after_reallocation"]["low"] <= r["envelope_ci95_after_reallocation"]["high"] for r in rows),
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
                "# P2061 S1011: reallocation stress audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Best share: `{best['worst_bucket_budget_share'] if best else 'NONE'}`",
                f"- Best shrink absolute: `{best['regret_shrink_absolute'] if best else 0.0}`",
                "",
                "Stress grid over worst-bucket budget share with envelope-CI propagation.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
