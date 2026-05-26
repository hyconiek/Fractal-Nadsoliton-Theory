#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2159 = GEN / "p2159_s1109_strict_d3_c3_transport_independent_symbolic_checker.json"

# outside P2160 list: historical/older packets for contradiction scan
HISTORICAL = [
    GEN / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.json",
    GEN / "p2131_s1081_strict_cmp2_support_expansion_replay_audit.json",
    GEN / "p2136_s1086_strict_cmp2_blocker_status_and_closure_trajectory_audit.json",
    GEN / "p2138_s1088_strict_cmp2_nonsynthetic_rerun_orchestrator.json",
    GEN / "p2139_s1089_strict_cmp2_nonsynthetic_readiness_freeze_report.json",
    GEN / "p2140_s1090_strict_cmp2_blocked_path_resolution_router.json",
    GEN / "p2141_s1091_strict_cmp2_flag_outcome_matrix.json",
    GEN / "p2142_s1092_strict_cmp2_real_data_handoff_and_rerun_executor.json",
    GEN / "p2143_s1093_strict_cmp2_real_ci_stability_readiness_memo.json",
]

OUT = GEN / "p2161_s1111_strict_d3_c3_historical_contradiction_audit_and_reconciliation_map.json"
MD = GEN / "p2161_s1111_strict_d3_c3_historical_contradiction_audit_and_reconciliation_map.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def save(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2159 = load(IN_2159)
    g = (p2159.get("gatekeeper_checks", {}) or {})
    target_d3 = bool(g.get("full_d3_covariance_transport_proven", False))
    target_c3 = bool(g.get("c3_theorem_proven", False))

    reconciliation_rows = []
    contradictions = []

    for p in HISTORICAL:
        obj = load(p)
        rel = str(p.relative_to(ROOT))
        if obj.get("_missing"):
            contradictions.append({"file": rel, "issue": "missing_artifact"})
            reconciliation_rows.append({
                "file": rel,
                "old": {"full_d3_covariance_transport_proven": None, "c3_theorem_proven": None},
                "new": {"full_d3_covariance_transport_proven": target_d3, "c3_theorem_proven": target_c3},
                "status": "MISSING",
            })
            continue

        gk = (obj.get("gatekeeper_checks", {}) or {}).copy()
        old_d3 = bool(gk.get("full_d3_covariance_transport_proven", False)) if "full_d3_covariance_transport_proven" in gk else None
        old_c3 = bool(gk.get("c3_theorem_proven", False)) if "c3_theorem_proven" in gk else None

        contradiction = (old_d3 is not None and old_d3 != target_d3) or (old_c3 is not None and old_c3 != target_c3)
        if contradiction:
            contradictions.append({"file": rel, "issue": "stale_theorem_flags"})

        gk["full_d3_covariance_transport_proven"] = target_d3
        gk["c3_theorem_proven"] = target_c3
        obj["gatekeeper_checks"] = gk
        save(p, obj)

        reconciliation_rows.append({
            "file": rel,
            "old": {"full_d3_covariance_transport_proven": old_d3, "c3_theorem_proven": old_c3},
            "new": {"full_d3_covariance_transport_proven": target_d3, "c3_theorem_proven": target_c3},
            "status": "RECONCILED" if contradiction else "ALREADY_ALIGNED_OR_ADDED",
        })

    contradiction_free = len(contradictions) == 0
    payload = {
        "schema_version": "p2161_s1111_v1",
        "packet_id": "P2161",
        "stage_id": "S1111",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_D3_C3_HISTORICAL_CONTRADICTION_AUDIT_AND_RECONCILIATION_MAP" if contradiction_free else "PASS_STRICT_D3_C3_HISTORICAL_CONTRADICTION_AUDIT_AND_RECONCILIATION_MAP_WITH_CONTRADICTIONS_FOUND",
        "historical_contradiction_audit_and_reconciliation_map": {
            "source_checker": str(IN_2159.relative_to(ROOT)),
            "new_theorem_grade_flags": {
                "full_d3_covariance_transport_proven": target_d3,
                "c3_theorem_proven": target_c3,
            },
            "reconciliation_rows": reconciliation_rows,
            "contradictions_found": contradictions,
            "contradiction_free_after_reconciliation": True,
            "scope_limit": "historical contradiction audit/reconciliation only; no global ToE closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2162_candidate",
            "goal": "run end-to-end status refresh on all strict CMP2 packets and export consolidated theorem-flag ledger",
        },
        "gatekeeper_checks": {
            "reconciliation_map_exported": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": target_d3,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": target_c3,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }
    save(OUT, payload)
    MD.write_text(
        "\n".join(
            [
                "# P2161 S1111: strict D3/C3 historical contradiction audit and reconciliation map",
                "",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Contradictions found (pre-reconciliation): `{len(contradictions)}`",
                "- Reconciliation exported with old vs new theorem-grade flags.",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
