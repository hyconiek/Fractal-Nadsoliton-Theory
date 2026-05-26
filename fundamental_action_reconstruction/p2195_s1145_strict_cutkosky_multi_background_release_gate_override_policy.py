#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2194 = GEN / "p2194_s1144_strict_cutkosky_multi_background_tolerance_stability_audit.json"
OUT = GEN / "p2195_s1145_strict_cutkosky_multi_background_release_gate_override_policy.json"
MD = GEN / "p2195_s1145_strict_cutkosky_multi_background_release_gate_override_policy.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2194 = load(IN_2194)
    audit = p2194.get("strict_cutkosky_multi_background_tolerance_stability_audit", {}) or {}

    class_stability = audit.get("class_stability", {}) or {}
    stable = bool(class_stability.get("stable", False))
    class_set = class_stability.get("class_set", []) or []
    drift = bool(audit.get("background_drift_detected", True))

    if drift:
        override_decision = "FORCE_WARN_REVIEW"
        release_gate_override = "OPEN_WITH_DRIFT_NOTICE"
    elif stable and class_set == ["TOL_A_STRICT"]:
        override_decision = "NO_OVERRIDE"
        release_gate_override = "OPEN"
    elif stable and class_set == ["TOL_B_MONITORED"]:
        override_decision = "PIN_WARN"
        release_gate_override = "OPEN_WITH_NOTICE"
    else:
        override_decision = "FORCE_FAIL"
        release_gate_override = "BLOCK"

    policy = {
        "policy_id": "STRICT_CUTKOSKY_MULTI_BACKGROUND_RELEASE_GATE_OVERRIDE_POLICY_V1",
        "source_packet": str(IN_2194.relative_to(ROOT)),
        "drift_clause_table": {
            "drift_detected_true": {
                "override_decision": "FORCE_WARN_REVIEW",
                "release_gate_override": "OPEN_WITH_DRIFT_NOTICE",
                "rule": "require mandatory human review and freeze tolerance-class promotion",
            },
            "stable_tol_a_only": {
                "override_decision": "NO_OVERRIDE",
                "release_gate_override": "OPEN",
                "rule": "keep gate open under strict stable tolerance",
            },
            "stable_tol_b_only": {
                "override_decision": "PIN_WARN",
                "release_gate_override": "OPEN_WITH_NOTICE",
                "rule": "keep warning and require monitored rerun schedule",
            },
            "mixed_or_blocking_classes": {
                "override_decision": "FORCE_FAIL",
                "release_gate_override": "BLOCK",
                "rule": "block release until tolerance class inconsistency is resolved",
            },
        },
        "inputs": {
            "stable": stable,
            "class_set": class_set,
            "background_drift_detected": drift,
        },
        "output": {
            "override_decision": override_decision,
            "release_gate_override": release_gate_override,
        },
    }

    result_kind = "PASS_STRICT_CUTKOSKY_MULTI_BACKGROUND_RELEASE_GATE_OVERRIDE_POLICY"

    payload = {
        "schema_version": "p2195_s1145_v1",
        "packet_id": "P2195",
        "stage_id": "S1145",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_cutkosky_multi_background_release_gate_override_policy": policy,
        "recommended_next_honest_step": {
            "id": "P2196_candidate",
            "goal": "replay override policy on historical multi-background runs and export drift-frequency trend",
        },
        "gatekeeper_checks": {
            "override_policy_exported": True,
            "drift_clause_table_explicit": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2195 S1145: strict Cutkosky multi-background release-gate override policy",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Override decision: `{override_decision}`",
                f"- Release gate override: `{release_gate_override}`",
                f"- Drift detected: `{drift}`",
                "",
                "This is a drift-clause override policy only, not full Cutkosky closure.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
