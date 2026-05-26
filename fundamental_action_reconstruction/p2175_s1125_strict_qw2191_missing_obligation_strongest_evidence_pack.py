#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2174 = GEN / "p2174_s1124_strict_qw2191_evidentiary_trace_normalizer_and_missing_obligation_witness_pack.json"
IN_2173 = GEN / "p2173_s1123_strict_qw2191_obligation_independent_crosschecker.json"
OUT = GEN / "p2175_s1125_strict_qw2191_missing_obligation_strongest_evidence_pack.json"
MD = GEN / "p2175_s1125_strict_qw2191_missing_obligation_strongest_evidence_pack.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def build_strongest_evidence_artifact(obligation: dict[str, Any]) -> dict[str, Any]:
    oid = obligation.get("id", "unknown")
    return {
        "obligation_id": oid,
        "artifact_id": f"EVID_{oid}_STRONGEST_V1",
        "artifact_type": "structured_evidence_packet",
        "evidence_claim": "obligation pass is backed by explicit witness+validator chain",
        "evidence_basis": obligation.get("evidence_items", []),
        "premise_basis": obligation.get("premise_items", []),
        "evidence_strength": obligation.get("evidence_strength", "missing_or_partial"),
        "admission_rule": "strongest only if status is PASS_WITH_TRACE and cross-check remains consistent",
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2174 = load(IN_2174)
    p2173 = load(IN_2173)

    pack = p2174.get("strict_qw2191_evidentiary_trace_normalizer_and_missing_obligation_witness_pack", {}) or {}
    normalized = pack.get("normalized_obligation_traces", []) or []
    missing = pack.get("missing_or_partial_obligations", []) or []

    strongest_for_passed = [build_strongest_evidence_artifact(o) for o in normalized if o.get("evidence_strength") == "strong_evidence"]
    remediation_for_missing = [
        {
            "obligation_id": o.get("id", "unknown"),
            "required_next_artifact": f"EVID_{o.get('id', 'unknown')}_STRONGEST_V1",
            "reason": "still classified as missing_or_partial in P2174",
        }
        for o in missing
    ]

    all_consistent = bool((p2173.get("gatekeeper_checks", {}) or {}).get("all_consistent", False))
    no_missing = len(missing) == 0

    result_kind = (
        "PASS_STRICT_QW2191_STRONGEST_EVIDENCE_PACK_COMPLETE_WITH_TRACE"
        if no_missing and all_consistent
        else "OPEN_STRICT_QW2191_STRONGEST_EVIDENCE_PACK_PARTIAL_WITH_TRACE"
    )

    payload = {
        "schema_version": "p2175_s1125_v1",
        "packet_id": "P2175",
        "stage_id": "S1125",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_missing_obligation_strongest_evidence_pack": {
            "source_trace_normalizer": str(IN_2174.relative_to(ROOT)),
            "source_crosschecker": str(IN_2173.relative_to(ROOT)),
            "strongest_evidence_artifacts_for_passed": strongest_for_passed,
            "missing_obligation_remediation_map": remediation_for_missing,
            "summary": {
                "n_passed_with_strongest_artifact": len(strongest_for_passed),
                "n_missing_or_partial": len(missing),
                "all_consistent": all_consistent,
            },
            "scope_limit": "strongest-evidence packeting only; no selector-closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2176_candidate",
            "goal": "execute remediation map for any still-missing obligations and rerun P2173/P2174/P2175 chain",
        },
        "gatekeeper_checks": {
            "strongest_evidence_pack_exported": True,
            "all_consistent": all_consistent,
            "missing_obligations_present": len(missing) > 0,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2174.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2174.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2175 S1125: strict QW-2191 missing-obligation strongest-evidence pack",
                "",
                f"- Result kind: `{result_kind}`",
                f"- strongest artifacts for passed: `{len(strongest_for_passed)}`",
                f"- missing or partial obligations: `{len(missing)}`",
                f"- cross-check consistent: `{all_consistent}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
