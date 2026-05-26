#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2175 = GEN / "p2175_s1125_strict_qw2191_missing_obligation_strongest_evidence_pack.json"
IN_2174 = GEN / "p2174_s1124_strict_qw2191_evidentiary_trace_normalizer_and_missing_obligation_witness_pack.json"
OUT = GEN / "p2176_s1126_strict_qw2191_remediation_execution_and_chain_rerun.json"
MD = GEN / "p2176_s1126_strict_qw2191_remediation_execution_and_chain_rerun.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def build_remediation_artifact(rem: dict[str, Any]) -> dict[str, Any]:
    oid = rem.get("obligation_id", "unknown")
    return {
        "obligation_id": oid,
        "artifact_id": rem.get("required_next_artifact", f"EVID_{oid}_STRONGEST_V1"),
        "execution_status": "EXECUTED_WITH_TRACE",
        "notes": "Generated remediation artifact skeleton for machine-traceable follow-up; no closure claim.",
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2175 = load(IN_2175)
    p2174 = load(IN_2174)

    pack_2175 = p2175.get("strict_qw2191_missing_obligation_strongest_evidence_pack", {}) or {}
    rem_map = pack_2175.get("missing_obligation_remediation_map", []) or []

    executed = [build_remediation_artifact(r) for r in rem_map]

    # chain rerun summary references prior chain statuses instead of claiming new theorem closure
    prior_consistent = bool((p2175.get("gatekeeper_checks", {}) or {}).get("all_consistent", False))
    missing_present_after_pack = bool((p2175.get("gatekeeper_checks", {}) or {}).get("missing_obligations_present", False))

    # If there were missing obligations, this step executes remediation skeletons;
    # still open unless no missing obligations remain.
    remaining_missing = len(rem_map)
    all_resolved = remaining_missing == 0

    result_kind = (
        "PASS_STRICT_QW2191_REMEDIATION_CHAIN_RERUN_READY_WITH_TRACE"
        if all_resolved and prior_consistent
        else "OPEN_STRICT_QW2191_REMEDIATION_CHAIN_PARTIAL_WITH_TRACE"
    )

    payload = {
        "schema_version": "p2176_s1126_v1",
        "packet_id": "P2176",
        "stage_id": "S1126",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_remediation_execution_and_chain_rerun": {
            "source_strongest_pack": str(IN_2175.relative_to(ROOT)),
            "source_trace_normalizer": str(IN_2174.relative_to(ROOT)),
            "executed_remediation_artifacts": executed,
            "chain_rerun_summary": {
                "prior_chain_consistent": prior_consistent,
                "missing_present_after_p2175": missing_present_after_pack,
                "remaining_missing_count": remaining_missing,
                "all_resolved": all_resolved,
            },
            "scope_limit": "remediation execution and chain-rerun readiness only; no selector-closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2177_candidate",
            "goal": "convert remediation skeletons into obligation-specific strongest evidence artifacts and re-evaluate missing count",
        },
        "gatekeeper_checks": {
            "remediation_execution_exported": True,
            "prior_chain_consistent": prior_consistent,
            "all_resolved": all_resolved,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2175.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2175.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2176 S1126: strict QW-2191 remediation execution + chain rerun",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Executed remediation artifacts: `{len(executed)}`",
                f"- Remaining missing count: `{remaining_missing}`",
                f"- Prior chain consistent: `{prior_consistent}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
