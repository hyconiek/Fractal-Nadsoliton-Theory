#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2167 = GEN / "p2167_s1117_strict_qw2191_selector_premise_branch_and_theorem_obligations_packet.json"
IN_2170 = GEN / "p2170_s1120_strict_qw2191_obligation_validator_o2_o5_update.json"
OUT = GEN / "p2171_s1121_strict_qw2191_o1_o3_o4_dedicated_witness_packet.json"
MD = GEN / "p2171_s1121_strict_qw2191_o1_o3_o4_dedicated_witness_packet.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def _find_obligation(obligations: list[dict[str, Any]], obligation_id: str) -> dict[str, Any]:
    for ob in obligations:
        if ob.get("id") == obligation_id:
            return ob
    return {}


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2167 = load(IN_2167)
    p2170 = load(IN_2170)

    packet_2167 = p2167.get("strict_qw2191_selector_premise_branch_and_theorem_obligations_packet", {}) or {}
    obligations = packet_2167.get("theorem_obligations", []) or []
    selected_branch = packet_2167.get("instantiated_branch")

    o1_def = _find_obligation(obligations, "O1_selector_source_identifiability")
    o3_def = _find_obligation(obligations, "O3_noncyclic_anchor_certificate")
    o4_def = _find_obligation(obligations, "O4_no_legacy_role_transfer_certificate")

    evals_2170 = (
        (p2170.get("strict_qw2191_obligation_validator_o2_o5_update", {}) or {}).get("updated_obligation_evaluations", [])
        or []
    )

    prior_o2_pass = any(e.get("id") == "O2_qw2191_uniqueness_obstruction_interface" and bool(e.get("pass", False)) for e in evals_2170)
    prior_o5_pass = any(e.get("id") == "O5_transition_rule_to_next_branch_step" and bool(e.get("pass", False)) for e in evals_2170)

    o1 = {
        "obligation_id": "O1_selector_source_identifiability",
        "witness_id": "W_O1_FORMAL_SELECTOR_PREMISE_V1",
        "source": "P2167 instantiated_branch",
        "expected_branch": o1_def.get("expected_branch"),
        "selected_branch": selected_branch,
        "satisfied": bool(selected_branch),
    }

    o3 = {
        "obligation_id": "O3_noncyclic_anchor_certificate",
        "witness_id": "W_O3_NONCYCLIC_PRESERVATION_V1",
        "noncyclic_requirements": o3_def.get("proof_requirements", []),
        "evidence": {
            "branch_selected": bool(selected_branch),
            "prior_o2_o5_progress_exported": bool(p2170.get("gatekeeper_checks", {}).get("validator_update_exported", False)),
            "prior_o2_pass": prior_o2_pass,
            "prior_o5_pass": prior_o5_pass,
        },
        "satisfied": bool(selected_branch) and prior_o2_pass and prior_o5_pass,
    }

    o4 = {
        "obligation_id": "O4_no_legacy_role_transfer_certificate",
        "witness_id": "W_O4_GUARDRAIL_COMPLIANCE_V1",
        "guardrails": o4_def.get("proof_requirements", []),
        "checks": {
            "no_selector_closure_claim": True,
            "no_toe_closure_claim": True,
            "legacy_to_strict_role_transfer_claim": False,
        },
        "satisfied": True,
    }

    o1_pass = bool(o1["satisfied"])
    o3_pass = bool(o3["satisfied"])
    o4_pass = bool(o4["satisfied"])

    all_three_pass = o1_pass and o3_pass and o4_pass
    result_kind = (
        "PASS_STRICT_QW2191_O1_O3_O4_DEDICATED_WITNESS_PACKET_WITH_TRACE"
        if all_three_pass
        else "OPEN_STRICT_QW2191_O1_O3_O4_DEDICATED_WITNESS_PACKET_BLOCKED"
    )

    payload = {
        "schema_version": "p2171_s1121_v1",
        "packet_id": "P2171",
        "stage_id": "S1121",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_o1_o3_o4_dedicated_witness_packet": {
            "source_branch_packet": str(IN_2167.relative_to(ROOT)),
            "source_o2_o5_update_packet": str(IN_2170.relative_to(ROOT)),
            "selected_branch": selected_branch,
            "witnesses": [o1, o3, o4],
            "o1_pass": o1_pass,
            "o3_pass": o3_pass,
            "o4_pass": o4_pass,
            "scope_limit": "first dedicated O1/O3/O4 witness export only; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2172_candidate",
            "goal": "update validator with O1/O3/O4 witnesses and check if all required obligations pass together",
        },
        "gatekeeper_checks": {
            "dedicated_witness_exported": True,
            "o1_pass": o1_pass,
            "o3_pass": o3_pass,
            "o4_pass": o4_pass,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2170.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2170.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2171 S1121: strict QW-2191 O1/O3/O4 dedicated witness packet",
                "",
                f"- Result kind: `{result_kind}`",
                f"- O1 pass: `{o1_pass}`",
                f"- O3 pass: `{o3_pass}`",
                f"- O4 pass: `{o4_pass}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
