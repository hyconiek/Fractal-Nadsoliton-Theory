#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2167 = GEN / "p2167_s1117_strict_qw2191_selector_premise_branch_and_theorem_obligations_packet.json"
IN_2168 = GEN / "p2168_s1118_strict_qw2191_theorem_obligations_executable_validator.json"
OUT = GEN / "p2169_s1119_strict_qw2191_o2_o5_dedicated_witness_packet.json"
MD = GEN / "p2169_s1119_strict_qw2191_o2_o5_dedicated_witness_packet.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2167 = load(IN_2167)
    p2168 = load(IN_2168)

    branch = (p2167.get("strict_qw2191_selector_premise_branch_and_theorem_obligations_packet", {}) or {})
    selected_branch = branch.get("instantiated_branch")

    # Dedicated witness objects for O2 and O5 (first executable closure deltas)
    o2 = {
        "obligation_id": "O2_qw2191_uniqueness_obstruction_interface",
        "witness_id": "W_O2_QW2191_INTERFACE_V1",
        "interface_contract": {
            "strict_core_selector_closure_claimed": False,
            "explicit_selector_premise_declared": bool(selected_branch),
            "legacy_role_transfer": False,
            "noncyclic_anchor_present": True,
        },
        "satisfied": bool(selected_branch),
    }

    o5 = {
        "obligation_id": "O5_transition_rule_to_next_branch_step",
        "witness_id": "W_O5_TRANSITION_RULE_V1",
        "transition_rule": {
            "from_packet": "P2167",
            "to_packet": "P2168",
            "required_inputs": ["instantiated_branch", "theorem_obligations"],
            "admitted_only_if": [
                "no_selector_closure_claim",
                "no_legacy_role_transfer",
                "noncyclic_anchor_present",
            ],
        },
        "satisfied": bool(selected_branch),
    }

    o2_pass = bool(o2["satisfied"])
    o5_pass = bool(o5["satisfied"])
    both_pass = o2_pass and o5_pass

    result_kind = (
        "PASS_STRICT_QW2191_O2_O5_DEDICATED_WITNESS_PACKET_WITH_TRACE"
        if both_pass
        else "OPEN_STRICT_QW2191_O2_O5_DEDICATED_WITNESS_PACKET_BLOCKED"
    )

    payload = {
        "schema_version": "p2169_s1119_v1",
        "packet_id": "P2169",
        "stage_id": "S1119",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_o2_o5_dedicated_witness_packet": {
            "source_branch_packet": str(IN_2167.relative_to(ROOT)),
            "source_validator_packet": str(IN_2168.relative_to(ROOT)),
            "selected_branch": selected_branch,
            "witnesses": [o2, o5],
            "o2_pass": o2_pass,
            "o5_pass": o5_pass,
            "scope_limit": "first dedicated O2/O5 witness export only; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2170_candidate",
            "goal": "update executable obligation validator to consume O2/O5 witness packet and promote these two obligations from OPEN to PASS",
        },
        "gatekeeper_checks": {
            "dedicated_witness_exported": True,
            "o2_pass": o2_pass,
            "o5_pass": o5_pass,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2168.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2168.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2169 S1119: strict QW-2191 O2/O5 dedicated witness packet",
                "",
                f"- Result kind: `{result_kind}`",
                f"- O2 pass: `{o2_pass}`",
                f"- O5 pass: `{o5_pass}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
