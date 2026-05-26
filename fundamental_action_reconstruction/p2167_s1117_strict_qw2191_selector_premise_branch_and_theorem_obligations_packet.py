#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2166 = GEN / "p2166_s1116_strict_qw2191_noncyclic_selector_source_witness.json"
OUT = GEN / "p2167_s1117_strict_qw2191_selector_premise_branch_and_theorem_obligations_packet.json"
MD = GEN / "p2167_s1117_strict_qw2191_selector_premise_branch_and_theorem_obligations_packet.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2166 = load(IN_2166)
    w = (p2166.get("strict_qw2191_noncyclic_selector_source_witness", {}) or {})

    selected_branch = w.get("selected_admissible_branch")
    branch_instantiated = bool(selected_branch)

    theorem_obligations = [
        {
            "id": "O1_selector_source_identifiability",
            "description": "Prove explicit strict-lane selector-source map is well-defined and non-empty on declared domain.",
            "status": "OPEN",
            "required_for_legal_progress": True,
        },
        {
            "id": "O2_qw2191_uniqueness_obstruction_interface",
            "description": "Show branch provides admissible interface to QW-2191 without claiming strict-core closure.",
            "status": "OPEN",
            "required_for_legal_progress": True,
        },
        {
            "id": "O3_noncyclic_anchor_certificate",
            "description": "Export certificate that route uses noncyclic anchor/provider class and does not repeat L5/L12 same blocker-cut.",
            "status": "OPEN",
            "required_for_legal_progress": True,
        },
        {
            "id": "O4_no_legacy_role_transfer_certificate",
            "description": "Certify no legacy-role transfer into strict selector branch assumptions.",
            "status": "OPEN",
            "required_for_legal_progress": True,
        },
        {
            "id": "O5_transition_rule_to_next_branch_step",
            "description": "Define legal transition rule from premise branch to executable selector witness validator packet.",
            "status": "OPEN",
            "required_for_legal_progress": True,
        },
    ]

    open_obligations = [o["id"] for o in theorem_obligations if o["status"] != "CLOSED"]
    legal_progress_ready = branch_instantiated and len(open_obligations) > 0

    result_kind = (
        "PASS_STRICT_QW2191_SELECTOR_PREMISE_BRANCH_AND_THEOREM_OBLIGATIONS_PACKET_WITH_TRACE"
        if branch_instantiated
        else "OPEN_STRICT_QW2191_SELECTOR_PREMISE_BRANCH_AND_THEOREM_OBLIGATIONS_PACKET_BLOCKED"
    )

    payload = {
        "schema_version": "p2167_s1117_v1",
        "packet_id": "P2167",
        "stage_id": "S1117",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_selector_premise_branch_and_theorem_obligations_packet": {
            "source_witness": str(IN_2166.relative_to(ROOT)),
            "instantiated_branch": selected_branch,
            "branch_instantiated": branch_instantiated,
            "theorem_obligations": theorem_obligations,
            "open_obligations": open_obligations,
            "legal_progress_ready": legal_progress_ready,
            "scope_limit": "branch/theorem-obligation mapping only; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2168_candidate",
            "goal": "build executable obligation validator for O1-O5 and export first closure deltas against QW-2191",
        },
        "gatekeeper_checks": {
            "branch_packet_exported": True,
            "branch_instantiated": branch_instantiated,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2166.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2166.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2167 S1117: strict QW-2191 selector-premise branch and theorem-obligations packet",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Instantiated branch: `{selected_branch}`",
                f"- Open obligations: `{open_obligations}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
