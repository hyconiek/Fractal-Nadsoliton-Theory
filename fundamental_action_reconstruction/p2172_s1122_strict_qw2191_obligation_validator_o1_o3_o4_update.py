#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2170 = GEN / "p2170_s1120_strict_qw2191_obligation_validator_o2_o5_update.json"
IN_2171 = GEN / "p2171_s1121_strict_qw2191_o1_o3_o4_dedicated_witness_packet.json"
OUT = GEN / "p2172_s1122_strict_qw2191_obligation_validator_o1_o3_o4_update.json"
MD = GEN / "p2172_s1122_strict_qw2191_obligation_validator_o1_o3_o4_update.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2170 = load(IN_2170)
    p2171 = load(IN_2171)

    base = (p2170.get("strict_qw2191_obligation_validator_o2_o5_update", {}) or {})
    evals = (base.get("updated_obligation_evaluations", []) or [])

    w = (p2171.get("strict_qw2191_o1_o3_o4_dedicated_witness_packet", {}) or {})
    o1_pass = bool(w.get("o1_pass", False))
    o3_pass = bool(w.get("o3_pass", False))
    o4_pass = bool(w.get("o4_pass", False))

    updated = []
    for e in evals:
        ee = dict(e)
        if ee.get("id") == "O1_selector_source_identifiability" and o1_pass:
            ee["status"] = "PASS_WITH_TRACE"
            ee["pass"] = True
        if ee.get("id") == "O3_noncyclic_anchor_certificate" and o3_pass:
            ee["status"] = "PASS_WITH_TRACE"
            ee["pass"] = True
        if ee.get("id") == "O4_no_legacy_role_transfer_certificate" and o4_pass:
            ee["status"] = "PASS_WITH_TRACE"
            ee["pass"] = True
        updated.append(ee)

    required = [e for e in updated if e.get("required_for_legal_progress", True)]
    n_required = len(required)
    n_pass = sum(1 for e in required if bool(e.get("pass", False)))
    n_open = n_required - n_pass
    all_required_pass = n_open == 0 and n_required > 0

    result_kind = (
        "PASS_STRICT_QW2191_ALL_REQUIRED_OBLIGATIONS_DISCHARGED_WITH_TRACE"
        if all_required_pass
        else "OPEN_STRICT_QW2191_REMAINING_REQUIRED_OBLIGATIONS_WITH_TRACE"
    )

    payload = {
        "schema_version": "p2172_s1122_v1",
        "packet_id": "P2172",
        "stage_id": "S1122",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_obligation_validator_o1_o3_o4_update": {
            "source_validator": str(IN_2170.relative_to(ROOT)),
            "source_witness": str(IN_2171.relative_to(ROOT)),
            "updated_obligation_evaluations": updated,
            "required_summary": {
                "n_required": n_required,
                "n_pass": n_pass,
                "n_open": n_open,
                "all_required_pass": all_required_pass,
            },
            "scope_limit": "validator update over O1/O3/O4 after O2/O5 progression; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2173_candidate",
            "goal": "build independent cross-checker for O1..O5 evaluator and evidentiary trace normalization",
        },
        "gatekeeper_checks": {
            "validator_update_exported": True,
            "o1_promoted": o1_pass,
            "o3_promoted": o3_pass,
            "o4_promoted": o4_pass,
            "all_required_pass": all_required_pass,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2171.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2171.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2172 S1122: strict QW-2191 obligation validator O1/O3/O4 update",
                "",
                f"- Result kind: `{result_kind}`",
                f"- O1 promoted: `{o1_pass}`",
                f"- O3 promoted: `{o3_pass}`",
                f"- O4 promoted: `{o4_pass}`",
                f"- all_required_pass: `{all_required_pass}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
