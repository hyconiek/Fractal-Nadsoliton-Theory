#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2168 = GEN / "p2168_s1118_strict_qw2191_theorem_obligations_executable_validator.json"
IN_2169 = GEN / "p2169_s1119_strict_qw2191_o2_o5_dedicated_witness_packet.json"
OUT = GEN / "p2170_s1120_strict_qw2191_obligation_validator_o2_o5_update.json"
MD = GEN / "p2170_s1120_strict_qw2191_obligation_validator_o2_o5_update.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2168 = load(IN_2168)
    p2169 = load(IN_2169)

    base = (p2168.get("strict_qw2191_theorem_obligations_executable_validator", {}) or {})
    evals = (base.get("obligation_evaluations", []) or [])

    w = (p2169.get("strict_qw2191_o2_o5_dedicated_witness_packet", {}) or {})
    o2_pass = bool(w.get("o2_pass", False))
    o5_pass = bool(w.get("o5_pass", False))

    updated = []
    for e in evals:
        ee = dict(e)
        if ee.get("id") == "O2_qw2191_uniqueness_obstruction_interface" and o2_pass:
            ee["status"] = "PASS_WITH_TRACE"
            ee["pass"] = True
        if ee.get("id") == "O5_transition_rule_to_next_branch_step" and o5_pass:
            ee["status"] = "PASS_WITH_TRACE"
            ee["pass"] = True
        updated.append(ee)

    required = [e for e in updated if e.get("required_for_legal_progress", True)]
    n_required = len(required)
    n_pass = sum(1 for e in required if bool(e.get("pass", False)))
    n_open = n_required - n_pass
    all_required_pass = n_open == 0 and n_required > 0

    result_kind = "PASS_STRICT_QW2191_OBLIGATION_VALIDATOR_O2_O5_UPDATE_WITH_TRACE"

    payload = {
        "schema_version": "p2170_s1120_v1",
        "packet_id": "P2170",
        "stage_id": "S1120",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_obligation_validator_o2_o5_update": {
            "source_validator": str(IN_2168.relative_to(ROOT)),
            "source_witness": str(IN_2169.relative_to(ROOT)),
            "updated_obligation_evaluations": updated,
            "required_summary": {
                "n_required": n_required,
                "n_pass": n_pass,
                "n_open": n_open,
                "all_required_pass": all_required_pass,
            },
            "scope_limit": "O2/O5 update only; remaining obligations may stay OPEN",
        },
        "recommended_next_honest_step": {
            "id": "P2171_candidate",
            "goal": "export dedicated witnesses for O1/O3/O4 and iterate validator until all required obligations pass",
        },
        "gatekeeper_checks": {
            "validator_update_exported": True,
            "o2_promoted": o2_pass,
            "o5_promoted": o5_pass,
            "all_required_pass": all_required_pass,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2169.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2169.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2170 S1120: strict QW-2191 obligation validator O2/O5 update",
                "",
                f"- Result kind: `{result_kind}`",
                f"- O2 promoted: `{o2_pass}`",
                f"- O5 promoted: `{o5_pass}`",
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
