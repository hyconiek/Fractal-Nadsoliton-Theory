#!/usr/bin/env python3
"""P1818 S768 strict verdict->state-vector consistency lock checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"


def load(name: str) -> dict:
    path = GENERATED / name
    if not path.exists():
        return {"_missing": name}
    return json.loads(path.read_text(encoding="utf-8"))


def normalize_sv_status(status: str) -> str:
    if status.startswith("PASS"):
        return "PASS"
    if status.startswith("OPEN"):
        return "OPEN"
    return "OTHER"


def main() -> None:
    p1793 = load("p1793_s743_sv1_sv5_evidence_intake_validation_result.json")
    p1794 = load("p1794_s744_strict_sv1_to_sv5_verdict_to_state_update_protocol_checkpoint.json")
    p1817 = load("p1817_s767_strict_intake_verdict_hardening_checkpoint.json")

    verdict = p1793.get("verdict", "OPEN_OBSTRUCTION_WITH_TRACE")
    sv = p1794.get("updated_state_vector", {})

    sv_keys = [
        "SV1_E_A_mu_nonproxy_covariant_explicit",
        "SV2_E_H_nonproxy_covariant_explicit",
        "SV3_EL_g_nonproxy_explicit",
        "SV4_boundary_term_control",
        "SV5_H1_4D_weak_form",
    ]

    violations = []
    allowed_prefix = "PASS" if verdict == "PASS_ZERO" else "OPEN"

    for key in sv_keys:
        status = sv.get(key, "MISSING")
        lane = normalize_sv_status(status)
        if lane == "OTHER" or (allowed_prefix == "OPEN" and lane == "PASS"):
            violations.append({
                "state_vector_key": key,
                "status": status,
                "reason": "STATE_VECTOR_STATUS_STRONGER_THAN_SOURCE_VERDICT",
            })

    out_status = "LOCK_CONSISTENT" if not violations else "OPEN_OBSTRUCTION_WITH_TRACE"

    out = {
        "packet_id": "P1818",
        "stage_id": "S768",
        "status": out_status,
        "route": "strict_only",
        "lock_name": "STRICT_VERDICT_TO_STATEVECTOR_CONSISTENCY_LOCK",
        "source_verdict": verdict,
        "hardened_verdict": p1817.get("post_fix_validation", {}).get("p1793_verdict", verdict),
        "checked_state_vector_keys": sv_keys,
        "consistency_rule": "SV1_to_SV5_must_not_exceed_source_verdict_strength",
        "violations": violations,
        "technical_progress": "Executable lock added to prevent false-pass drift between intake verdict and state-vector SV1..SV5.",
        "proven": "When source verdict is OPEN, state-vector cannot report PASS across SV1..SV5.",
        "open": "Real P1806 evidence injection and rerun chain still required for any theorem-gate PASS_ZERO promotion.",
        "false_pass_risk": "Manual edits or stale generated artifacts could bypass intent without this machine lock.",
        "next_honest_step": "Bind this lock as mandatory pre-gate check before BW/BRST/Cutkosky theorem-gate evaluation.",
        "lay_explanation": "To nie pozwala systemowi 'udawać sukcesu': jeśli źródłowy werdykt jest OPEN, wektor stanu też musi być OPEN.",
    }

    out_path = GENERATED / "p1818_s768_strict_verdict_statevector_consistency_lock_checkpoint.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
