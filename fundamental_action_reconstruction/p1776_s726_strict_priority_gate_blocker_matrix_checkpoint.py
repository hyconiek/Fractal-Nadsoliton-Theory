#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1764 = GEN / "p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json"
IN1765 = GEN / "p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json"
IN1775 = GEN / "p1775_s725_strict_w1_hr2_full_export_delivery_attempt_checkpoint.json"
OUT = GEN / "p1776_s726_strict_priority_gate_blocker_matrix_checkpoint.json"


def main() -> None:
    p1764 = json.loads(IN1764.read_text(encoding="utf-8"))
    p1765 = json.loads(IN1765.read_text(encoding="utf-8"))
    p1775 = json.loads(IN1775.read_text(encoding="utf-8"))

    ea_state = p1764["d1_d2_upgrade"]["E_A_mu_nonproxy_explicit_v2"]["componentwise_state"]
    eh_state = p1764["d1_d2_upgrade"]["E_H_nonproxy_explicit_v2"]["componentwise_state"]
    elg_state = p1765["explicit_nonproxy_EL_g_munu_v1"]["classification"]["strict_local"]
    w1_verdict = p1775["acceptance_verdict"]

    blockers = {
        "explicit_covariant_nonproxy_E_A_mu": {
            "status": "EXPORTED_OPERATOR_OPEN_COMPONENTWISE" if "OPEN" in ea_state else "READY",
            "evidence": ea_state,
        },
        "explicit_covariant_nonproxy_E_H": {
            "status": "EXPORTED_OPERATOR_OPEN_COMPONENTWISE" if "OPEN" in eh_state else "READY",
            "evidence": eh_state,
        },
        "metric_EL_g_export": {
            "status": "EXPORTED_OPERATOR_OPEN_COMPONENTWISE" if "OPEN" in elg_state else "READY",
            "evidence": elg_state,
        },
        "boundary_term_control": {
            "status": "CONTRACT_FINALIZED_REUSED",
            "evidence": p1765["consistency_gates_update"]["boundary_term_control"],
        },
        "H1_4D_weak_form_readiness": {
            "status": "OPEN_COMPONENTWISE_EXECUTION_REQUIRED",
            "evidence": p1764["readiness_update"]["H1_4D_strict_local_readiness"],
        },
        "Bianchi_Ward_consistency": {
            "status": "BLOCKED_BY_W1_NOT_FULL_EXPORT" if "OBSTRUCTION" in w1_verdict else "OPEN",
            "evidence": w1_verdict,
        },
        "BRST_Cutkosky_theorem_gates": {
            "status": "BLOCKED_BY_BW_GATE_AND_MISSING_FULL_EXPORTS",
            "evidence": p1775["gate_policy_effect"],
        },
    }

    payload = {
        "checkpoint": "P1776_S726",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchors": ["p1764_s714", "p1765_s715", "p1775_s725"],
        "priority_blocker_matrix": blockers,
        "proofs_and_nonproofs": {
            "proven_now": [
                "E_A^mu and E_H are explicitly exported at covariant operator level (nonproxy).",
                "EL_g^{mu nu} is explicitly exported at covariant operator level (nonproxy).",
                "W1 acceptance remains obstruction; GBW rerun is not admissible.",
            ],
            "still_open": [
                "Componentwise H1 residual witness on shared background family.",
                "Componentwise EL_g-E_munu residual basis closure (B1/B2/B3/C1/C2).",
                "Bianchi/Ward divergence pass-zero or obstruction trace at theorem gate quality.",
                "BRST nilpotency and Cutkosky unitarity theorem-level witnesses.",
            ],
        },
        "no_false_pass_claim": True,
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "next_honest_step": "Dowieźć brakujące współczynniki projekcji W1 (B2/B3 + dywergencja H(R2)), zamknąć W1 jako FULL_EXPORT, a następnie uruchomić komponentowy H1 i EL_g-E_munu wraz z raportem Bianchi/Ward.",
        "lay_summary": "Mamy już jawne równania, ale jeszcze nie końcowe dowody. Najpierw trzeba dokończyć jeden blok techniczny W1, żeby legalnie odblokować testy zgodności kwantowej.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
