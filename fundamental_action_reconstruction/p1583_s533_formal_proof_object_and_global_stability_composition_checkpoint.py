#!/usr/bin/env python3
"""P1583/S533: formal proof-object skeleton for T1582 + composition with global stability requirement."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1582 = GEN / "p1582_s532_strict_selector_uniqueness_theorem_bridge_to_full_lagrangian_summary.json"


def main() -> None:
    if not IN1582.exists():
        raise FileNotFoundError("Run P1582 before P1583")

    s1582 = json.loads(IN1582.read_text(encoding="utf-8"))
    branch_gap = float(s1582["theorem_candidate"]["branch_gap"])
    mono = float(s1582["theorem_candidate"]["monotonic_energy_ratio"])

    lemma_L1 = branch_gap < 0.08
    lemma_L2 = mono >= 0.55
    lemma_L3 = False  # still missing theorem-level global stability object from P1571-P1577 sequence

    proof_object = {
        "theorem": "T1583_formal_composition_candidate_for_strict_core_progress",
        "premises": [
            "L1_selector_branch_gap_bound",
            "L2_monotonic_selector_energy_segment",
            "L3_global_stability_link_to_SM_GR_bundle",
        ],
        "checks": {
            "L1_selector_branch_gap_bound": lemma_L1,
            "L2_monotonic_selector_energy_segment": lemma_L2,
            "L3_global_stability_link_to_SM_GR_bundle": lemma_L3,
        },
        "derivation_note": "Strict chain is preserved from K_strict to EOM; closure blocked by unmet L1/L3.",
    }

    all_ok = all(proof_object["checks"].values())
    status = "PASS_T1583_FORMAL_PROOF_OBJECT_CANDIDATE" if all_ok else "FAIL_T1583_FORMAL_PROOF_OBJECT_CANDIDATE"

    summary = {
        "checkpoint": "P1583_S533",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": s1582["strict_chain"],
        "proof_object": proof_object,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_core_closure": {
            "status": "OPEN",
            "qw2191_closed": False,
            "remaining_missing": [
                "tight_selector_branch_gap_witness (L1)",
                "global_stability_theorem_object (L3)",
                "final_ToE_composition_theorem",
            ],
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1584_construct_L1_witness_refinement_and_L3_global_stability_object",
        "lay_summary": "Mamy szkic formalnego dowodu, ale dwa warunki krytyczne (L1 i L3) nie są jeszcze spełnione, więc domknięcie strict-core pozostaje otwarte."
    }

    out = GEN / "p1583_s533_formal_proof_object_and_global_stability_composition_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
