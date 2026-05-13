#!/usr/bin/env python3
"""P1507 S4.57: classify physical gap vs release-8.1 and set strict theorem-draft step."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1506 = GEN / "p1506_s456_f_nadsoliton_to_lsm_plus_lgr_coupled_operator_consistency_summary.json"
SUMMARY = GEN / "p1507_s457_physical_gap_vs_release81_and_strict_f_to_lsm_lgr_theorem_draft_summary.json"


def main() -> None:
    s1506 = json.loads(P1506.read_text(encoding="utf-8"))

    coupled_pass = s1506.get("status") == "PASS_COUPLED_OPERATOR_CONSISTENCY_STRICT_ONLY"
    checks = s1506.get("checks", {})

    strong_strict_consistency = coupled_pass and all(
        bool(checks.get(k, False)) for k in [
            "selector_present_and_positive",
            "f_weight_normalization",
            "shared_selector_orientation",
        ]
    ) and (not bool(checks.get("legacy_bridge_used", True)))

    release81_like_theorem_closure_ready = False

    summary = {
        "packet": "P1507",
        "status": "PASS_PHYSICAL_GAP_CLASSIFIED_AND_THEOREM_STEP_SET" if strong_strict_consistency else "FAIL_PHYSICAL_GAP_CLASSIFIED_AND_THEOREM_STEP_SET",
        "scope": "STRICT_ONLY_NO_LEGACY_BRIDGE",
        "physical_gap_vs_release81": {
            "current_level": "strong_strict_consistency_witness",
            "release81_target_level": "theorem_level_physical_closure",
            "gap_items": [
                "quantified_coupled_theorem_not_exported_yet",
                "full_boundary_noncontradiction_chain_not_exported_yet",
                "perturbation_stability_at_theorem_level_not_exported_yet",
            ],
            "release81_like_theorem_closure_ready": release81_like_theorem_closure_ready,
        },
        "f_direction": {
            "track": "F(Nadsoliton)=>LSM+LGR",
            "next_step": "Export strict quantified coupled theorem draft with explicit premises and falsifier-ready contradiction branch.",
        },
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "Implement P1508 strict quantified coupled theorem draft checkpoint for F->LSM+LGR.",
        "layman_explanation": "Dziś wiemy, że dwie części modelu są zgodne, ale to jeszcze nie pełny dowód. Następny krok to formalne twierdzenie matematyczno-fizyczne, które da się próbować obalić.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1507] status={summary['status']} strong_strict_consistency={strong_strict_consistency}")


if __name__ == "__main__":
    main()
