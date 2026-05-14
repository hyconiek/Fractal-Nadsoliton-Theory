#!/usr/bin/env python3
"""P1585/S535: promote L1/L3 candidates to formal theorem-object candidates and compose ToE closure gate."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1584 = GEN / "p1584_s534_l1_witness_refinement_and_l3_global_stability_object_summary.json"


def main() -> None:
    if not IN1584.exists():
        raise FileNotFoundError("Run P1584 before P1585")

    s = json.loads(IN1584.read_text(encoding="utf-8"))
    l1_ok = bool(s["l1_refinement"]["pass"])
    l3_ok = bool(s["l3_global_stability_object_candidate"]["pass"])

    theorem_objects = {
        "T1585_L1_formal_object_candidate": {
            "premise": "trimmed antisymmetry witness + strict selector source continuity",
            "validated": l1_ok,
        },
        "T1585_L3_formal_object_candidate": {
            "premise": "coercivity + bounded curvature imply stability-window object",
            "validated": l3_ok,
        },
    }

    toe_composition_gate = l1_ok and l3_ok
    status = "PASS_P1585_FORMAL_OBJECTS_AND_TOE_COMPOSITION_GATE" if toe_composition_gate else "FAIL_P1585_FORMAL_OBJECTS_AND_TOE_COMPOSITION_GATE"

    summary = {
        "checkpoint": "P1585_S535",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": s["strict_chain"],
        "theorem_objects": theorem_objects,
        "toe_composition_gate": {
            "ready": toe_composition_gate,
            "required": [
                "T1585_L1_formal_object_candidate.validated",
                "T1585_L3_formal_object_candidate.validated",
            ],
        },
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_core_closure": {
            "status": "OPEN",
            "qw2191_closed": False,
            "remaining_missing": [
                "untrimmed_domain_L1_formalization",
                "full_global_proof_for_L3_beyond_window",
                "final_ToE_composition_theorem_export",
            ],
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1586_export_final_ToE_composition_theorem_or_nonclosure_certificate",
        "lay_summary": "Mamy formalne obiekty-kandydaty L1 i L3, ale brakuje pełnej wersji globalnej, więc strict-core closure nadal pozostaje otwarty."
    }

    out = GEN / "p1585_s535_formal_l1_l3_theorem_objects_and_toe_composition_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
