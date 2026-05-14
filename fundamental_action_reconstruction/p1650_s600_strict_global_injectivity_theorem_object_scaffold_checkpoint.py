#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1649 = GEN / "p1649_s599_strict_regional_projector_injectivity_witness_summary.json"


def main() -> None:
    s49 = json.loads(IN1649.read_text(encoding="utf-8"))

    theorem_object = {
        "id": "T_STRICT_INJ_BOX_01",
        "goal": "Sufficient conditions for injectivity of projector map {I0,I2,I_mix,I_curv} on parameter box B.",
        "box_B": s49["region"],
        "assumptions": [
            "A1: projector family is C1 on B and extends continuously to boundary",
            "A2: det(J(theta)) keeps constant sign on B",
            "A3: overlap patch maps between local charts are C1 and orientation-preserving",
            "A4: selector branch fixed by strict internal witness class (no legacy bridge)",
        ],
        "obligations": [
            "O1: export analytic lower bound det(J)>=eps>0 on B",
            "O2: export overlap cocycle consistency theorem",
            "O3: export nullspace exclusion witness compatible with QW-2191",
        ],
        "conclusion_template": "If A1..A4 and O1..O3 hold, then coefficients->kernel map is injective on B up to admitted gauge-equivalence.",
        "status": "OPEN",
    }

    summary = {
        "checkpoint": "P1650_S600",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1650_GLOBAL_INJECTIVITY_THEOREM_OBJECT_SCAFFOLD_EXPORTED",
        "strict_only": True,
        "legacy_bridge_used": False,
        "kernel_definition": s49["kernel_definition"],
        "full_chain": "K_strict -> coefficients -> L_total -> EOM and reverse obligations",
        "regional_evidence_reference": {
            "status": s49["status"],
            "sample_count": s49["sample_count"],
            "min_abs_detJ": s49["min_abs_detJ"],
        },
        "global_injectivity_theorem_object": theorem_object,
        "strict_core_closure": {
            "status": "OPEN",
            "reason": "Theorem object exported but obligations O1..O3 not yet discharged",
        },
        "next_honest_step": "Udowodnić O1: jawna dolna granica det(J)>=eps na całym B i dołączyć certyfikat obliczeniowy.",
        "lay_summary": "Mamy teraz szkic formalnego twierdzenia, które tłumaczy co dokładnie trzeba jeszcze pokazać, by przejść z mocnych testów lokalnych do pełnego dowodu globalnego.",
    }

    out = GEN / "p1650_s600_strict_global_injectivity_theorem_object_scaffold_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
