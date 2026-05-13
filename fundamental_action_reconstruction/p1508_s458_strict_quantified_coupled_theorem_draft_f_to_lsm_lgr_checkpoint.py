#!/usr/bin/env python3
"""P1508 S4.58: export strict quantified coupled theorem draft for F=>LSM+LGR."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1500 = GEN / "p1500_s450_qw2191_selector_source_and_f_mapping_witness_summary.json"
P1507 = GEN / "p1507_s457_physical_gap_vs_release81_and_strict_f_to_lsm_lgr_theorem_draft_summary.json"
SUMMARY = GEN / "p1508_s458_strict_quantified_coupled_theorem_draft_f_to_lsm_lgr_summary.json"


def main() -> None:
    s1500 = json.loads(P1500.read_text(encoding="utf-8"))
    s1507 = json.loads(P1507.read_text(encoding="utf-8"))

    s_internal = s1500.get("objects", {}).get("S_internal_v1", {})
    fmap = s1500.get("objects", {}).get("W_Fmap_v1", {})

    assumptions = {
        "A1_selector_positive": float(s_internal.get("value", 0.0)) > 0.0,
        "A2_selector_strict_internal": bool(s_internal.get("strict_internal", False)),
        "A3_weight_normalization": abs(float(fmap.get("F_to_LSM_weight", 0.0)) + float(fmap.get("F_to_LGR_weight", 0.0)) - 1.0) <= 1e-12,
        "A4_shared_orientation": str(s_internal.get("orientation", "")) == str(fmap.get("selector_orientation", "")) and bool(str(s_internal.get("orientation", ""))),
        "A5_no_legacy_bridge": not bool(s1507.get("legacy_bridge_used", True)),
    }

    theorem_draft_ready = all(assumptions.values())

    summary = {
        "packet": "P1508",
        "status": "PASS_STRICT_QUANTIFIED_COUPLED_THEOREM_DRAFT_EXPORTED" if theorem_draft_ready else "FAIL_STRICT_QUANTIFIED_COUPLED_THEOREM_DRAFT_EXPORTED",
        "scope": "STRICT_ONLY_F_NADSOLITON_TO_LSM_PLUS_LGR",
        "theorem_draft": {
            "label": "T_QW2191_strict_coupled_draft_v1",
            "assumptions": assumptions,
            "thesis": "Under A1..A5, F admits a strict-side coupled consistency realization across LSM and LGR channels with shared selector orientation.",
            "falsifier_ready_branch": "If any of A1..A5 fails under admissible perturbation families, the draft theorem is rejected (no closure upgrade).",
        },
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "Implement P1509 contradiction-branch stress checkpoint for theorem draft rejection under admissible perturbation families.",
        "layman_explanation": "Spisaliśmy formalny szkic dowodu: jeśli kilka konkretnych warunków jest spełnionych, model łączy część cząstkową i grawitacyjną spójnie. Ale jeśli którykolwiek warunek padnie, dowód odpada.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1508] status={summary['status']} theorem_draft_ready={theorem_draft_ready}")


if __name__ == "__main__":
    main()
