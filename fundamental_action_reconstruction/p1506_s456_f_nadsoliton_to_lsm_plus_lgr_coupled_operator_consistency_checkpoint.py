#!/usr/bin/env python3
"""P1506 S4.56: coupled operator-consistency check for F(Nadsoliton)=>LSM+LGR."""

from __future__ import annotations

import json
from math import isclose
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1500 = GEN / "p1500_s450_qw2191_selector_source_and_f_mapping_witness_summary.json"
P1505 = GEN / "p1505_s455_internal_selector_source_and_release81_comparison_to_f_lsm_lgr_summary.json"
SUMMARY = GEN / "p1506_s456_f_nadsoliton_to_lsm_plus_lgr_coupled_operator_consistency_summary.json"


def main() -> None:
    s1500 = json.loads(P1500.read_text(encoding="utf-8"))
    s1505 = json.loads(P1505.read_text(encoding="utf-8"))

    s_internal = s1500.get("objects", {}).get("S_internal_v1", {})
    fmap = s1500.get("objects", {}).get("W_Fmap_v1", {})

    selector_present = bool(s_internal.get("strict_internal", False)) and float(s_internal.get("value", 0.0)) > 0.0
    w_lsm = float(fmap.get("F_to_LSM_weight", 0.0))
    w_lgr = float(fmap.get("F_to_LGR_weight", 0.0))
    normalized = isclose(w_lsm + w_lgr, 1.0, rel_tol=0.0, abs_tol=1e-12)

    selector_orientation = str(s_internal.get("orientation", ""))
    fmap_orientation = str(fmap.get("selector_orientation", ""))
    shared_orientation = bool(selector_orientation) and selector_orientation == fmap_orientation

    legacy_bridge_used = bool(s1505.get("legacy_bridge_used", False))
    coupled_consistency_pass = selector_present and normalized and shared_orientation and (not legacy_bridge_used)

    summary = {
        "packet": "P1506",
        "status": "PASS_COUPLED_OPERATOR_CONSISTENCY_STRICT_ONLY" if coupled_consistency_pass else "FAIL_COUPLED_OPERATOR_CONSISTENCY_STRICT_ONLY",
        "scope": "STRICT_ONLY_F_NADSOLITON_TO_LSM_PLUS_LGR",
        "checks": {
            "selector_present_and_positive": selector_present,
            "f_weight_normalization": normalized,
            "shared_selector_orientation": shared_orientation,
            "legacy_bridge_used": legacy_bridge_used,
        },
        "computed": {
            "F_to_LSM_weight": w_lsm,
            "F_to_LGR_weight": w_lgr,
            "selector_orientation": selector_orientation,
            "fmap_orientation": fmap_orientation,
        },
        "qw2191_closed": False,
        "next_honest_step": "Promote to strict-side quantified coupled theorem draft tying F->LSM and F->LGR under one selector orientation, still without legacy bridge transfer.",
        "layman_explanation": "Sprawdziliśmy, czy model mówi jednym głosem: część cząstkowa i grawitacyjna mają wspólne ustawienie kierunku oraz poprawne proporcje. To ważny krok do domknięcia, ale jeszcze nie finał.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1506] status={summary['status']} coupled_consistency_pass={coupled_consistency_pass}")


if __name__ == "__main__":
    main()
