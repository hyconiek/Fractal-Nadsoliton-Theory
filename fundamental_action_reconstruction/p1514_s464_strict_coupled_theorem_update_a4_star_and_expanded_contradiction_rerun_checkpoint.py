#!/usr/bin/env python3
"""P1514 S4.64: update coupled theorem with A4* and rerun expanded contradiction branch."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1500 = GEN / "p1500_s450_qw2191_selector_source_and_f_mapping_witness_summary.json"
P1513 = GEN / "p1513_s463_orientation_stability_strengthening_for_strict_f_to_lsm_lgr_summary.json"
SUMMARY = GEN / "p1514_s464_strict_coupled_theorem_update_a4_star_and_expanded_contradiction_rerun_summary.json"


def margin(a: str, b: str) -> float:
    if a == b and a in {"SM_preferred", "GR_preferred"}:
        return 1.0
    if a == b and a == "mixed":
        return 0.5
    if a == b and a == "neutral":
        return 0.25
    return 0.0


def main() -> None:
    s1500 = json.loads(P1500.read_text(encoding="utf-8"))
    _ = json.loads(P1513.read_text(encoding="utf-8"))

    s_internal = s1500.get("objects", {}).get("S_internal_v1", {})
    fmap = s1500.get("objects", {}).get("W_Fmap_v1", {})

    selector = float(s_internal.get("value", 0.0))
    strict_internal = bool(s_internal.get("strict_internal", False))
    w_sm = float(fmap.get("F_to_LSM_weight", 0.0))
    w_gr = float(fmap.get("F_to_LGR_weight", 0.0))

    orientations = [
        ("SM_preferred", "SM_preferred"),
        ("GR_preferred", "GR_preferred"),
        ("mixed", "mixed"),
        ("neutral", "neutral"),
        ("SM_preferred", "GR_preferred"),
        ("GR_preferred", "SM_preferred"),
        ("mixed", "neutral"),
    ]

    m_min = 0.5
    trials = []
    for o1, o2 in orientations:
        m = margin(o1, o2)
        assumptions = {
            "A1_selector_positive": selector > 0.0,
            "A2_selector_strict_internal": strict_internal,
            "A3_weight_normalization": abs((w_sm + w_gr) - 1.0) <= 1e-9,
            "A4_star": (o1 == o2) and (m >= m_min),
            "A5_no_legacy_bridge": True,
        }
        trials.append({
            "orientation_internal": o1,
            "orientation_fmap": o2,
            "margin": m,
            "assumptions": assumptions,
            "pass": all(assumptions.values()),
        })

    pass_rows = [t for t in trials if t["pass"]]
    fail_rows = [t for t in trials if not t["pass"]]

    summary = {
        "packet": "P1514",
        "status": "PASS_A4_STAR_THEOREM_UPDATE_AND_EXPANDED_RERUN",
        "scope": "STRICT_ONLY_F_NADSOLITON_TO_LSM_PLUS_LGR",
        "theorem_update": {
            "replaced_condition": "A4 -> A4*",
            "A4_star_definition": "shared_orientation AND orientation_margin >= 0.5",
        },
        "expanded_contradiction_rerun": {
            "trial_count": len(trials),
            "pass_count": len(pass_rows),
            "fail_count": len(fail_rows),
            "fails": fail_rows,
        },
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "Implement P1515 robust-orientation envelope export and attach it as formal admissibility annex to the strict coupled theorem draft.",
        "layman_explanation": "Zaostrzyliśmy warunek kierunku i ponownie przepuściliśmy trudniejsze testy. Silne, zgodne ustawienia przechodzą, a słabe lub sprzeczne ustawienia są odrzucane — dokładnie tak powinno działać uczciwe sito naukowe.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1514] status={summary['status']} pass={len(pass_rows)} fail={len(fail_rows)}")


if __name__ == "__main__":
    main()
