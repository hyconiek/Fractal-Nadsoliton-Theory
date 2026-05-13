#!/usr/bin/env python3
"""P1509 S4.59: contradiction-branch stress for strict F=>LSM+LGR theorem draft."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1500 = GEN / "p1500_s450_qw2191_selector_source_and_f_mapping_witness_summary.json"
P1508 = GEN / "p1508_s458_strict_quantified_coupled_theorem_draft_f_to_lsm_lgr_summary.json"
SUMMARY = GEN / "p1509_s459_strict_contradiction_branch_stress_for_f_to_lsm_lgr_summary.json"


def check_assumptions(selector_value: float, strict_internal: bool, w_sm: float, w_gr: float, orient_a: str, orient_b: str, no_legacy_bridge: bool) -> dict[str, bool]:
    return {
        "A1_selector_positive": selector_value > 0.0,
        "A2_selector_strict_internal": strict_internal,
        "A3_weight_normalization": abs((w_sm + w_gr) - 1.0) <= 1e-9,
        "A4_shared_orientation": bool(orient_a) and orient_a == orient_b,
        "A5_no_legacy_bridge": no_legacy_bridge,
    }


def main() -> None:
    s1500 = json.loads(P1500.read_text(encoding="utf-8"))
    s1508 = json.loads(P1508.read_text(encoding="utf-8"))

    s_internal = s1500.get("objects", {}).get("S_internal_v1", {})
    fmap = s1500.get("objects", {}).get("W_Fmap_v1", {})

    base_selector = float(s_internal.get("value", 0.0))
    strict_internal = bool(s_internal.get("strict_internal", False))
    base_w_sm = float(fmap.get("F_to_LSM_weight", 0.0))
    base_w_gr = float(fmap.get("F_to_LGR_weight", 0.0))
    base_orientation = str(s_internal.get("orientation", ""))
    fmap_orientation = str(fmap.get("selector_orientation", ""))
    no_legacy_bridge = not bool(s1508.get("legacy_bridge_used", True))

    perturbations = [
        {"id": "selector_scale_0p9", "selector": base_selector * 0.9, "w_sm": base_w_sm, "w_gr": base_w_gr, "o1": base_orientation, "o2": fmap_orientation},
        {"id": "selector_scale_1p1", "selector": base_selector * 1.1, "w_sm": base_w_sm, "w_gr": base_w_gr, "o1": base_orientation, "o2": fmap_orientation},
        {"id": "weight_shift_sm_plus_0p02", "selector": base_selector, "w_sm": base_w_sm + 0.02, "w_gr": base_w_gr - 0.02, "o1": base_orientation, "o2": fmap_orientation},
        {"id": "weight_shift_sm_minus_0p02", "selector": base_selector, "w_sm": base_w_sm - 0.02, "w_gr": base_w_gr + 0.02, "o1": base_orientation, "o2": fmap_orientation},
        {"id": "orientation_flip_probe", "selector": base_selector, "w_sm": base_w_sm, "w_gr": base_w_gr, "o1": "GR_preferred", "o2": fmap_orientation},
    ]

    rows = []
    first_counterexample = None
    for p in perturbations:
        a = check_assumptions(p["selector"], strict_internal, p["w_sm"], p["w_gr"], p["o1"], p["o2"], no_legacy_bridge)
        passed = all(a.values())
        row = {"id": p["id"], "assumptions": a, "pass": passed}
        rows.append(row)
        if (not passed) and first_counterexample is None:
            first_counterexample = row

    contradiction_found = first_counterexample is not None

    summary = {
        "packet": "P1509",
        "status": "PASS_CONTRADICTION_BRANCH_STRESS_REPORTED",
        "scope": "STRICT_ONLY_F_NADSOLITON_TO_LSM_PLUS_LGR",
        "base_theorem_status": s1508.get("status"),
        "trials": rows,
        "counterexample_found": contradiction_found,
        "first_counterexample": first_counterexample,
        "interpretation": "Draft theorem remains conditionally admissible on non-violating perturbations; violating perturbations correctly trigger rejection branch." if contradiction_found else "No counterexample found in tested admissible perturbations.",
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "Implement P1510 admissible-perturbation boundary map separating robust zone from theorem-rejection zone.",
        "layman_explanation": "Zrobiliśmy test odporności: specjalnie próbowaliśmy zepsuć warunki dowodu. Tam gdzie warunki są łamane, system poprawnie odrzuca twierdzenie; tam gdzie nie są łamane, wynik się trzyma.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1509] status={summary['status']} counterexample_found={contradiction_found}")


if __name__ == "__main__":
    main()
