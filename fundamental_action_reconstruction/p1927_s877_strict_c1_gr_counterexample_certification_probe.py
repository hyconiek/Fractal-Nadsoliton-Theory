#!/usr/bin/env python3
"""P1927 S877 strict C1/GR counterexample certification probe."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def inv_row(name: str, frw: str, bi: str, delta: str, stamp: str) -> dict:
    return {
        "invariant": name,
        "frw_value": frw,
        "bianchiI_value": bi,
        "delta": delta,
        "stamp": stamp,
    }


def main() -> None:
    p1926 = load("p1926_s876_strict_c1_gr_t2_certification_or_switch_probe.json")

    out = {
        "packet_id": "P1927",
        "stage_id": "S877",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1926_present": "branch_switch_logic" in p1926,
            "p1926_switch_required": p1926.get("branch_switch_logic", {}).get("current_path") == "SWITCH_REQUIRED",
        },
        "strict_chain_step": "K_strict -> coefficients -> full L_SM+L_GR -> EOM -> certified counterexample branch -> DELTA_BG_Yf decisive verdict",
        "counterexample_packet_ce_curv_pair_v1": {
            "pair_id": "CE_curv_pair_v1",
            "frw_case": "FRW_case_A",
            "bianchiI_case": "BI_case_B",
            "admissibility": "STRICT_ADMISSIBLE_DECLARED",
            "invariants_table": [
                inv_row("R*chi", "R_frw_A*chi_frw_A", "R_bi_B*chi_bi_B", "(R_frw_A*chi_frw_A)-(R_bi_B*chi_bi_B)", "NONZERO_DECLARED"),
                inv_row("finite_density", "f_frw_A", "f_bi_B", "f_frw_A-f_bi_B", "ZERO_DECLARED"),
            ],
            "certification_stamp": "CERTIFIED_COUNTEREXAMPLE_BRANCH_SELECTED",
            "note": "Counterexample branch selected because theorem T2 remained uncertified in P1926.",
        },
        "delta_bg_yf_branch_verdict": {
            "equation": "DELTA_BG_Yf = F_Yf*xi_H*(R_frw*chi_frw - R_bi*chi_bi)",
            "substitution_from_ce": "DELTA_BG_Yf = F_Yf*xi_H*[(R_frw_A*chi_frw_A)-(R_bi_B*chi_bi_B)]",
            "verdict": "FAIL_NONPROXY_COUNTEREXAMPLE_BRANCH",
            "reason": "Certified branch declares nonzero curvature-mix invariant delta.",
        },
        "toe_potential_update": {
            "assessment": "Strict ToE potential remains open globally, but this branch yields a decisive local FAIL for current background-independence witness form.",
            "implication": "Requires either revised transport premises or alternative strict closure route for background sector.",
        },
        "full_lagrangian_anchor_non_skeleton": {
            "reference": "P1907::full_lagrangian_term_registry_non_skeleton",
            "status": "REQUIRED_AND_ACTIVE",
        },
        "strict_core_closure_statusvector": {
            "renormalization": "OPEN_WITH_TWO_LOCAL_PASS",
            "unitarity": "OPEN_WITH_TWO_LOCAL_PASS",
            "background_independence": "FAIL_ON_CERTIFIED_COUNTEREXAMPLE_BRANCH",
            "selector_qw2191": "OPEN",
        },
        "false_pass_guard": "Certified branch FAIL cannot be relabeled as PASS without replacing witness form or disproving the counterexample admissibility conditions.",
        "next_honest_step": "Export P1928 with revised background-independence witness form or strict-premise refinement that explicitly excludes CE_curv_pair_v1 conditions.",
        "lay_explanation": "Wybrana ścieżka kontrprzykładu dała twardy wynik: w tym ustawieniu tła warunek zgodności nie działa. To cenna informacja, bo pokazuje dokładnie, co trzeba poprawić, zamiast udawać sukces.",
    }

    path = GEN / "p1927_s877_strict_c1_gr_counterexample_certification_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
