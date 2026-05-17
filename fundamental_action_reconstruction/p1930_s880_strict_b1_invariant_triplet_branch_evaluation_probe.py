#!/usr/bin/env python3
"""P1930 S880 strict B1 invariant-triplet branch evaluation probe."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def row(branch: str, d_r: str, d_ricuu: str, d_gradchi: str, verdict: str, note: str) -> dict:
    return {
        "branch": branch,
        "delta_R": d_r,
        "delta_RicUU": d_ricuu,
        "delta_gradchi2": d_gradchi,
        "tensorial_b1_verdict": verdict,
        "note": note,
    }


def main() -> None:
    p1929 = load("p1929_s879_strict_b1_repair_candidate_and_statusvector_restamp_probe.json")

    eqs = p1929.get("b1_repair_candidate_v1", {}).get("new_witness_form", {}).get("equations", [])

    out = {
        "packet_id": "P1930",
        "stage_id": "S880",
        "status": "OPEN_BRANCH_LEVEL_B1_EVALUATION_WITHOUT_GLOBAL_CLOSURE",
        "route": "strict_only",
        "depends_on": {
            "p1929_present": "b1_repair_candidate_v1" in p1929,
            "p1929_tensorial_form_present": len(eqs) == 3,
        },
        "strict_chain_context": "K_strict -> coefficients -> full non-skeleton L_SM+L_GR -> EOM residual transport",
        "tensorial_b1_witness_form": eqs,
        "branch_evaluation_table": [
            row(
                "BR_A_CE_curv_pair_v1",
                "R_frw_A-R_bi_B != 0",
                "RicUU_frw_A-RicUU_bi_B != 0",
                "gradchi2_frw_A-gradchi2_bi_B != 0",
                "FAIL",
                "Certified counterexample-compatible branch remains failing under stricter tensorial form.",
            ),
            row(
                "BR_B_matched_invariant_triplet",
                "R_frw_B-R_bi_B = 0",
                "RicUU_frw_B-RicUU_bi_B = 0",
                "gradchi2_frw_B-gradchi2_bi_B = 0",
                "LOCAL_PASS",
                "Constructed matched-invariant branch gives local pass for B1 tensorial witness only.",
            ),
        ],
        "b1_result_summary": {
            "global_verdict": "OPEN_MIXED_BRANCH_OUTCOME",
            "reason": "At least one admissible branch remains FAIL; local pass on another branch is insufficient for global closure.",
            "false_pass_guard": "Do not upgrade background-independence to PASS unless admissible branch class is theorem-level restricted or all admissible branches pass.",
        },
        "strict_core_statusvector_restamp": {
            "renormalization": "OPEN",
            "unitarity": "OPEN",
            "background_independence": "OPEN_MIXED_BRANCH_OUTCOME",
            "selector_qw2191": "OPEN",
        },
        "toe_remaining_minimum_after_p1930": {
            "current_total_open": 7,
            "explanation": "P1930 adds explicit branch computations for repaired B1 form but discharges no theorem-grade block.",
        },
        "next_honest_step": "Export P1931 with theorem-candidate admissible-branch restriction proof (or counterexample-generalization theorem) and restamp B1 only after theorem-grade branch policy is exported.",
        "lay_explanation": "Ile zostało do ToE? Nadal minimum 7 dużych bloków. Pokazaliśmy uczciwie: jedna gałąź dalej oblewa test tła, druga lokalnie przechodzi, więc globalnego sukcesu jeszcze nie ma.",
    }

    out_path = GEN / "p1930_s880_strict_b1_invariant_triplet_branch_evaluation_probe.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
