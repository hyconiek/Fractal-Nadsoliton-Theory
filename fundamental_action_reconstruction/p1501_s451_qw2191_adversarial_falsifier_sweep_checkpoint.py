#!/usr/bin/env python3
"""P1501 S4.51: adversarial falsifier sweep for QW-2191 R2/R5 witness objects."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1500 = GEN / "p1500_s450_qw2191_selector_source_and_f_mapping_witness_summary.json"
P1491 = GEN / "p1491_s441_qw2191_kappa_robustness_sweep_summary.json"

SUMMARY = GEN / "p1501_s451_qw2191_adversarial_falsifier_sweep_summary.json"


def main() -> None:
    s1500 = json.loads(P1500.read_text(encoding="utf-8"))
    s1491 = json.loads(P1491.read_text(encoding="utf-8"))

    sm_w = float(s1500["objects"]["W_Fmap_v1"]["F_to_LSM_weight"])
    gr_w = float(s1500["objects"]["W_Fmap_v1"]["F_to_LGR_weight"])
    s_src = float(s1500["objects"]["S_internal_v1"]["value"])
    margin = float(s1491["safety_margin"])

    adversarial_grid = [-1.0, -0.5, 0.0, 0.5, 1.0]
    rows = []
    for a in adversarial_grid:
        # bounded perturbation in policy-scale
        d = a * margin
        sm_t = sm_w + d
        gr_t = gr_w - d
        gap_before = abs(sm_w - gr_w)
        gap_after = abs((sm_t - gr_t) - s_src * gap_before)

        cond_positive = sm_t > 0 and gr_t > 0
        cond_orientation = (sm_t - gr_t) > 0
        cond_gap = gap_after < gap_before

        fail = not (cond_positive and cond_orientation and cond_gap)
        rows.append({
            "adv": a,
            "delta": d,
            "sm_t": sm_t,
            "gr_t": gr_t,
            "gap_before": gap_before,
            "gap_after": gap_after,
            "cond_positive": cond_positive,
            "cond_orientation": cond_orientation,
            "cond_gap": cond_gap,
            "status": "FAIL" if fail else "PASS",
        })

    falsifiers = [r for r in rows if r["status"] == "FAIL"]

    summary = {
        "packet": "P1501",
        "status": "PASS_ADVERSARIAL_FALSIFIER_SWEEP_LOCAL_ONLY" if not falsifiers else "FAIL_ADVERSARIAL_FALSIFIER_SWEEP_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "rows": rows,
        "falsifier_count": len(falsifiers),
        "falsifiers": falsifiers,
        "global_closure_candidate_strengthened": len(falsifiers) == 0,
        "qw2191_closed": False,
        "next_step_recommendation": "S4.52: publish global-closure candidate note with explicit caveat boundaries and request independent external replication.",
        "layman_explanation": "To test 'na złość': specjalnie szukamy ustawień, które mogłyby obalić mechanizm. Jeśli nie znajdujemy obalenia w dopuszczalnym zakresie, rośnie wiarygodność rozwiązania.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1501] status={summary['status']} falsifier_count={len(falsifiers)}")


if __name__ == "__main__":
    main()
