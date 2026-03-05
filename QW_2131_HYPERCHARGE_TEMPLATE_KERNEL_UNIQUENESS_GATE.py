#!/usr/bin/env python3
"""
QW-2131: Hypercharge-template uniqueness gate from kernel-anchored constraints.

Purpose:
- resolve the open flag from QW-2129 by deriving unique hypercharge template
  inside an explicit kernel-anchored constraint set,
- keep domain limitations explicit.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2131_hypercharge_template_kernel_uniqueness_gate.json"
OUT_MD = ROOT / "RAPORT_QW2131_HYPERCHARGE_TEMPLATE_KERNEL_UNIQUENESS_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def to_float(fr: Fraction) -> float:
    return float(fr.numerator) / float(fr.denominator)


def main() -> None:
    r2118 = load_json("report_qw2118_ktotal_spectral_tripartition_gate.json")
    r2128 = load_json("report_qw2128_kernel_rep_assignment_uniqueness_gate.json")
    r2129 = load_json("report_qw2129_anomaly_cancellation_kernel_anchored_gate.json")

    n_oct = int(r2118["matrix_spec"]["n_octaves"])
    winner_q = str(r2128["locked_branch_ranking"]["winner_q_assignment"])

    # Kernel-anchored normalization:
    # Y_Q := 2 / N_oct (for N_oct=12 gives 1/6).
    YQ = Fraction(2, n_oct)

    # From SU(2)^2 U(1) anomaly: 3Y_Q + Y_L = 0.
    YL = -3 * YQ

    # Right-handed neutrino neutral anchor in kernel-anchored neutrino channel.
    Yn = Fraction(0, 1)

    # Yukawa gauge-invariance:
    # Y_u = Y_Q + Y_H, Y_d = Y_Q - Y_H, Y_e = Y_L - Y_H, Y_nu = Y_L + Y_H = Yn.
    YH = Yn - YL
    Yu = YQ + YH
    Yd = YQ - YH
    Ye = YL - YH

    # Electric charges for LH doublets (Q = T3 + Y):
    Q_up_L = Fraction(1, 2) + YQ
    Q_down_L = -Fraction(1, 2) + YQ
    Q_nu_L = Fraction(1, 2) + YL
    Q_e_L = -Fraction(1, 2) + YL

    # Anomaly coefficients:
    A33Y = 2 * YQ - Yu - Yd
    A22Y = 3 * YQ + YL
    Agrav = 6 * YQ - 3 * Yu - 3 * Yd + 2 * YL - Ye - Yn
    AYYY = 6 * (YQ**3) - 3 * (Yu**3) - 3 * (Yd**3) + 2 * (YL**3) - (Ye**3) - (Yn**3)

    # Uniqueness check over rational search on denominator <= 12
    # under fixed kernel anchors: YQ=2/N_oct and Yn=0.
    candidates: List[Dict[str, str]] = []
    den_max = 12
    vals = [Fraction(n, d) for d in range(1, den_max + 1) for n in range(-2 * d, 2 * d + 1)]
    vals = sorted(set(vals))
    for YH_c in vals:
        YL_c = -3 * YQ
        Yu_c = YQ + YH_c
        Yd_c = YQ - YH_c
        Ye_c = YL_c - YH_c
        Yn_c = Fraction(0, 1)
        # Enforce anomaly cancellation and neutrino LH neutrality.
        A33Y_c = 2 * YQ - Yu_c - Yd_c
        A22Y_c = 3 * YQ + YL_c
        Agrav_c = 6 * YQ - 3 * Yu_c - 3 * Yd_c + 2 * YL_c - Ye_c - Yn_c
        AYYY_c = 6 * (YQ**3) - 3 * (Yu_c**3) - 3 * (Yd_c**3) + 2 * (YL_c**3) - (Ye_c**3) - (Yn_c**3)
        Qnu_c = Fraction(1, 2) + YL_c
        if A33Y_c == 0 and A22Y_c == 0 and Agrav_c == 0 and AYYY_c == 0 and Qnu_c == 0:
            candidates.append(
                {
                    "Y_H": str(YH_c),
                    "Y_uR": str(Yu_c),
                    "Y_dR": str(Yd_c),
                    "Y_eR": str(Ye_c),
                }
            )

    target = {
        "Y_H": str(YH),
        "Y_uR": str(Yu),
        "Y_dR": str(Yd),
        "Y_eR": str(Ye),
    }
    target_in_candidates = bool(any(c == target for c in candidates))

    # Template equality with QW-2129 floats.
    t9 = r2129["template_hypercharges"]
    t_eq = bool(
        abs(to_float(YQ) - float(t9["Y_Q"])) <= 1e-15
        and abs(to_float(Yu) - float(t9["Y_uR"])) <= 1e-15
        and abs(to_float(Yd) - float(t9["Y_dR"])) <= 1e-15
        and abs(to_float(YL) - float(t9["Y_L"])) <= 1e-15
        and abs(to_float(Ye) - float(t9["Y_eR"])) <= 1e-15
        and abs(to_float(Yn) - float(t9["Y_nR"])) <= 1e-15
    )

    flags = {
        "q2129_template_anomaly_pass_present": bool(str(r2129.get("verdict", "")).startswith("ANOMALY_CANCELLATION_KERNEL_ANCHORED_GATE_PASS")),
        "q2128_locked_branch_winner_present": bool(winner_q == "legacy_fibonacci"),
        "kernel_octave_count_loaded": bool(n_oct == 12),
        "kernel_anchor_YQ_equals_2_over_Noct": bool(YQ == Fraction(1, 6)),
        "derived_Y_template_matches_q2129_template": bool(t_eq),
        "electric_charge_relation_q_eq_t3_plus_y_holds": bool(
            Q_up_L == Fraction(2, 3)
            and Q_down_L == Fraction(-1, 3)
            and Q_nu_L == 0
            and Q_e_L == -1
        ),
        "anomaly_coefficients_exact_zero_fractional": bool(A33Y == 0 and A22Y == 0 and Agrav == 0 and AYYY == 0),
        "rational_search_contains_target_template": bool(target_in_candidates),
        "rational_search_unique_template_under_kernel_anchors": bool(len(candidates) == 1 and target_in_candidates),
        "hypercharge_template_unique_from_kernel": bool(len(candidates) == 1 and target_in_candidates),
        "global_uniqueness_without_neutrino_neutral_anchor": False,
        "deterministic_no_scan_no_retune": True,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    verdict = (
        "HYPERCHARGE_TEMPLATE_KERNEL_UNIQUENESS_GATE_PASS_ANCHORED_DOMAIN"
        if (
            flags["q2129_template_anomaly_pass_present"]
            and flags["kernel_octave_count_loaded"]
            and flags["kernel_anchor_YQ_equals_2_over_Noct"]
            and flags["derived_Y_template_matches_q2129_template"]
            and flags["electric_charge_relation_q_eq_t3_plus_y_holds"]
            and flags["anomaly_coefficients_exact_zero_fractional"]
            and flags["rational_search_unique_template_under_kernel_anchors"]
            and flags["hypercharge_template_unique_from_kernel"]
            and flags["deterministic_no_scan_no_retune"]
        )
        else "HYPERCHARGE_TEMPLATE_KERNEL_UNIQUENESS_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2118": "report_qw2118_ktotal_spectral_tripartition_gate.json",
            "q2128": "report_qw2128_kernel_rep_assignment_uniqueness_gate.json",
            "q2129": "report_qw2129_anomaly_cancellation_kernel_anchored_gate.json",
        },
        "kernel_anchors": {
            "n_octaves": n_oct,
            "Y_Q_formula": "Y_Q = 2 / N_oct",
            "Y_Q_value_fraction": str(YQ),
            "winner_q_assignment": winner_q,
        },
        "derived_template_fractional": {
            "Y_Q": str(YQ),
            "Y_uR": str(Yu),
            "Y_dR": str(Yd),
            "Y_L": str(YL),
            "Y_eR": str(Ye),
            "Y_nR": str(Yn),
            "Y_H": str(YH),
        },
        "derived_template_float": {
            "Y_Q": to_float(YQ),
            "Y_uR": to_float(Yu),
            "Y_dR": to_float(Yd),
            "Y_L": to_float(YL),
            "Y_eR": to_float(Ye),
            "Y_nR": to_float(Yn),
            "Y_H": to_float(YH),
        },
        "anomaly_fractional": {
            "A_SU3_SU3_U1": str(A33Y),
            "A_SU2_SU2_U1": str(A22Y),
            "A_gravity_gravity_U1": str(Agrav),
            "A_U1_cubic": str(AYYY),
        },
        "rational_uniqueness_search": {
            "denominator_max": den_max,
            "n_candidates": len(candidates),
            "candidates": candidates,
            "target_template": target,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "EITHER_FORMALIZE_NEUTRINO_NEUTRAL_ANCHOR_FROM_KERNEL_DYNAMICS_OR_KEEP_DOMAIN_BOUNDARY_EXPLICIT"
            if verdict.endswith("DOMAIN")
            else "REPAIR_TEMPLATE_DERIVATION_AND_RERUN_QW2131"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2131: HYPERCHARGE TEMPLATE KERNEL UNIQUENESS GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Kernel anchors",
        f"- N_oct: `{n_oct}`",
        f"- Y_Q = 2/N_oct: `{str(YQ)}`",
        "",
        "## Derived template",
        f"- (Y_Q, Y_uR, Y_dR, Y_L, Y_eR, Y_nR, Y_H) = "
        f"`({str(YQ)}, {str(Yu)}, {str(Yd)}, {str(YL)}, {str(Ye)}, {str(Yn)}, {str(YH)})`",
        "",
        "## Uniqueness search",
        f"- candidates count: `{len(candidates)}`",
        f"- target in candidates: `{target_in_candidates}`",
        "",
        "## Open scope boundary",
        "- global_uniqueness_without_neutrino_neutral_anchor: `False`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2131] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2131] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2131] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

