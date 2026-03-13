#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_R18 = (
    GENERATED
    / "r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route.json"
)

OUT_JSON = (
    GENERATED
    / "p433_current_strict_canonical_local_diagonal_mode2_defect_via_r18_pair1_entry_system_consistency_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p433_current_strict_canonical_local_diagonal_mode2_defect_via_r18_pair1_entry_system_consistency_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b, strict=True))


def compute_re_im_f2_from_sigma(sigmas: list[float]) -> tuple[float, float]:
    # N467/P426 coefficients on Sigma_psi{k}_psi{k+6} (k=0..5).
    re_coeffs = [1.0, 0.5, -0.5, -1.0, -0.5, 0.5]
    s3 = math.sqrt(3.0) / 2.0
    im_coeffs = [0.0, s3, s3, 0.0, -s3, -s3]
    return dot(re_coeffs, sigmas), dot(im_coeffs, sigmas)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r18 = load_json(IN_R18)
    classes = {row["class_symbol"]: row for row in r18["pair1_coefficient_class_reduction"]["coefficient_classes"]}

    class_order = [
        "Sigma_psi0_psi6",
        "Sigma_psi1_psi7",
        "Sigma_psi2_psi8",
        "Sigma_psi3_psi9",
        "Sigma_psi4_psi10",
        "Sigma_psi5_psi11",
    ]

    coeffs = {
        "c1c1": [float(classes[s]["signature_on_pair1_entries"]["c1c1"]) for s in class_order],
        "c1s1": [float(classes[s]["signature_on_pair1_entries"]["c1s1"]) for s in class_order],
        "s1s1": [float(classes[s]["signature_on_pair1_entries"]["s1s1"]) for s in class_order],
    }

    rng = random.Random(7331)
    num_trials = 24
    # R18 stores several transport coefficients rounded to ~12 decimal digits, so we audit at a tolerance that is
    # robust under that truncation when Sigma-values are O(1).
    tol = 1e-10

    trials: list[dict[str, Any]] = []
    max_abs_re_err = 0.0
    max_abs_im_err = 0.0

    for _ in range(num_trials):
        sigmas = [rng.uniform(-3.0, 3.0) for _ in range(6)]

        a1 = dot(coeffs["c1c1"], sigmas)
        b1 = dot(coeffs["c1s1"], sigmas)
        d1 = dot(coeffs["s1s1"], sigmas)

        re1, im1 = compute_re_im_f2_from_sigma(sigmas)
        re2 = 6.0 * (a1 - d1)
        im2 = 12.0 * b1

        re_err = re2 - re1
        im_err = im2 - im1
        max_abs_re_err = max(max_abs_re_err, abs(re_err))
        max_abs_im_err = max(max_abs_im_err, abs(im_err))

        trials.append(
            {
                "Sigma_values": sigmas,
                "pair1_entries_from_R18": {"a1_c1c1": a1, "b1_c1s1": b1, "d1_s1s1": d1},
                "Re_Im_F2_six_class_formula": {"Re_F2": re1, "Im_F2": im1},
                "Re_Im_F2_from_pair1_entries": {"Re_F2": re2, "Im_F2": im2},
                "errors": {"Re_F2_error": re_err, "Im_F2_error": im_err},
            }
        )

    artifact = {
        "stage": "P433",
        "goal": "audit_numeric_consistency_of_the_N473_identity_recovering_the_mode2_defect_F2_from_the_R18_pair1_entry_system",
        "inputs": {
            "r18_pair1_coefficient_class_reduction": str(IN_R18.relative_to(REPO)),
        },
        "coeffs_used": coeffs,
        "trial_count": num_trials,
        "tolerance": tol,
        "max_abs_errors": {"Re_F2": max_abs_re_err, "Im_F2": max_abs_im_err},
        "trials": trials,
        "result": {
            "identity_holds_within_tolerance": (max_abs_re_err <= tol and max_abs_im_err <= tol),
            "strict_decision_of_F2_for_canonical_profile_obtained": False,
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P433",
        "status": "PASS_IDENTITY_AUDITED"
        if artifact["result"]["identity_holds_within_tolerance"]
        else "FAIL_IDENTITY_MISMATCH",
        "trial_count": num_trials,
        "max_abs_errors": artifact["max_abs_errors"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
