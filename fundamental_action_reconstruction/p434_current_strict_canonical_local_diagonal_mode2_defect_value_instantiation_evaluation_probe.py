#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_VALUES = GENERATED / "p434_input_sigma_opposite_pair_sum_values_candidate.json"

OUT_JSON = (
    GENERATED
    / "p434_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_evaluation_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p434_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_evaluation_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = [
        "Sigma_psi0_psi6",
        "Sigma_psi1_psi7",
        "Sigma_psi2_psi8",
        "Sigma_psi3_psi9",
        "Sigma_psi4_psi10",
        "Sigma_psi5_psi11",
    ]

    vals = load_json(IN_VALUES)
    missing = [k for k in required if not is_number(vals.get(k))]

    artifact: dict[str, Any] = {
        "stage": "P434",
        "goal": "evaluate_the_mode2_defect_F2_and_the_induced_pair1_O2_cut_and_axis_data_given_a_numeric_instantiation_of_the_six_opposite_pair_sums",
        "inputs": {
            "sigma_opposite_pair_sums_candidate_values": str(IN_VALUES.relative_to(REPO)),
        },
        "required_keys": required,
        "missing_or_non_numeric_keys": missing,
        "result": {
            "computable": False,
            "cuts_O2_on_pair1_by_N466": None,
            "F2": None,
            "theta_star_by_N468": None,
        },
        "no_false_pass": True,
    }

    if missing:
        summary = {
            "stage": "P434",
            "status": "NOT_COMPUTABLE_MISSING_INPUT_VALUES",
            "missing_or_non_numeric_keys": missing,
            "theorem_level_pass": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_JSON)
        return

    sigmas = [float(vals[k]) for k in required]

    re_coeffs = [1.0, 0.5, -0.5, -1.0, -0.5, 0.5]
    s3 = math.sqrt(3.0) / 2.0
    im_coeffs = [0.0, s3, s3, 0.0, -s3, -s3]

    re_f2 = sum(c * s for c, s in zip(re_coeffs, sigmas, strict=True))
    im_f2 = sum(c * s for c, s in zip(im_coeffs, sigmas, strict=True))
    abs_f2 = math.hypot(re_f2, im_f2)

    # N466: (a1-d1, b1) = ( (1/6)Re(F2), (1/12)Im(F2) ) for n=12.
    a1_minus_d1 = (1.0 / 6.0) * re_f2
    b1 = (1.0 / 12.0) * im_f2

    # Mean diagonal value on pair1:
    # sum_{i=0..11} d_i = sum_{k=0..5} (d_k + d_{k+6}) = sum sigmas.
    mean = (1.0 / 12.0) * sum(sigmas)

    # Pair1 block entries (optional explicit):
    a1 = mean + (1.0 / 12.0) * re_f2
    d1 = mean - (1.0 / 12.0) * re_f2

    # Pair1 eigenvalues:
    lambda_plus = mean + (1.0 / 12.0) * abs_f2
    lambda_minus = mean - (1.0 / 12.0) * abs_f2

    tol = 1e-12
    cuts = abs_f2 > tol

    theta_star = None
    if cuts:
        # N468: theta_* = (1/2) atan2(Im(F2), Re(F2)).
        theta_star = 0.5 * math.atan2(im_f2, re_f2)

    artifact["result"] = {
        "computable": True,
        "Sigma_values": {k: float(vals[k]) for k in required},
        "F2": {"Re": re_f2, "Im": im_f2, "abs": abs_f2},
        "pair1_signature_by_N466": {"a1_minus_d1": a1_minus_d1, "b1": b1},
        "pair1_block_entries": {"a1_c1c1": a1, "b1_c1s1": b1, "d1_s1s1": d1},
        "pair1_eigenvalues": {"lambda_plus": lambda_plus, "lambda_minus": lambda_minus},
        "cuts_O2_on_pair1_by_N466": cuts,
        "theta_star_by_N468": theta_star,
        "tolerance": tol,
    }

    summary = {
        "stage": "P434",
        "status": "PASS_EVALUATED",
        "cuts_O2_on_pair1_by_N466": cuts,
        "F2_abs": abs_f2,
        "theta_star_by_N468": theta_star,
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

