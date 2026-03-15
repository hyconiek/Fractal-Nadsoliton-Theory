#!/usr/bin/env python3
from __future__ import annotations

import argparse
import cmath
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P437_JSON = (
    GENERATED
    / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe.json"
)
P437_SUMMARY = (
    GENERATED
    / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe_summary.json"
)

OUT_JSON = (
    GENERATED
    / "p449_current_strict_canonical_local_diagonal_multi_pair_o2_cut_defect_evaluation_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p449_current_strict_canonical_local_diagonal_multi_pair_o2_cut_defect_evaluation_probe_summary.json"
)


def is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "P449: evaluate multi-pair diagonal/local Fourier defects F_{2m}(d) for m=1..5 from the P437-computed "
            "canonical diagonal residual profile d_k (n=12)."
        )
    )
    p.add_argument(
        "--in-json",
        dest="in_json",
        default=str(P437_JSON),
        help="Input JSON path (default: P437 output artifact containing d_local_residual_profile).",
    )
    p.add_argument("--out-json", dest="out_json", default=str(OUT_JSON), help="Output artifact JSON path.")
    p.add_argument("--out-summary", dest="out_summary", default=str(OUT_SUMMARY), help="Output summary JSON path.")
    p.add_argument(
        "--tolerance",
        dest="tolerance",
        type=float,
        default=1e-12,
        help="Numeric nonzero tolerance for |F_{2m}| (default: 1e-12).",
    )
    return p.parse_args()


def fourier_coeff(d: list[float], l: int) -> complex:
    n = len(d)
    s = 0j
    for k, x in enumerate(d):
        s += float(x) * cmath.exp(1j * 2.0 * math.pi * l * k / n)
    return s


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    args = parse_args()
    in_json = Path(args.in_json)
    out_json = Path(args.out_json)
    out_summary = Path(args.out_summary)
    tol = float(args.tolerance)

    missing_files: list[str] = []
    if not in_json.exists():
        missing_files.append(str(in_json))
    if not P437_SUMMARY.exists():
        missing_files.append(str(P437_SUMMARY.relative_to(REPO)))

    if missing_files:
        summary = {
            "stage": "P449",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FILES",
            "missing_files": missing_files,
            "theorem_level_pass": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }
        out_summary.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(out_summary)
        return

    p437 = load_json(in_json)
    p437_summary = load_json(P437_SUMMARY)

    d = (((p437.get("computed") or {}) if isinstance(p437.get("computed"), dict) else {}).get("d_local_residual_profile"))
    if not (isinstance(d, list) and len(d) == 12 and all(is_number(x) for x in d)):
        summary = {
            "stage": "P449",
            "status": "NOT_COMPUTABLE_MISSING_INPUT_VALUES",
            "missing_or_invalid": ["P437.computed.d_local_residual_profile (length-12 finite numeric list)"],
            "theorem_level_pass": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }
        out_summary.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(out_summary)
        return

    d_v = [float(x) for x in d]
    n = 12
    mean_d = float(sum(d_v) / n)

    pairs: dict[str, Any] = {}
    all_cut = True
    for m in range(1, 6):
        l = 2 * m
        F = fourier_coeff(d_v, l)
        absF = float(abs(F))
        cuts = absF > tol
        all_cut = all_cut and cuts
        theta_star = 0.5 * math.atan2(float(F.imag), float(F.real)) if cuts else None
        pairs[f"pair{m}"] = {
            "m": m,
            "l": l,
            "F": {"Re": float(F.real), "Im": float(F.imag), "abs": absF},
            "cuts_O2_by_absF_nonzero": bool(cuts),
            "theta_star": (float(theta_star) if theta_star is not None else None),
            "eigenvalues": {"lambda_plus": mean_d + absF / n, "lambda_minus": mean_d - absF / n} if cuts else None,
        }

    theorem_level_pass = bool(p437_summary.get("theorem_level_pass"))

    artifact = {
        "stage": "P449",
        "goal": "evaluate_multi_pair_fourier_defects_F_2m_for_canonical_local_diagonal_residual_profile_from_P437",
        "inputs": {
            "p437_output_artifact": str(in_json),
            "p437_summary": str(P437_SUMMARY.relative_to(REPO)),
            "tolerance": tol,
        },
        "computed": {
            "n": n,
            "mean_d": mean_d,
            "pairs": pairs,
            "all_pairs_cut": bool(all_cut),
        },
        "provenance_note": "This probe evaluates Fourier defects from the P437-computed diagonal residual profile. It does not itself supply new strict inputs.",
        "no_false_pass": True,
    }

    summary = {
        "stage": "P449",
        "status": "PASS_EVALUATED",
        "all_pairs_cut": bool(all_cut),
        "pair_cut_flags": {k: bool(v["cuts_O2_by_absF_nonzero"]) for k, v in pairs.items()},
        "F_abs": {k: float(v["F"]["abs"]) for k, v in pairs.items()},
        "input_theorem_level_pass_from_P437": bool(theorem_level_pass),
        "theorem_level_pass": bool(theorem_level_pass),
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    out_json.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    out_summary.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out_json)


if __name__ == "__main__":
    main()

