#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_DIAGONAL_ASSIGNMENT = GENERATED / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json"
IN_SHANNON_ASSIGNMENT = GENERATED / "mode_index_assignment_shannon_element_order_reference_strict_core_v1.json"

OUT_JSON = GENERATED / "p455_current_strict_mode_index_assignment_shannon_vs_diagonal_alignment_audit_probe.json"
OUT_SUMMARY = GENERATED / "p455_current_strict_mode_index_assignment_shannon_vs_diagonal_alignment_audit_probe_summary.json"


def load_json_path(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def angle_mod_pi_distance(a: float, b: float) -> float:
    # Treat angles as elements of R / (pi Z) since the diagonalization convention is only meaningful mod pi.
    aa = float(a) % math.pi
    bb = float(b) % math.pi
    d = abs(aa - bb) % math.pi
    return float(min(d, math.pi - d))


def best_pair_basis_match(
    uplus_a: np.ndarray,
    uminus_a: np.ndarray,
    uplus_b: np.ndarray,
    uminus_b: np.ndarray,
) -> dict[str, Any]:
    ma = np.column_stack([uplus_a, uminus_a])
    mb = np.column_stack([uplus_b, uminus_b])
    m = ma.T @ mb

    # Two possibilities: aligned (u+->u+, u-->u-) or swapped (u+->u-, u-->u+), each up to independent Z2 signs.
    score_aligned = abs(float(m[0, 0])) + abs(float(m[1, 1]))
    score_swapped = abs(float(m[0, 1])) + abs(float(m[1, 0]))
    if score_swapped > score_aligned:
        mapping = "swap"
        primary = {"u_plus": {"target": "u_minus", "dot": float(m[0, 1])}, "u_minus": {"target": "u_plus", "dot": float(m[1, 0])}}
        cross = {"u_plus_to_u_plus": float(m[0, 0]), "u_minus_to_u_minus": float(m[1, 1])}
    else:
        mapping = "aligned"
        primary = {"u_plus": {"target": "u_plus", "dot": float(m[0, 0])}, "u_minus": {"target": "u_minus", "dot": float(m[1, 1])}}
        cross = {"u_plus_to_u_minus": float(m[0, 1]), "u_minus_to_u_plus": float(m[1, 0])}

    diag_abs = [abs(float(primary["u_plus"]["dot"])), abs(float(primary["u_minus"]["dot"]))]
    cross_abs = [abs(float(v)) for v in cross.values()]
    return {
        "mapping": mapping,
        "primary_matches": primary,
        "cross_dots": cross,
        "min_primary_abs_dot": float(min(diag_abs)),
        "max_cross_abs_dot": float(max(cross_abs)),
        "alignment_matrix_aT_b": [[float(m[0, 0]), float(m[0, 1])], [float(m[1, 0]), float(m[1, 1])]],
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing: list[str] = []
    for p in (IN_DIAGONAL_ASSIGNMENT, IN_SHANNON_ASSIGNMENT):
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
    if missing:
        summary = {
            "stage": "P455",
            "status": "NOT_COMPUTABLE_MISSING_ASSIGNMENT_OBJECTS",
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    diag = load_json_path(IN_DIAGONAL_ASSIGNMENT)
    sha = load_json_path(IN_SHANNON_ASSIGNMENT)

    n_diag = int(diag["outputs"]["n"])
    n_sha = int(sha["outputs"]["n"])
    if n_diag != n_sha:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "P455",
                    "status": "INCOMPATIBLE_N",
                    "n_diagonal": n_diag,
                    "n_shannon": n_sha,
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    n = n_diag
    if n != 12:
        raise SystemExit(json.dumps({"stage": "P455", "status": "UNSUPPORTED_N_EXPECTED_12", "n": n, "no_false_pass": True}, ensure_ascii=True))

    tol_theta = 1e-12
    tol_dot = 1e-12

    e0_diag = np.array([float(x) for x in diag["outputs"]["e0"]], dtype=float)
    e0_sha = np.array([float(x) for x in sha["outputs"]["e0"]], dtype=float)

    e6_diag = np.array([float(x) for x in diag["outputs"]["e6"]], dtype=float)
    e6_sha = np.array([float(x) for x in sha["outputs"]["e6"]], dtype=float)

    audits: dict[str, Any] = {
        "e0_dot": float(np.dot(e0_diag, e0_sha)),
        "e6_dot": float(np.dot(e6_diag, e6_sha)),
        "e0_l2_diff": float(np.linalg.norm(e0_diag - e0_sha)),
        "e6_l2_diff": float(np.linalg.norm(e6_diag - e6_sha)),
    }

    pairs: dict[str, Any] = {}
    all_pairs_aligned = True
    max_theta_star_mod_pi_diff = 0.0
    min_primary_abs_dot = 1.0
    max_cross_abs_dot = 0.0

    for m in range(1, n // 2):
        pd = diag["outputs"]["pairs"][f"pair{m}"]
        ps = sha["outputs"]["pairs"][f"pair{m}"]

        td = float(pd["theta_star"])
        ts = float(ps["theta_star"])
        theta_diff = angle_mod_pi_distance(td, ts)
        max_theta_star_mod_pi_diff = max(max_theta_star_mod_pi_diff, float(theta_diff))

        uplus_d = np.array([float(x) for x in pd["u_plus"]], dtype=float)
        uminus_d = np.array([float(x) for x in pd["u_minus"]], dtype=float)
        uplus_s = np.array([float(x) for x in ps["u_plus"]], dtype=float)
        uminus_s = np.array([float(x) for x in ps["u_minus"]], dtype=float)

        match = best_pair_basis_match(uplus_d, uminus_d, uplus_s, uminus_s)

        min_primary_abs_dot = min(min_primary_abs_dot, float(match["min_primary_abs_dot"]))
        max_cross_abs_dot = max(max_cross_abs_dot, float(match["max_cross_abs_dot"]))

        aligned = bool(theta_diff <= tol_theta and match["min_primary_abs_dot"] >= (1.0 - tol_dot) and match["max_cross_abs_dot"] <= tol_dot)
        all_pairs_aligned = all_pairs_aligned and aligned

        pairs[f"pair{m}"] = {
            "m": m,
            "theta_star": {"diagonal": td, "shannon": ts, "mod_pi_abs_diff": float(theta_diff)},
            "basis_match": match,
            "aligned_up_to_residual_z2": aligned,
        }

    status = "PASS_ALIGNED_UP_TO_RESIDUAL_Z2" if all_pairs_aligned else "FAIL_ALIGNMENT_MISMATCH"

    artifact = {
        "stage": "P455",
        "goal": (
            "audit_alignment_of_two_independently_exported_strict_mode_index_assignments_on_n12_fourier_scaffold_"
            "diagonal_local_strict_derived_vs_shannon_element_order_reference"
        ),
        "status": status,
        "inputs": {
            "diagonal_local_assignment": str(IN_DIAGONAL_ASSIGNMENT.relative_to(REPO)),
            "shannon_assignment": str(IN_SHANNON_ASSIGNMENT.relative_to(REPO)),
        },
        "tolerances": {"theta_mod_pi_abs_diff": tol_theta, "basis_dot_residual": tol_dot},
        "audits": {
            **audits,
            "max_theta_star_mod_pi_abs_diff": float(max_theta_star_mod_pi_diff),
            "min_primary_abs_dot_across_pairs": float(min_primary_abs_dot),
            "max_cross_abs_dot_across_pairs": float(max_cross_abs_dot),
        },
        "pairs": pairs,
        "meaning_no_false_pass": [
            "This is an alignment audit on the current exported n=12 instantiations only.",
            "It does not claim theorem-level identification of the diagonal/local residual profile with ord_Z12(x).",
            "It does not claim strict-core selector closure nor global discharge of QW-2191.",
            "Residual Z2 sign remains a convention unless separately proven gauge-irrelevant for a given downstream observable.",
        ],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "no_false_pass": True,
    }

    summary = {
        "stage": "P455",
        "status": status,
        "aligned_all_pairs_up_to_residual_z2": bool(all_pairs_aligned),
        "max_theta_star_mod_pi_abs_diff": float(max_theta_star_mod_pi_diff),
        "min_primary_abs_dot_across_pairs": float(min_primary_abs_dot),
        "max_cross_abs_dot_across_pairs": float(max_cross_abs_dot),
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

