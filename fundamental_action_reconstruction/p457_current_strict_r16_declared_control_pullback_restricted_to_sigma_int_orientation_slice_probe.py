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

IN_R15 = GENERATED / "r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"
IN_R16 = GENERATED / "r16_explicit_residual_local_diagonal_declared_control_pullback_packet_for_host_matching_route.json"
IN_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
IN_QW2190 = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"

OUT_RESTRICTED = GENERATED / "m_control_residual_diag_declared_restricted_to_sigma_int_orientation_slice_v1.json"
OUT_JSON = GENERATED / "p457_current_strict_r16_declared_control_pullback_restricted_to_sigma_int_orientation_slice_probe.json"
OUT_SUMMARY = (
    GENERATED / "p457_current_strict_r16_declared_control_pullback_restricted_to_sigma_int_orientation_slice_probe_summary.json"
)


def load_json_path(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def clean_scalar(value: float) -> float:
    if abs(value) < 1e-15:
        return 0.0
    return round(float(value), 15)


def build_entry_expression(coefficients: list[float], residual_terms: list[str]) -> str:
    pieces: list[str] = []
    for coeff, residual in zip(coefficients, residual_terms):
        cleaned = clean_scalar(float(coeff))
        if cleaned == 0.0:
            continue
        if cleaned == 1.0:
            pieces.append(f"({residual})")
        elif cleaned == -1.0:
            pieces.append(f"-({residual})")
        else:
            pieces.append(f"({cleaned})*({residual})")
    if not pieces:
        return "0.0"
    expression = pieces[0]
    for piece in pieces[1:]:
        if piece.startswith("-"):
            expression += " - " + piece[1:]
        else:
            expression += " + " + piece
    return expression


def real_fourier_basis(n: int) -> dict[str, np.ndarray]:
    x = np.arange(n, dtype=float)
    basis: dict[str, np.ndarray] = {"e0": np.ones(n, dtype=float) / math.sqrt(n)}
    for m in range(1, n // 2):
        basis[f"c{m}"] = math.sqrt(2.0 / n) * np.cos(2.0 * math.pi * m * x / n)
        basis[f"s{m}"] = math.sqrt(2.0 / n) * np.sin(2.0 * math.pi * m * x / n)
    if n % 2 == 0:
        basis[f"e{n//2}"] = ((-1.0) ** x) / math.sqrt(n)
    return basis


def wrap_angle_mod_2pi(theta: float) -> float:
    t = float(theta) % (2.0 * math.pi)
    if t < 0.0:
        t += 2.0 * math.pi
    return t


def min_angle_diff_mod_2pi(a: float, b: float) -> float:
    da = wrap_angle_mod_2pi(a) - wrap_angle_mod_2pi(b)
    da = (da + math.pi) % (2.0 * math.pi) - math.pi
    return float(abs(da))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing: list[str] = []
    for p in (IN_R15, IN_R16, IN_THETA_PAIR, IN_QW2190):
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
    if missing:
        raise SystemExit(json.dumps({"status": "MISSING_INPUT_ARTIFACTS", "missing": missing}, ensure_ascii=True))

    r15 = load_json_path(IN_R15)
    r16 = load_json_path(IN_R16)
    theta_pair = load_json_path(IN_THETA_PAIR)
    q2190 = load_json_path(IN_QW2190)

    n = int(q2190["mode_mapping"]["n_octaves"])
    if n != 12:
        raise SystemExit(json.dumps({"status": "UNEXPECTED_N", "expected": 12, "actual": n}, ensure_ascii=True))

    expected_theta_pair_obj = "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1"
    if theta_pair.get("object") != expected_theta_pair_obj:
        raise SystemExit(
            json.dumps(
                {
                    "status": "UNEXPECTED_THETA_PAIR_OBJECT",
                    "expected": expected_theta_pair_obj,
                    "actual": theta_pair.get("object"),
                },
                ensure_ascii=True,
            )
        )

    residual_terms = [str(row["residual_local_diagonal_entry"]) for row in r15["diagonal_decomposition"]["entrywise_rows"]]
    if len(residual_terms) != n:
        raise SystemExit(
            json.dumps(
                {
                    "status": "UNEXPECTED_RESIDUAL_TERMS_LENGTH",
                    "expected": n,
                    "actual": len(residual_terms),
                },
                ensure_ascii=True,
            )
        )

    control_basis = r16["declared_control_pullback_packet"]["control_basis"]
    expected_control_basis = ["c1", "s1", "c2", "s2"]
    if control_basis != expected_control_basis:
        raise SystemExit(
            json.dumps(
                {
                    "status": "UNEXPECTED_CONTROL_BASIS",
                    "expected": expected_control_basis,
                    "actual": control_basis,
                },
                ensure_ascii=True,
            )
        )

    coeff_tensor = np.array(r16["declared_control_pullback_packet"]["coefficient_tensor_by_carrier_slot"], dtype=float)
    if coeff_tensor.shape != (4, 4, n):
        raise SystemExit(
            json.dumps(
                {
                    "status": "UNEXPECTED_COEFFICIENT_TENSOR_SHAPE",
                    "expected": [4, 4, n],
                    "actual": list(coeff_tensor.shape),
                },
                ensure_ascii=True,
            )
        )

    theta_1 = float(theta_pair["outputs"]["pair1"]["theta_1"])
    theta_2 = float(theta_pair["outputs"]["pair2"]["theta_2"])
    u1 = np.array([float(x) for x in theta_pair["outputs"]["pair1"]["u_1"]], dtype=float)
    u2 = np.array([float(x) for x in theta_pair["outputs"]["pair2"]["u_2"]], dtype=float)

    fourier = real_fourier_basis(n)
    c1, s1, c2, s2 = fourier["c1"], fourier["s1"], fourier["c2"], fourier["s2"]

    x1 = float(np.dot(c1, u1))
    y1 = float(np.dot(s1, u1))
    u1_outside = float(np.linalg.norm(u1 - (x1 * c1 + y1 * s1)))

    x2 = float(np.dot(c2, u2))
    y2 = float(np.dot(s2, u2))
    u2_outside = float(np.linalg.norm(u2 - (x2 * c2 + y2 * s2)))

    theta_1_recovered = float(math.atan2(y1, x1))
    theta_2_recovered = float(math.atan2(y2, x2))

    # Control coordinates in basis (c1,s1,c2,s2)
    u1_ctrl = np.array([x1, y1, 0.0, 0.0], dtype=float)
    u2_ctrl = np.array([0.0, 0.0, x2, y2], dtype=float)
    U = np.column_stack([u1_ctrl, u2_ctrl])  # 4x2

    gram = U.T @ U
    gram_res = float(np.linalg.norm(gram - np.eye(2)))

    # Build restricted coefficient-vectors by carrier slot: for each slot i, R_i = U^T C_i U
    restricted_by_slot = np.zeros((n, 2, 2), dtype=float)
    for i in range(n):
        restricted_by_slot[i] = U.T @ coeff_tensor[:, :, i] @ U

    coeff_00 = [clean_scalar(float(restricted_by_slot[i, 0, 0])) for i in range(n)]
    coeff_01 = [clean_scalar(float(restricted_by_slot[i, 0, 1])) for i in range(n)]
    coeff_10 = [clean_scalar(float(restricted_by_slot[i, 1, 0])) for i in range(n)]
    coeff_11 = [clean_scalar(float(restricted_by_slot[i, 1, 1])) for i in range(n)]

    sym_res = float(max(abs(float(a - b)) for a, b in zip(coeff_01, coeff_10)))

    expr_00 = build_entry_expression(coeff_00, residual_terms)
    expr_01 = build_entry_expression(coeff_01, residual_terms)
    expr_10 = build_entry_expression(coeff_10, residual_terms)
    expr_11 = build_entry_expression(coeff_11, residual_terms)

    tol = 1e-12
    checks = [
        {
            "id": "u1_in_span_c1_s1",
            "pass": bool(u1_outside <= tol),
            "expected": f"<= {tol:g}",
            "actual": u1_outside,
            "meaning": "sigma-int u_1 lies in span{c1,s1} (numerical tolerance)",
        },
        {
            "id": "u2_in_span_c2_s2",
            "pass": bool(u2_outside <= tol),
            "expected": f"<= {tol:g}",
            "actual": u2_outside,
            "meaning": "sigma-int u_2 lies in span{c2,s2} (numerical tolerance)",
        },
        {
            "id": "recovered_theta_1_matches_exported_mod_2pi",
            "pass": bool(min_angle_diff_mod_2pi(theta_1_recovered, theta_1) <= tol),
            "expected": f"<= {tol:g}",
            "actual": min_angle_diff_mod_2pi(theta_1_recovered, theta_1),
            "meaning": "atan2(<s1,u1>,<c1,u1>) matches exported theta_1 modulo 2π",
        },
        {
            "id": "recovered_theta_2_matches_exported_mod_2pi",
            "pass": bool(min_angle_diff_mod_2pi(theta_2_recovered, theta_2) <= tol),
            "expected": f"<= {tol:g}",
            "actual": min_angle_diff_mod_2pi(theta_2_recovered, theta_2),
            "meaning": "atan2(<s2,u2>,<c2,u2>) matches exported theta_2 modulo 2π",
        },
        {
            "id": "slice_basis_is_orthonormal_in_control_coords",
            "pass": bool(gram_res <= tol),
            "expected": f"<= {tol:g}",
            "actual": gram_res,
            "meaning": "u1_ctrl,u2_ctrl are orthonormal in the control coordinate basis",
        },
        {
            "id": "restricted_coefficients_symmetric",
            "pass": bool(sym_res <= tol),
            "expected": f"<= {tol:g}",
            "actual": sym_res,
            "meaning": "restricted coefficient decomposition is symmetric: coeff_01 == coeff_10",
        },
    ]

    passed = bool(all(bool(c["pass"]) for c in checks))

    restricted = {
        "object": "M_control_residual_diag_declared_restricted_to_sigma_int_orientation_slice_v1",
        "status": "actual_exported_restriction_artifact__declared_residual_diagonal_sector_only__no_false_pass",
        "as_of": "2026-03-15",
        "intent": (
            "Restrict the declared residual local diagonal control pullback (R16: M_control_residual_diag_declared = T_control^T D_local_residual T_control) "
            "to the sigma-int orientation slice basis (u_1,u_2) exported by the slot-free sigma-int theta-pair source (F451/N489). "
            "This yields an explicit 2x2 linear coefficient decomposition over the residual local diagonal entries from R15. "
            "No vanishing/cancellation witness is claimed."
        ),
        "inputs": {
            "r15_residual_local_diagonal_entries": str(IN_R15.relative_to(REPO)),
            "r16_declared_control_pullback": str(IN_R16.relative_to(REPO)),
            "sigma_int_theta_pair": str(IN_THETA_PAIR.relative_to(REPO)),
        },
        "control_basis_order": expected_control_basis,
        "slice_basis": {
            "theta_1": theta_1,
            "theta_2": theta_2,
            "u1_ctrl": [float(x) for x in u1_ctrl.tolist()],
            "u2_ctrl": [float(x) for x in u2_ctrl.tolist()],
        },
        "restricted_coefficients_by_carrier_slot": {
            "entry_00": coeff_00,
            "entry_01": coeff_01,
            "entry_10": coeff_10,
            "entry_11": coeff_11,
        },
        "restricted_matrix_rows": [[expr_00, expr_01], [expr_10, expr_11]],
        "hard_limits": [
            "Restriction is only for the declared residual local diagonal sector (D_local_residual), not for the full Psi-sector Hessian block.",
            "Does not claim the restricted matrix vanishes.",
            "Does not claim host-side cancellation or matching.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    OUT_RESTRICTED.write_text(json.dumps(restricted, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P457",
        "goal": "restrict_R16_declared_residual_diagonal_control_pullback_to_sigma_int_orientation_slice_basis",
        "inputs": {
            "r15": str(IN_R15.relative_to(REPO)),
            "r16": str(IN_R16.relative_to(REPO)),
            "theta_pair": str(IN_THETA_PAIR.relative_to(REPO)),
            "qw2190_scaffold": str(IN_QW2190.relative_to(REPO)),
        },
        "constructed_artifact": str(OUT_RESTRICTED.relative_to(REPO)),
        "audits": {
            "u1_outside_span_c1_s1": u1_outside,
            "u2_outside_span_c2_s2": u2_outside,
            "theta_1_recovered": theta_1_recovered,
            "theta_2_recovered": theta_2_recovered,
            "theta_1_exported": theta_1,
            "theta_2_exported": theta_2,
            "theta_1_mod_2pi_abs_diff": min_angle_diff_mod_2pi(theta_1_recovered, theta_1),
            "theta_2_mod_2pi_abs_diff": min_angle_diff_mod_2pi(theta_2_recovered, theta_2),
            "slice_gram_residual": gram_res,
            "restricted_symmetry_coeff_max_abs_diff": sym_res,
            "tolerance": tol,
        },
        "checks": checks,
        "verdict": {
            "audits_pass": passed,
            "statement": (
                "Restriction artifact is constructible: R16 declared residual-diagonal control pullback can be restricted to the sigma-int orientation slice "
                "basis (u_1,u_2) exported by the slot-free sigma-int theta-pair supply. This does not claim any cancellation, host matching, or global closure."
            ),
        },
        "hard_limits": restricted["hard_limits"],
        "no_false_pass": True,
    }

    summary = {
        "stage": report["stage"],
        "status": "PASS_RESTRICTION_ARTIFACT_EXPORTED" if passed else "FAIL_RESTRICTION_ARTIFACT_EXPORTED",
        "audits_pass": passed,
        "constructed_artifact": report["constructed_artifact"],
        "theta_1": theta_1,
        "theta_2": theta_2,
        "slice_gram_residual": gram_res,
        "restricted_symmetry_coeff_max_abs_diff": sym_res,
        "tolerance": tol,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()

