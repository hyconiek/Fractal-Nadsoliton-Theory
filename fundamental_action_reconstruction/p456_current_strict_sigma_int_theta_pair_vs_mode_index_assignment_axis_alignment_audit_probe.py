#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
IN_ASSIGNMENT_DIAGONAL = GENERATED / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json"
IN_ASSIGNMENT_SHANNON = GENERATED / "mode_index_assignment_shannon_element_order_reference_strict_core_v1.json"

OUT_JSON = GENERATED / "p456_current_strict_sigma_int_theta_pair_vs_mode_index_assignment_axis_alignment_audit_probe.json"
OUT_SUMMARY = GENERATED / "p456_current_strict_sigma_int_theta_pair_vs_mode_index_assignment_axis_alignment_audit_probe_summary.json"


def load_json_path(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def as_vec(xs: list[Any]) -> np.ndarray:
    return np.array([float(x) for x in xs], dtype=float)


def abs_dot(a: np.ndarray, b: np.ndarray) -> float:
    return float(abs(float(np.dot(a, b))))


def norm_residual(v: np.ndarray) -> float:
    return float(abs(float(np.linalg.norm(v)) - 1.0))


def projector(columns: list[np.ndarray]) -> np.ndarray:
    u = np.column_stack(columns)
    return u @ u.T


def axis_alignment_record(u: np.ndarray, u_plus: np.ndarray, u_minus: np.ndarray) -> dict[str, Any]:
    d_plus = abs_dot(u, u_plus)
    d_minus = abs_dot(u, u_minus)
    best = "u_plus" if d_plus >= d_minus else "u_minus"
    best_abs_dot = float(max(d_plus, d_minus))
    cross_abs_dot = float(min(d_plus, d_minus))
    return {
        "abs_dot_u_plus": float(d_plus),
        "abs_dot_u_minus": float(d_minus),
        "best_match": best,
        "best_abs_dot": best_abs_dot,
        "cross_abs_dot": cross_abs_dot,
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing: list[str] = []
    for p in (IN_THETA_PAIR, IN_ASSIGNMENT_DIAGONAL, IN_ASSIGNMENT_SHANNON):
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
    if missing:
        raise SystemExit(json.dumps({"status": "MISSING_INPUT_ARTIFACTS", "missing": missing}, ensure_ascii=True))

    theta_pair = load_json_path(IN_THETA_PAIR)
    diag = load_json_path(IN_ASSIGNMENT_DIAGONAL)
    shannon = load_json_path(IN_ASSIGNMENT_SHANNON)

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

    u1 = as_vec(theta_pair["outputs"]["pair1"]["u_1"])
    u2 = as_vec(theta_pair["outputs"]["pair2"]["u_2"])

    u1_norm_res = norm_residual(u1)
    u2_norm_res = norm_residual(u2)
    u1_dot_u2 = float(np.dot(u1, u2))

    diag_p1 = diag["outputs"]["pairs"]["pair1"]
    diag_p2 = diag["outputs"]["pairs"]["pair2"]
    sha_p1 = shannon["outputs"]["pairs"]["pair1"]
    sha_p2 = shannon["outputs"]["pairs"]["pair2"]

    records = {
        "diagonal_local": {
            "pair1": axis_alignment_record(u1, as_vec(diag_p1["u_plus"]), as_vec(diag_p1["u_minus"])),
            "pair2": axis_alignment_record(u2, as_vec(diag_p2["u_plus"]), as_vec(diag_p2["u_minus"])),
        },
        "shannon_ord_reference": {
            "pair1": axis_alignment_record(u1, as_vec(sha_p1["u_plus"]), as_vec(sha_p1["u_minus"])),
            "pair2": axis_alignment_record(u2, as_vec(sha_p2["u_plus"]), as_vec(sha_p2["u_minus"])),
        },
    }

    tol = 1e-12

    max_projector_sign_dev = 0.0
    base_proj = projector([u1, u2])
    for s1 in (-1.0, 1.0):
        for s2 in (-1.0, 1.0):
            alt_proj = projector([s1 * u1, s2 * u2])
            max_projector_sign_dev = max(max_projector_sign_dev, float(np.linalg.norm(alt_proj - base_proj)))

    checks: list[dict[str, Any]] = [
        {
            "id": "sigma_int_theta_pair_vectors_unit_norm",
            "pass": bool(u1_norm_res <= tol and u2_norm_res <= tol),
            "expected": f"norm residual <= {tol:g}",
            "actual": {"u1_norm_residual": u1_norm_res, "u2_norm_residual": u2_norm_res},
            "meaning": "u_1, u_2 exported by the sigma-int theta-pair artifact are unit vectors (numerical tolerance)",
        },
        {
            "id": "sigma_int_theta_pair_vectors_orthogonal",
            "pass": bool(abs(u1_dot_u2) <= tol),
            "expected": f"|dot| <= {tol:g}",
            "actual": {"u1_dot_u2": u1_dot_u2},
            "meaning": "u_1 ⟂ u_2 as required by the R1 target-slot span semantics (numerical tolerance)",
        },
        {
            "id": "r1_span_projector_invariant_under_sign_flips",
            "pass": bool(max_projector_sign_dev <= tol),
            "expected": f"<= {tol:g}",
            "actual": {"max_projector_sign_dev": max_projector_sign_dev},
            "meaning": "the R1 target-slot object class depends only on span{u_1,u_2}, hence its projector is invariant under residual Z2 sign flips",
        },
    ]

    for lane_id, lane in records.items():
        for pair_id, rec in lane.items():
            checks.append(
                {
                    "id": f"axis_alignment_{lane_id}_{pair_id}",
                    "pass": bool(rec["best_abs_dot"] >= 1.0 - tol and rec["cross_abs_dot"] <= tol),
                    "expected": {"best_abs_dot": f">= {1.0 - tol:g}", "cross_abs_dot": f"<= {tol:g}"},
                    "actual": rec,
                    "meaning": "sigma-int theta-pair axis aligns with one of the two exported lane basis axes up to residual Z2 sign",
                }
            )

    passed = bool(all(bool(c["pass"]) for c in checks))
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P456",
        "goal": (
            "audit_axis_alignment_between_the_exported_slot_free_sigma_int_theta_pair_source_(pair1/pair2) "
            "and_the_exported_mode_index_assignment_basis_objects_(diagonal_local_and_shannon_ord_reference), "
            "and_confirm_span_semantics_sign_invariance_for_R1"
        ),
        "inputs": {
            "theta_pair_sigma_int": str(IN_THETA_PAIR.relative_to(REPO)),
            "mode_index_assignment_diagonal_local": str(IN_ASSIGNMENT_DIAGONAL.relative_to(REPO)),
            "mode_index_assignment_shannon_ord_reference": str(IN_ASSIGNMENT_SHANNON.relative_to(REPO)),
        },
        "records": records,
        "audits": {
            "u1_norm_residual": u1_norm_res,
            "u2_norm_residual": u2_norm_res,
            "u1_dot_u2": u1_dot_u2,
            "max_projector_sign_dev": max_projector_sign_dev,
            "tolerance": tol,
        },
        "checks": checks,
        "verdict": {
            "audits_pass": passed,
            "statement": (
                "On the current exported strict objects, the sigma-int slot-free theta-pair axes (pair1/pair2) align "
                "with the exported mode-index assignment axes on both the diagonal/local lane and the Shannon ord-reference lane "
                "up to residual Z2 sign. In addition, the R1 target-slot span projector is invariant under sign flips. "
                "This is an internal consistency audit and does not promote any global closure claim."
            ),
        },
        "hard_limits": [
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191 (kernel-alone obstruction remains true).",
            "Does not claim sign-sensitive physical orientation datum (only span semantics and axis alignment up to sign).",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": report["stage"],
        "status": "PASS_AXIS_ALIGNMENT_UP_TO_RESIDUAL_Z2" if passed else "FAIL_AXIS_ALIGNMENT_UP_TO_RESIDUAL_Z2",
        "audits_pass": passed,
        "records": records,
        "u1_norm_residual": u1_norm_res,
        "u2_norm_residual": u2_norm_res,
        "u1_dot_u2": u1_dot_u2,
        "max_projector_sign_dev": max_projector_sign_dev,
        "tolerance": tol,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()

